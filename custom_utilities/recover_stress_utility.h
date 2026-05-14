//
//   Project Name:        KratosStructuralApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 12 May 2026 $
//
//


#if !defined(KRATOS_RECOVER_STRESS_UTILITY_H_INCLUDED)
#define  KRATOS_RECOVER_STRESS_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cmath>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"


namespace Kratos
{

/*
 * Utility for operation with stress recovery
 */
class RecoverStressUtility
{
public:

    /**
     * Type Definitions
     */

    KRATOS_CLASS_POINTER_DEFINITION(RecoverStressUtility);

    typedef Element::GeometryType GeometryType;

    typedef typename GeometryType::IntegrationMethod IntegrationMethod;

    typedef typename GeometryType::LocalCoordinatesArrayType LocalCoordinatesArrayType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename GeometryType::JacobiansType JacobiansType;

    typedef typename GeometryType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    struct HalfFace
    {
        GeometryType face;                      // the representative geometry of the face, used for searching, indexing
        GeometryType::Pointer left  = nullptr;  // the real geometry interface (relatively left)
        GeometryType::Pointer right = nullptr;  // the real geometry interface (relatively right)
        Element::Pointer first  = nullptr;      // the element on the left side
        Element::Pointer second = nullptr;      // the element on the right side
        mutable double jump = 0.0;
    };

    struct Comparator
    {
        bool operator()(const HalfFace& a, const HalfFace& b) const
        {
            return a.face.IsLess(b.face);
        }
    };

    /**
     * Operations
     */

    template<typename TEntityType>
    static typename TEntityType::DataType ComputeZZErrorEstimation(TEntityType& rElement, const ProcessInfo& rCurrentProcessInfo)
    {
        typedef typename TEntityType::DataType DataType;

        typedef typename TEntityType::ValueType ValueType;

        typedef typename TEntityType::VectorType VectorType;

        typedef typename TEntityType::MatrixType MatrixType;

        typedef typename TEntityType::GeometryType GeometryType;

        GeometryType& rGeometry = rElement.GetGeometry();

        const IntegrationMethod ThisIntegrationMethod = rElement.GetIntegrationMethod();

        #ifdef ENABLE_BEZIER_GEOMETRY
        rGeometry.Initialize(ThisIntegrationMethod);
        #endif

        const typename GeometryType::IntegrationPointsArrayType& integration_points =
                    rGeometry.IntegrationPoints( ThisIntegrationMethod );

        std::vector<DataType> E;
        rElement.CalculateOnIntegrationPoints( VARSEL(DataType, YOUNG_MODULUS), E, rCurrentProcessInfo );

        std::vector<DataType> NU;
        rElement.CalculateOnIntegrationPoints( VARSEL(DataType, POISSON_RATIO), NU, rCurrentProcessInfo );

        std::vector<DataType> J0;
        rElement.CalculateOnIntegrationPoints( VARSEL(DataType, JACOBIAN_0), J0, rCurrentProcessInfo );

        std::vector<VectorType> stress;
        rElement.CalculateOnIntegrationPoints( VARSEL(DataType, STRESSES), stress, rCurrentProcessInfo );

        std::vector<VectorType> rstress;
        rElement.CalculateOnIntegrationPoints( VARSEL(DataType, RECOVERY_STRESSES), rstress, rCurrentProcessInfo );

        DataType result = 0.0;
        MatrixType stress_tensor(3, 3), rstress_tensor(3, 3), ds(3, 3), Eds(3, 3);
        typename SD_MathUtils<DataType>::Fourth_Order_Tensor Ci;

        for(std::size_t PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            ValueType weight = integration_points[PointNumber].Weight();

            SD_MathUtils<DataType>::StressVectorToTensor(stress[PointNumber], stress_tensor);
            SD_MathUtils<DataType>::StressVectorToTensor(rstress[PointNumber], rstress_tensor);

            SD_MathUtils<DataType>::CalculateInversedElasticTensor(Ci, E[PointNumber], NU[PointNumber]);

            noalias(ds) = (rstress_tensor - stress_tensor);
            Eds.clear();
            SD_MathUtils<DataType>::ContractFourthOrderTensor(1.0, Ci, ds, Eds);

            result += SD_MathUtils<DataType>::mat_inner_prod(ds, Eds) * weight * J0[PointNumber];
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        rGeometry.Clean();
        #endif

        return std::sqrt(result);
    }

    /// Compute the Kelly error estimation across all the elements in the model_part. On output,
    /// returns the sum of them
    template<typename TVariableType>
    static double ComputeKellyErrorEstimation(ModelPart& rModelPart, const TVariableType& rVariable);

private:

    static std::set<HalfFace, Comparator> ConstructHalfFaceStructure(const ElementsContainerType& rElements);

    static Vector ComputeGradient(const GeometryType& rGeometry, const LocalCoordinatesArrayType& rPoint, const Variable<double>& rVariable);

    static Matrix ComputeGradient(const GeometryType& rGeometry, const LocalCoordinatesArrayType& rPoint, const Variable<array_1d<double, 3> >& rVariable);

    static double ComputeJump(const GeometryType& rGeometry1, const GeometryType& rGeometry2,
            const LocalCoordinatesArrayType& p1, const LocalCoordinatesArrayType& p2,
            const Vector& n, const Variable<double>& rVariable);

    static double ComputeJump(const GeometryType& rGeometry1, const GeometryType& rGeometry2,
            const LocalCoordinatesArrayType& p1, const LocalCoordinatesArrayType& p2,
            const Vector& n, const Variable<array_1d<double, 3> >& rVariable);

};

}  // namespace Kratos.

#endif // KRATOS_RECOVER_STRESS_UTILITY_H_INCLUDED  defined
