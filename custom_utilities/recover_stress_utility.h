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
#include "includes/element.h"
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

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    KRATOS_CLASS_POINTER_DEFINITION(RecoverStressUtility);

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

};

}  // namespace Kratos.

#endif // KRATOS_RECOVER_STRESS_UTILITY_H_INCLUDED  defined
