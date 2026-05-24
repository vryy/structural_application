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
#include "containers/interface_container.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"


namespace Kratos
{

/*
 * Utility for error estimation and stress recovery. It used LOCAL_ERROR exclusively for the error estimation
 * at the element level
 */
class RecoverStressUtility
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(RecoverStressUtility);

    typedef Element::GeometryType GeometryType;

    typedef typename GeometryType::IntegrationMethod IntegrationMethod;

    typedef typename GeometryType::LocalCoordinatesArrayType LocalCoordinatesArrayType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename GeometryType::JacobiansType JacobiansType;

    typedef typename GeometryType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    ///@}
    ///@name Operations
    ///@{

    /// Reset the elemental values, especially the local error so that it can be computed
    // by the estimator
    template<typename TContainerType>
    static void ResetLocalError(TContainerType& rElements)
    {
        for (auto it = rElements.begin(); it != rElements.end(); ++it)
            it->SetValue(LOCAL_ERROR, 0.0);
    }

    /// Compute the error estimation at the element as sqrt( \int_Omega sigma : Ce^-1 : sigma )
    template<typename TEntityType>
    static typename TEntityType::DataType ComputeZZErrorEstimation(TEntityType& rElement, const ProcessInfo& rCurrentProcessInfo);

    /// Compute the Kelly error estimation across all the elements in the model_part. On output,
    /// returns the sum of them
    template<int TDim, typename TVariableType>
    static double ComputeKellyErrorEstimation(ModelPart& rModelPart, const TVariableType& rVariable)
    {
        ResetLocalError(rModelPart.Elements());
        InterfaceContainer<ModelPart> icon;
        icon.template ConstructInterfaces<TDim>(rModelPart);
        double result = ComputeKellyErrorEstimation(icon, rVariable);
        return result;
    }

    /// Compute the Kelly error estimation across all the half-faces. On output,
    /// returns the sum of them
    template<typename TVariableType>
    static double ComputeKellyErrorEstimation(const InterfaceContainer<ModelPart>& icon, const TVariableType& rVariable);

    ///@}

private:

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
