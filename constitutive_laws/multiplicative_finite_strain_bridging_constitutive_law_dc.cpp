/*
LICENSE: see material_point_application/LICENSE.txt
*/
/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 20 Sep 2023 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/eigen_utility.h"
#include "constitutive_laws/multiplicative_finite_strain_bridging_constitutive_law_dc.h"
#include "structural_application_variables.h"

// #define DEBUG_CONSTITUTIVE_LAW
#define DEBUG_ELEMENT_ID 1

#ifdef DEBUG_CONSTITUTIVE_LAW
#include <iomanip>
#endif

namespace Kratos
{

//**********************************************************************
template<int TStressType>
MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>::MultiplicativeFiniteStrainBridgingConstitutiveLawDC()
    : BaseType()
{
}

//**********************************************************************
template<int TStressType>
MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>::MultiplicativeFiniteStrainBridgingConstitutiveLawDC(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : BaseType(pConstitutiveLaw)
{
}

//**********************************************************************
template<int TStressType>
MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>::~MultiplicativeFiniteStrainBridgingConstitutiveLawDC()
{
}

//**********************************************************************
template<int TStressType>
bool MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>::Has( const Variable<Matrix>& rThisVariable )
{
    if (rThisVariable == CURRENT_DEFORMATION_GRADIENT)
        return true;
    return BaseType::mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
template<int TStressType>
Matrix& MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    if (rThisVariable == THREED_ALGORITHMIC_TANGENT)
    {
        if (rValue.size1() != 9 || rValue.size2() != 9)
            rValue.resize(9, 9, false);
        this->ComputeTangent( rValue );
        return rValue;
    }
    else if (rThisVariable == CAUCHY_STRESS_TENSOR)
    {
        if (rValue.size1() != 3 || rValue.size2() != 3)
            rValue.resize(3, 3, false);
        noalias(rValue) = BaseType::m_stress_n1;
        return rValue;
    }

    return BaseType::mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == CURRENT_DEFORMATION_GRADIENT_DETERMINANT)
    {
        BaseType::m_J_n1 = rValue;
        return;
    }

    BaseType::mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == CURRENT_DEFORMATION_GRADIENT)
    {
        const Matrix& F = rValue;
        if (F.size1() == 3)
        {
            noalias(BaseType::m_F_n1) = F;
        }
        else if (F.size1() == 2)
        {
            BaseType::m_F_n1.clear();
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    BaseType::m_F_n1(i, j) = F(i, j);
            BaseType::m_F_n1(2, 2) = 1.0;
        }
        BaseType::m_J_n1 = MathUtils<double>::Det(BaseType::m_F_n1);
        return;
    }

    BaseType::mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
template<>
void MultiplicativeFiniteStrainBridgingConstitutiveLawDC<1>::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    // compute the equivalent (small) strain vector
    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetStrainSize(dim);

    Vector StrainVector(strain_size);
    BaseType::ComputeLogarithmicStrain(StrainVector, BaseType::m_left_elastic_cauchy_green_tensor_trial, BaseType::m_F_n1);

    // integrate the (small strain) constitutive law, obtaining Cauchy stress
    BaseType::mpConstitutiveLaw->SetValue(CURRENT_STRAIN_VECTOR, StrainVector, CurrentProcessInfo);
    BaseType::mpConstitutiveLaw->FinalizeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
    if (!BaseType::mpConstitutiveLaw->Has(CAUCHY_STRESS_TENSOR))
        KRATOS_ERROR << "Constitutive law is not able to return CAUCHY_STRESS_TENSOR";
    BaseType::mpConstitutiveLaw->GetValue(CAUCHY_STRESS_TENSOR, BaseType::m_stress_n1);
}

//**********************************************************************
template<>
void MultiplicativeFiniteStrainBridgingConstitutiveLawDC<2>::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    #ifdef DEBUG_CONSTITUTIVE_LAW
    int ElemId, GaussId;
    BaseType::mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
    BaseType::mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, GaussId);
    std::cout << std::setprecision(16);
    #endif

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == DEBUG_ELEMENT_ID && GaussId == 0)
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(BaseType::m_F_n1)
    }
    #endif

    // compute the equivalent (small) strain vector
    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetStrainSize(dim);

    Vector StrainVector(strain_size);
    BaseType::ComputeLogarithmicStrain(StrainVector, BaseType::m_left_elastic_cauchy_green_tensor_trial, BaseType::m_F_n1);

    // integrate the (small strain) constitutive law, obtaining Kirchhoff stress
    BaseType::mpConstitutiveLaw->SetValue(CURRENT_STRAIN_VECTOR, StrainVector, CurrentProcessInfo);
    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        // KRATOS_WATCH(StrainVector)
    }
    #endif
    BaseType::mpConstitutiveLaw->FinalizeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
    if (!BaseType::mpConstitutiveLaw->Has(CAUCHY_STRESS_TENSOR))
        KRATOS_ERROR << "Constitutive law is not able to return CAUCHY_STRESS_TENSOR";
    BaseType::mpConstitutiveLaw->GetValue(CAUCHY_STRESS_TENSOR, BaseType::m_stress_n1);
    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(BaseType::m_J_n1)
    }
    #endif
    // transform the stress to Cauchy stress
    BaseType::m_stress_n1 /= BaseType::m_J_n1;
    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == DEBUG_ELEMENT_ID && GaussId == 0)
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(BaseType::m_stress_n1)
        KRATOS_WATCH("-----------")
    }
    #endif
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
    if (CurrentProcessInfo[SET_CALCULATE_REACTION])
    {
        if (rValues.IsSetStressVector())
        {
            Vector& StressVector = rValues.GetStressVector();
            SD_MathUtils<double>::StressTensorToVector(BaseType::m_stress_n1, StressVector);
            return;
        }
    }

    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        this->ComputeTangent(AlgorithmicTangent);
    }
    if (rValues.IsSetStressVector())
    {
        Vector& StressVector = rValues.GetStressVector();
        SD_MathUtils<double>::StressTensorToVector(BaseType::m_stress_n1, StressVector);
    }
}

//**********************************************************************

template class MultiplicativeFiniteStrainBridgingConstitutiveLawDC<1>;
template class MultiplicativeFiniteStrainBridgingConstitutiveLawDC<2>;

} // Namespace Kratos

#undef DEBUG_CONSTITUTIVE_LAW
#undef DEBUG_ELEMENT_ID
