/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 22 Apr 2024 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "constitutive_laws/hypoelastic_finite_strain_bridging_constitutive_law_dc.h"
#include "structural_application_variables.h"

// #define DEBUG_CONSTITUTIVE_LAW
#define DEBUG_ELEMENT_ID 1

#ifdef DEBUG_CONSTITUTIVE_LAW
#include <iomanip>
#endif

namespace Kratos
{

//**********************************************************************
template<int THWSchemeType, int TStressType>
bool HypoelasticFiniteStrainBridgingConstitutiveLawDC<THWSchemeType, TStressType>::Has( const Variable<Matrix>& rThisVariable )
{
    if (rThisVariable == CURRENT_DEFORMATION_GRADIENT)
        return true;
    return BaseType::mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLawDC<THWSchemeType, TStressType>::SetValue( const Variable<double>& rThisVariable, const double& rValue,
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
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLawDC<THWSchemeType, TStressType>::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
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
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLawDC<THWSchemeType, TStressType>::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    // compute the equivalent (small) strain vector
    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetStrainSize(dim);

    Matrix DDu_half(3, 3); // TODO moving this calculation to element
    BaseType::ComputeDDuMidpoint( DDu_half, geom, BaseType::mIntPoint );

    Vector StrainVector(strain_size);
    BaseType::ComputeStrain(StrainVector, DDu_half);

    // integrate the (small strain) constitutive law, obtaining Cauchy stress
    BaseType::mpConstitutiveLaw->SetValue(CURRENT_STRAIN_VECTOR, StrainVector, CurrentProcessInfo);
    BaseType::mpConstitutiveLaw->FinalizeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
    if (!BaseType::mpConstitutiveLaw->Has(CAUCHY_STRESS_TENSOR))
        KRATOS_ERROR << "Constitutive law is not able to return CAUCHY_STRESS_TENSOR";
    BaseType::mpConstitutiveLaw->GetValue(CAUCHY_STRESS_TENSOR, BaseType::m_stress_n1);

    if constexpr (TStressType == 2)
    {
        // transform the stress to Cauchy stress
        BaseType::m_stress_n1 /= BaseType::m_J_n1;
    }
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLawDC<THWSchemeType, TStressType>::CalculateMaterialResponseCauchy (Parameters& rValues)
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
        SuperType::ComputeTangent(rValues, AlgorithmicTangent);
    }
    if (rValues.IsSetStressVector())
    {
        Vector& StressVector = rValues.GetStressVector();
        SD_MathUtils<double>::StressTensorToVector(BaseType::m_stress_n1, StressVector);
    }
}

//**********************************************************************
//**********************************************************************

template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC<THWSchemeType, TStressType>::CalculateDu( const unsigned int dim,
        Matrix& DDu, const Matrix& G_Operator, const Matrix& CurrentDisp ) const
{
    SD_MathUtils<double>::CalculateFaxi<true>( DDu, G_Operator, CurrentDisp );
}

template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC<THWSchemeType, TStressType>::CalculateG( Matrix& G_Operator,
        const Vector& N, const Matrix& DN_DX, const GeometryType& rGeometry ) const
{
    const unsigned int number_of_nodes = N.size();

    if (G_Operator.size1() != 5 || G_Operator.size2() != 2*number_of_nodes)
        G_Operator.resize( 5, 2*number_of_nodes, false );

    SD_MathUtils<double>::CalculateGaxi( G_Operator, rGeometry, N, DN_DX );
}

//**********************************************************************

template class HypoelasticFiniteStrainBridgingConstitutiveLawDC<1, 1>;
template class HypoelasticFiniteStrainBridgingConstitutiveLawDC<1, 2>;
template class HypoelasticFiniteStrainBridgingConstitutiveLawDC<2, 1>;
template class HypoelasticFiniteStrainBridgingConstitutiveLawDC<2, 2>;

template class HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC<1, 1>;
template class HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC<1, 2>;
template class HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC<2, 1>;
template class HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC<2, 2>;

} // Namespace Kratos

#undef DEBUG_CONSTITUTIVE_LAW
#undef DEBUG_ELEMENT_ID
