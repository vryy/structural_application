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

#ifdef DEBUG_CONSTITUTIVE_LAW
#include <iomanip>
#endif

namespace Kratos
{

//**********************************************************************
MultiplicativeFiniteStrainBridgingConstitutiveLawDC::MultiplicativeFiniteStrainBridgingConstitutiveLawDC()
    : MultiplicativeFiniteStrainBridgingConstitutiveLaw()
{
}

//**********************************************************************
MultiplicativeFiniteStrainBridgingConstitutiveLawDC::MultiplicativeFiniteStrainBridgingConstitutiveLawDC(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : MultiplicativeFiniteStrainBridgingConstitutiveLaw(pConstitutiveLaw)
{
}

//**********************************************************************
MultiplicativeFiniteStrainBridgingConstitutiveLawDC::~MultiplicativeFiniteStrainBridgingConstitutiveLawDC()
{
}

//**********************************************************************
bool MultiplicativeFiniteStrainBridgingConstitutiveLawDC::Has( const Variable<Matrix>& rThisVariable )
{
    if (rThisVariable == CURRENT_DEFORMATION_GRADIENT)
        return true;
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
Matrix& MultiplicativeFiniteStrainBridgingConstitutiveLawDC::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    if (rThisVariable == THREED_ALGORITHMIC_TANGENT)
    {
        if (rValue.size1() != 9 || rValue.size2() != 9)
            rValue.resize(9, 9, false);
        ComputeTangent( rValue );
        return rValue;
    }
    else if (rThisVariable == CAUCHY_STRESS_TENSOR)
    {
        if (rValue.size1() != 3 || rValue.size2() != 3)
            rValue.resize(3, 3, false);
        noalias(rValue) = m_stress_n1;
        return rValue;
    }

    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLawDC::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == CURRENT_DEFORMATION_GRADIENT_DETERMINANT)
    {
        mCurrentDetF = rValue;
        return;
    }

    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLawDC::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == CURRENT_DEFORMATION_GRADIENT)
    {
        const Matrix& F = rValue;
        if (F.size1() == 3)
        {
            noalias(mCurrentF) = F;
        }
        else if (F.size1() == 2)
        {
            mCurrentF.clear();
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    mCurrentF(i, j) = F(i, j);
            mCurrentF(2, 2) = 1.0;
        }
        mCurrentDetF = MathUtils<double>::Det(mCurrentF);
        return;
    }

    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLawDC::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = dim*(dim+1) / 2;

    #ifdef DEBUG_CONSTITUTIVE_LAW
    int ElemId, GaussId;
    mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
    mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, GaussId);
    std::cout << std::setprecision(16);
    #endif

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == 1 && GaussId == 0)
    if (ElemId == 1)
    {
        KRATOS_WATCH(mCurrentF)
    }
    #endif

    /*
     * Integration algorithm for isotropic multiplicative finite strain elastoplasticity
     * Reference: Souza de Neto, Computational Plasticity, Box 14.3
     */

    // compute elastic left Cauchy-Green tensor from previous step
    // It is noted that the constitutive law must be able to provide ELASTIC_STRAIN_TENSOR
    Matrix elastic_strain_tensor(3, 3), elastic_strain_tensor_trial(3, 3), left_cauchy_green_tensor_n(3, 3);

    mpConstitutiveLaw->GetValue(ELASTIC_STRAIN_TENSOR, elastic_strain_tensor);

    EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2, left_cauchy_green_tensor_n, elastic_strain_tensor);

    // compute trial elastic left Cauchy-Green tensor
    Matrix Fincr(3, 3), invFn(3, 3);
    double detFn;
    MathUtils<double>::InvertMatrix(mLastF, invFn, detFn);
    noalias(Fincr) = prod(mCurrentF, invFn);

    noalias(m_left_cauchy_green_tensor_trial) = prod(Fincr, Matrix(prod(left_cauchy_green_tensor_n, trans(Fincr))));

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == 1 && GaussId == 0)
    if (ElemId == 1)
    {
        KRATOS_WATCH(elastic_strain_tensor)
        KRATOS_WATCH(left_cauchy_green_tensor_n)
        KRATOS_WATCH(Fincr)
        KRATOS_WATCH(mCurrentDetF)
        KRATOS_WATCH(m_left_cauchy_green_tensor_trial)
    }
    #endif

    // compute trial elastic logarithmic strain tensor
    std::vector<double> left_cauchy_green_tensor_trial_pri(3);
    std::vector<Matrix> left_cauchy_green_tensor_trial_eigprj(3);
    EigenUtility::calculate_principle_stresses(
            m_left_cauchy_green_tensor_trial(0, 0),
            m_left_cauchy_green_tensor_trial(1, 1),
            m_left_cauchy_green_tensor_trial(2, 2),
            m_left_cauchy_green_tensor_trial(0, 1),
            m_left_cauchy_green_tensor_trial(1, 2),
            m_left_cauchy_green_tensor_trial(0, 2),
            left_cauchy_green_tensor_trial_pri[0],
            left_cauchy_green_tensor_trial_pri[1],
            left_cauchy_green_tensor_trial_pri[2],
            left_cauchy_green_tensor_trial_eigprj[0],
            left_cauchy_green_tensor_trial_eigprj[1],
            left_cauchy_green_tensor_trial_eigprj[2]);
    SD_MathUtils<double>::ComputeIsotropicTensorFunction(EigenUtility::logd2, elastic_strain_tensor_trial, left_cauchy_green_tensor_trial_pri, left_cauchy_green_tensor_trial_eigprj);

    // create the strain vector as input to the small strain constitutive law
    Vector StrainVector(strain_size), IncrementalStrainVector(strain_size);
    SD_MathUtils<double>::StrainTensorToVector(elastic_strain_tensor_trial - elastic_strain_tensor, IncrementalStrainVector);
    mpConstitutiveLaw->GetValue(STRAIN, StrainVector);
    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        KRATOS_WATCH(StrainVector)
        KRATOS_WATCH(IncrementalStrainVector)
    }
    #endif
    noalias(StrainVector) += IncrementalStrainVector;

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == 1 && GaussId == 0)
    if (ElemId == 1)
    {
        KRATOS_WATCH(StrainVector)
    }
    #endif

    // integrate the (small strain) constitutive law, obtaining Kirchhoff stress
    mpConstitutiveLaw->SetValue(CURRENT_STRAIN_VECTOR, StrainVector, CurrentProcessInfo);
    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        // KRATOS_WATCH(StrainVector)
    }
    #endif
    mpConstitutiveLaw->FinalizeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
    // transform the stress to Cauchy stress
    mpConstitutiveLaw->GetValue(CAUCHY_STRESS_TENSOR, m_stress_n1);
    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        KRATOS_WATCH(mCurrentDetF)
    }
    #endif
    m_stress_n1 /= mCurrentDetF;
    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == 1 && GaussId == 0)
    if (ElemId == 1)
    {
        KRATOS_WATCH(m_stress_n1)
        KRATOS_WATCH("-----------")
    }
    #endif
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLawDC::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
    if (CurrentProcessInfo[SET_CALCULATE_REACTION])
    {
        if (rValues.IsSetStressVector())
        {
            Vector& StressVector = rValues.GetStressVector();
            SD_MathUtils<double>::StressTensorToVector(m_stress_n1, StressVector);
            return;
        }
    }

    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        ComputeTangent(AlgorithmicTangent);
    }
    if (rValues.IsSetStressVector())
    {
        Vector& StressVector = rValues.GetStressVector();
        SD_MathUtils<double>::StressTensorToVector(m_stress_n1, StressVector);
    }
}

} // Namespace Kratos

#ifdef DEBUG_CONSTITUTIVE_LAW
#undef DEBUG_CONSTITUTIVE_LAW
#endif
