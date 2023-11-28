/*
LICENSE: see material_point_application/LICENSE.txt
*/
/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 5 Oct 2022 $
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
#include "constitutive_laws/multiplicative_finite_strain_bridging_constitutive_law.h"
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
MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::MultiplicativeFiniteStrainBridgingConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//**********************************************************************
template<int TStressType>
MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::MultiplicativeFiniteStrainBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : ConstitutiveLaw(), mpConstitutiveLaw(pConstitutiveLaw)
{
}

//**********************************************************************
template<int TStressType>
MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::~MultiplicativeFiniteStrainBridgingConstitutiveLaw()
{
}

//**********************************************************************
template<int TStressType>
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::Has( const Variable<int>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
template<int TStressType>
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::Has( const Variable<double>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
template<int TStressType>
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::Has( const Variable<Vector>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
template<int TStressType>
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::Has( const Variable<Matrix>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
template<int TStressType>
int& MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
template<int TStressType>
double& MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
template<int TStressType>
Vector& MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
template<int TStressType>
Matrix& MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
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
        noalias(rValue) = m_stress_n1;
        return rValue;
    }

    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::SetValue( const Variable<ConstitutiveLaw::Pointer>& rThisVariable, ConstitutiveLaw::Pointer rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == CONSTITUTIVE_LAW)
        mpConstitutiveLaw = rValue;
    else
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetStrainSize(dim);

    mpConstitutiveLaw->InitializeMaterial(props, geom, ShapeFunctionsValues);

    m_F_n.resize(3, 3, false);
    m_F_n1.resize(3, 3, false);
    noalias(m_F_n) = IdentityMatrix(3);
    noalias(m_F_n1) = m_F_n;

    m_left_elastic_cauchy_green_tensor_trial.resize(3, 3, false);
    noalias(m_left_elastic_cauchy_green_tensor_trial) = IdentityMatrix(3);

    m_stress_n1.resize(3, 3, false);
    noalias(m_stress_n1) = ZeroMatrix(3, 3);

    mPrestressFactor = 1.0;
    mPrestress.resize(strain_size, false);
    noalias(mPrestress) = ZeroVector(strain_size);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->ResetMaterial(props, geom, ShapeFunctionsValues);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);

    if (mpConstitutiveLaw->GetStrainMeasure() == ConstitutiveLaw::StrainMeasure_Infinitesimal)
    {
        // reset the trial left Cauchy-Green tensor
        Matrix elastic_strain_tensor(3, 3), left_elastic_cauchy_green_tensor_n(3, 3);

        if (!mpConstitutiveLaw->Has(ELASTIC_STRAIN_TENSOR))
            KRATOS_ERROR << "Constitutive law is not able to return ELASTIC_STRAIN_TENSOR";
        mpConstitutiveLaw->GetValue(ELASTIC_STRAIN_TENSOR, elastic_strain_tensor);

        EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2, left_elastic_cauchy_green_tensor_n, elastic_strain_tensor);

        noalias(m_left_elastic_cauchy_green_tensor_trial) = left_elastic_cauchy_green_tensor_n;
    }
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);

    noalias(m_F_n) = m_F_n1;
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::UpdateDeformationGradient(Matrix& F_new, double& J_new, const Parameters& rValues) const
{
    const Matrix& F = rValues.GetDeformationGradientF();
    if (F.size1() == 3)
    {
        noalias(F_new) = F;
    }
    else if (F.size1() == 2)
    {
        F_new.clear();
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                F_new(i, j) = F(i, j);
        F_new(2, 2) = 1.0;
    }
    J_new = rValues.GetDeterminantF();
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeLogarithmicStrain(Vector& StrainVector,
        Matrix& left_elastic_cauchy_green_tensor_trial, const Matrix& F) const
{
    const unsigned int strain_size = StrainVector.size();

    /*
     * Integration algorithm for isotropic multiplicative finite strain elastoplasticity
     * Reference: Souza de Neto, Computational Plasticity, Box 14.3
     */

    // compute elastic left Cauchy-Green tensor from previous step
    // It is noted that the constitutive law must be able to provide ELASTIC_STRAIN_TENSOR
    Matrix elastic_strain_tensor(3, 3), elastic_strain_tensor_trial(3, 3), left_elastic_cauchy_green_tensor_n(3, 3);

    mpConstitutiveLaw->GetValue(ELASTIC_STRAIN_TENSOR, elastic_strain_tensor);
    // EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2, left_elastic_cauchy_green_tensor_n, elastic_strain_tensor);

    // special treatment when the elastic strain tensor is very small, e.g. at the beginning of the analysis
    const double norm_es = norm_frobenius(elastic_strain_tensor);
    if (norm_es > 1.0e-10)
    {
        EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2,
                left_elastic_cauchy_green_tensor_n, elastic_strain_tensor);
    }
    else
    {
        noalias(left_elastic_cauchy_green_tensor_n) = IdentityMatrix(3);
    }

    // compute trial elastic left Cauchy-Green tensor
    Matrix Fincr(3, 3), invFn(3, 3);
    double detFn;
    MathUtils<double>::InvertMatrix(m_F_n, invFn, detFn);
    noalias(Fincr) = prod(F, invFn);

    noalias(left_elastic_cauchy_green_tensor_trial) = prod(Fincr, Matrix(prod(left_elastic_cauchy_green_tensor_n, trans(Fincr))));

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == DEBUG_ELEMENT_ID && GaussId == 0)
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(elastic_strain_tensor)
        KRATOS_WATCH(left_elastic_cauchy_green_tensor_n)
        KRATOS_WATCH(Fincr)
        KRATOS_WATCH(m_J_n1)
        KRATOS_WATCH(left_elastic_cauchy_green_tensor_trial)
    }
    #endif

    // compute trial elastic logarithmic strain tensor
    std::vector<double> left_elastic_cauchy_green_tensor_trial_pri(3);
    std::vector<Matrix> left_elastic_cauchy_green_tensor_trial_eigprj(3);
    EigenUtility::calculate_principle_stresses(
            left_elastic_cauchy_green_tensor_trial(0, 0),
            left_elastic_cauchy_green_tensor_trial(1, 1),
            left_elastic_cauchy_green_tensor_trial(2, 2),
            left_elastic_cauchy_green_tensor_trial(0, 1),
            left_elastic_cauchy_green_tensor_trial(1, 2),
            left_elastic_cauchy_green_tensor_trial(0, 2),
            left_elastic_cauchy_green_tensor_trial_pri[0],
            left_elastic_cauchy_green_tensor_trial_pri[1],
            left_elastic_cauchy_green_tensor_trial_pri[2],
            left_elastic_cauchy_green_tensor_trial_eigprj[0],
            left_elastic_cauchy_green_tensor_trial_eigprj[1],
            left_elastic_cauchy_green_tensor_trial_eigprj[2]);
    SD_MathUtils<double>::ComputeIsotropicTensorFunction(EigenUtility::logd2, elastic_strain_tensor_trial, left_elastic_cauchy_green_tensor_trial_pri, left_elastic_cauchy_green_tensor_trial_eigprj);
    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == DEBUG_ELEMENT_ID && GaussId == 0)
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(m_F_n)
        KRATOS_WATCH(m_F_n1)
        KRATOS_WATCH(Fincr)
        KRATOS_WATCH(left_elastic_cauchy_green_tensor_n)
        KRATOS_WATCH(left_elastic_cauchy_green_tensor_trial)
        std::cout << "left_elastic_cauchy_green_tensor_trial_pri:";
        for (unsigned int i = 0; i < left_elastic_cauchy_green_tensor_trial_pri.size(); ++i)
        {
            std::cout << " " << left_elastic_cauchy_green_tensor_trial_pri[i];
        }
        std::cout << std::endl;
        KRATOS_WATCH(elastic_strain_tensor_trial)
    }
    #endif

    // create the strain vector as input to the small strain constitutive law
    Vector IncrementalStrainVector(strain_size);
    SD_MathUtils<double>::StrainTensorToVector(elastic_strain_tensor_trial - elastic_strain_tensor, IncrementalStrainVector);
    mpConstitutiveLaw->GetValue(STRAIN, StrainVector);
    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(StrainVector)
        KRATOS_WATCH(IncrementalStrainVector)
    }
    #endif
    noalias(StrainVector) += IncrementalStrainVector;

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == DEBUG_ELEMENT_ID && GaussId == 0)
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(StrainVector)
    }
    #endif
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::StressIntegration(const Parameters& rValues,
        const Matrix& F, Matrix& stress_tensor, Matrix& left_elastic_cauchy_green_tensor_trial) const
{
    const auto& geom = rValues.GetElementGeometry();
    const auto strain_measure = mpConstitutiveLaw->GetStrainMeasure();

    #ifdef DEBUG_CONSTITUTIVE_LAW
    int ElemId, GaussId;
    mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
    mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, GaussId);
    // std::cout << std::setprecision(16);
    // ElemId = 1;
    // GaussId = 0;
    #endif

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == DEBUG_ELEMENT_ID && GaussId == 0)
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(F)
    }
    #endif

    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetStrainSize(dim);

    Vector StrainVector(strain_size);
    if (strain_measure == ConstitutiveLaw::StrainMeasure_Infinitesimal)
    {
        // compute the equivalent (small) strain vector
        this->ComputeLogarithmicStrain(StrainVector, left_elastic_cauchy_green_tensor_trial, F);
    }
    else if (strain_measure == ConstitutiveLaw::StrainMeasure_GreenLagrange)
    {
        // compute the Green-Lagrange strain
        Matrix C;
        if (dim == 2)
        {
            Matrix F2d(2, 2);
            noalias(F2d) = subrange(F, 0, 2, 0, 2);
            C = prod( trans( F2d ), F2d );
        }
        else if (dim == 3)
        {
            C = prod( trans( F ), F );
        }
        this->ComputeGreenLagrangeStrain(StrainVector, C);
    }
    else
        KRATOS_ERROR << "Strain measure " << strain_measure << " is not supported";

    // integrate the constitutive law, obtaining stress
    Vector StressVector(strain_size);
    Matrix Dmat(strain_size, strain_size);
    //
    ConstitutiveLaw::Parameters const_params;
    const_params.SetStrainVector(StrainVector);
    const_params.SetDeformationGradientF(F);
    const_params.SetStressVector(StressVector);
    const_params.SetConstitutiveMatrix(Dmat);
    const_params.SetProcessInfo(rValues.GetProcessInfo());
    const_params.SetMaterialProperties(rValues.GetMaterialProperties());
    const_params.SetElementGeometry(rValues.GetElementGeometry());

    if (strain_measure == ConstitutiveLaw::StrainMeasure_Infinitesimal)
    {
        // for small strain we always integrate the Kirchhoff/Cauchy stress using Cauchy stress integration routine
        mpConstitutiveLaw->CalculateMaterialResponseCauchy(const_params);

        // obtain the stress (1:Cauchy stress; 2: Kirchhoff stress)
        if (!mpConstitutiveLaw->Has(CAUCHY_STRESS_TENSOR))
            KRATOS_ERROR << "Constitutive law is not able to return CAUCHY_STRESS_TENSOR";
        mpConstitutiveLaw->GetValue(CAUCHY_STRESS_TENSOR, stress_tensor);
    }
    else if (strain_measure == ConstitutiveLaw::StrainMeasure_GreenLagrange)
    {
        // integrate the PK2 stress
        mpConstitutiveLaw->CalculateMaterialResponsePK2(const_params);

        // obtain PK2 stress
        Matrix pk2_stress(3, 3);
        mpConstitutiveLaw->GetValue(PK2_STRESS_TENSOR, pk2_stress);

        // transform to Kirchhoff/Cauchy stress
        this->ComputeStress(stress_tensor, pk2_stress);
    }
}

//**********************************************************************
template<>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<1>::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    // integrate the Cauchy stress
    const Matrix& F = rValues.GetDeformationGradientF();
    this->UpdateDeformationGradient(m_F_n1, m_J_n1, rValues);
    this->StressIntegration(rValues, m_F_n1, m_stress_n1, m_left_elastic_cauchy_green_tensor_trial);

    // export the stress
    if (rValues.IsSetStressVector())
    {
        Vector& CauchyStressVector = rValues.GetStressVector();
        SD_MathUtils<double>::StressTensorToVector(m_stress_n1, CauchyStressVector);
    }

    // export the tangent
    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        this->ComputeTangent(AlgorithmicTangent);
    }
}

//**********************************************************************
template<>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<2>::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    // integrate the Kirchhoff stress
    const Matrix& F = rValues.GetDeformationGradientF();
    this->UpdateDeformationGradient(m_F_n1, m_J_n1, rValues);
    this->StressIntegration(rValues, m_F_n1, m_stress_n1, m_left_elastic_cauchy_green_tensor_trial);

    // transform to Cauchy stress
    m_stress_n1 /= m_J_n1;

    // export the stress
    if (rValues.IsSetStressVector())
    {
        Vector& CauchyStressVector = rValues.GetStressVector();
        SD_MathUtils<double>::StressTensorToVector(m_stress_n1, CauchyStressVector);

        #ifdef DEBUG_CONSTITUTIVE_LAW
        if (ElemId == 200 && GaussId == 0)
        {
            KRATOS_WATCH(CauchyStressVector)
        }
        #endif
    }

    // export the tangent
    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        this->ComputeTangent(AlgorithmicTangent);
    }

    // /// compute the numerical tangent (for debugging)
    // Fourth_Order_Tensor dTaudF;
    // SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dTaudF);
    // this->ComputeStressDerivatives(dTaudF, rValues, m_F_n1, 1e-8);
    // KRATOS_WATCH(dTaudF)

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == DEBUG_ELEMENT_ID && GaussId == 0)
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(m_stress_n1)
        KRATOS_WATCH("-----------")
    }
    #endif
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeTangent(Fourth_Order_Tensor& A) const
{
    const auto strain_measure = mpConstitutiveLaw->GetStrainMeasure();

    if (strain_measure == ConstitutiveLaw::StrainMeasure_Infinitesimal)
    {
        this->ComputeInfinitesimalTangent(A);
    }
    else if (strain_measure == ConstitutiveLaw::StrainMeasure_GreenLagrange)
    {
        this->ComputeGreenLagrangeTangent(A);
    }
}

//**********************************************************************
template<>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<1>::ComputeInfinitesimalTangent(Fourth_Order_Tensor& A) const
{
    Fourth_Order_Tensor D, L, B;
    this->ComputeInfinitesimalTangentTerms(D, L, B);

    Fourth_Order_Tensor DL;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(DL);
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0, D, L, DL);

    // compute tangent tensor A
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0, DL, B, A);

    const Matrix eye = IdentityMatrix(3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    A[i][j](k, l) += (m_stress_n1(i, j) * eye(k, l) - m_stress_n1(i, l) * eye(j, k));
}

template<>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<2>::ComputeInfinitesimalTangent(Fourth_Order_Tensor& A) const
{
    Fourth_Order_Tensor D, L, B;
    this->ComputeInfinitesimalTangentTerms(D, L, B);

    double J = MathUtils<double>::Det(m_F_n1);
    Fourth_Order_Tensor DL;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(DL);
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0/J, D, L, DL);

    // compute tangent tensor A
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0, DL, B, A);

    #ifdef DEBUG_CONSTITUTIVE_LAW
    const unsigned int g_size = 9;
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(g_size)
        KRATOS_WATCH(m_left_elastic_cauchy_green_tensor_trial)

        Matrix Du(g_size, g_size);
        SD_MathUtils<double>::TensorToUnsymmetricMatrix(D, Du);
        KRATOS_WATCH(Du)

        Matrix Ls(3, 3);
        SD_MathUtils<double>::TensorToMatrix(L, Ls);
        KRATOS_WATCH(Ls)

        Matrix Lu(g_size, g_size);
        SD_MathUtils<double>::TensorToUnsymmetricMatrix(L, Lu);
        KRATOS_WATCH(Lu/J)

        Matrix Bu(g_size, g_size);
        SD_MathUtils<double>::TensorToUnsymmetricMatrix(B, Bu);
        KRATOS_WATCH(Bu)

        Matrix DLB(g_size, g_size);
        SD_MathUtils<double>::TensorToUnsymmetricMatrix(A, DLB);
        KRATOS_WATCH(DLB)
    }
    #endif

    const Matrix eye = IdentityMatrix(3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    A[i][j](k, l) -= m_stress_n1(i, l) * eye(j, k);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeTangent(Matrix& AlgorithmicTangent) const
{
    Fourth_Order_Tensor A;
    this->ComputeTangent(A);
    SD_MathUtils<double>::TensorToUnsymmetricMatrix(A, AlgorithmicTangent);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeInfinitesimalTangentTerms(Fourth_Order_Tensor& D, Fourth_Order_Tensor& L, Fourth_Order_Tensor& B) const
{
    if (!mpConstitutiveLaw->Has(THREED_ALGORITHMIC_TANGENT))
        KRATOS_ERROR << "Constitutive law is not able to return THREED_ALGORITHMIC_TANGENT";

    // obtain the tangent from the small strain constitutive law. It must be from
    // a fourth order tensor to not missing the out-of-plane component in plane strain
    // analysis
    Matrix Dmat(6, 6);
    mpConstitutiveLaw->GetValue(THREED_ALGORITHMIC_TANGENT, Dmat);

    SD_MathUtils<double>::MatrixToTensor(Dmat, D);

    #ifdef DEBUG_CONSTITUTIVE_LAW
    int ElemId, GaussId;
    mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
    mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, GaussId);
    std::cout << std::setprecision(16);

    if (ElemId == DEBUG_ELEMENT_ID)
    {
        Matrix Dmat2d(3, 3);
        mpConstitutiveLaw->GetValue(ALGORITHMIC_TANGENT, Dmat2d);

        KRATOS_WATCH(Dmat2d)
        KRATOS_WATCH(Dmat)
        KRATOS_WATCH(D)
    }
    // double E, Nu;
    // mpConstitutiveLaw->GetValue(YOUNG_MODULUS, E);
    // mpConstitutiveLaw->GetValue(POISSON_RATIO, Nu);
    // Fourth_Order_Tensor De;
    // SD_MathUtils<double>::CalculateElasticTensor(De, E, Nu);
    // // KRATOS_WATCH(De)
    #endif

    // spectral decomposition for B trial
    std::vector<double> left_elastic_cauchy_green_tensor_trial_pri(3);
    std::vector<Matrix> left_elastic_cauchy_green_tensor_trial_eigprj(3);
    EigenUtility::calculate_principle_stresses(
            m_left_elastic_cauchy_green_tensor_trial(0, 0),
            m_left_elastic_cauchy_green_tensor_trial(1, 1),
            m_left_elastic_cauchy_green_tensor_trial(2, 2),
            m_left_elastic_cauchy_green_tensor_trial(0, 1),
            m_left_elastic_cauchy_green_tensor_trial(1, 2),
            m_left_elastic_cauchy_green_tensor_trial(0, 2),
            left_elastic_cauchy_green_tensor_trial_pri[0],
            left_elastic_cauchy_green_tensor_trial_pri[1],
            left_elastic_cauchy_green_tensor_trial_pri[2],
            left_elastic_cauchy_green_tensor_trial_eigprj[0],
            left_elastic_cauchy_green_tensor_trial_eigprj[1],
            left_elastic_cauchy_green_tensor_trial_eigprj[2]);

    // Perform extra kinematical operations required by materials of this
    // class at large strains for computation of the spatial modulus 'a'
    // Eq. (14.95), Computational Plasticity, de Souza Neto
    SD_MathUtils<double>::ComputeDerivativeIsotropicTensorFunction(EigenUtility::logd2, EigenUtility::dlogd2,
        L,
        m_left_elastic_cauchy_green_tensor_trial,
        left_elastic_cauchy_green_tensor_trial_pri,
        left_elastic_cauchy_green_tensor_trial_eigprj,
        1e-10); // tolerance to compare the eigenvalues

    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(B);
    const Matrix eye = IdentityMatrix(3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    B[i][j](k, l) = eye(i, k) * m_left_elastic_cauchy_green_tensor_trial(j, l)
                                  + eye(j, k) * m_left_elastic_cauchy_green_tensor_trial(i, l);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeGreenLagrangeTangent(Fourth_Order_Tensor& A) const
{
    if (!mpConstitutiveLaw->Has(THREED_ALGORITHMIC_TANGENT))
        KRATOS_ERROR << "Constitutive law is not able to return THREED_ALGORITHMIC_TANGENT";

    // obtain the tangent from the small strain constitutive law. It must be from
    // a fourth order tensor to not missing the out-of-plane component in plane strain
    // analysis
    Matrix Dmat(6, 6);
    mpConstitutiveLaw->GetValue(THREED_ALGORITHMIC_TANGENT, Dmat);

    Fourth_Order_Tensor D;
    SD_MathUtils<double>::MatrixToTensor(Dmat, D);

    // compute the tensor a
    const Matrix& F = m_F_n1;
    double J = m_J_n1;
    // KRATOS_WATCH(J)

    const Matrix& sigma = m_stress_n1;

    double aux;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    const Matrix eye = IdentityMatrix(3);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                for (int l = 0; l < 3; ++l)
                {
                    aux = eye(i, k) * sigma(l, j);

                    for (int m = 0; m < 3; ++m)
                        for (int u = 0; u < 3; ++u)
                            for (int p = 0; p < 3; ++p)
                                for (int q = 0; q < 3; ++q)
                                    aux += 0.5/J * F(i, m) * F(j, u) * D[m][u](p, q)
                                            * (F(l, p)*F(k, q) + F(l, q)*F(k, p));

                    A[i][j](k, l) = aux;
                }
            }
        }
    }
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeGreenLagrangeStrain( Vector& StrainVector, const Matrix& C ) const
{
    KRATOS_TRY

    const unsigned int dimension = C.size1();

    if ( dimension == 2 )
    {
        if ( StrainVector.size() != 3 ) StrainVector.resize( 3, false );

        StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        StrainVector[2] = C( 0, 1 );
    }

    if ( dimension == 3 )
    {
        if ( StrainVector.size() != 6 ) StrainVector.resize( 6, false );

        StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        StrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

        StrainVector[3] = C( 0, 1 ); // xy

        StrainVector[4] = C( 1, 2 ); // yz

        StrainVector[5] = C( 0, 2 ); // xz
    }

    KRATOS_CATCH( "" )
}

//**********************************************************************
template<>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<1>::ComputeStress(Matrix& stress_tensor, const Matrix& PK2_stress) const
{
    noalias(stress_tensor) = 1.0/m_J_n1 * prod(m_F_n1, Matrix(prod(PK2_stress, trans(m_F_n1))));
}

template<>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<2>::ComputeStress(Matrix& stress_tensor, const Matrix& PK2_stress) const
{
    noalias(stress_tensor) = prod(m_F_n1, Matrix(prod(PK2_stress, trans(m_F_n1))));
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeBetrialDerivatives(Fourth_Order_Tensor& B, const Matrix& F, double epsilon) const
{
    Matrix elastic_strain_tensor(3, 3), left_elastic_cauchy_green_tensor_n(3, 3);

    mpConstitutiveLaw->GetValue(ELASTIC_STRAIN_TENSOR, elastic_strain_tensor);

    // special treatment when the elastic strain tensor is very small, e.g. at the beginning of the analysis
    const double norm_es = norm_frobenius(elastic_strain_tensor);
    if (norm_es > 1.0e-10)
    {
        EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2,
                left_elastic_cauchy_green_tensor_n, elastic_strain_tensor);
    }
    else
    {
        noalias(left_elastic_cauchy_green_tensor_n) = IdentityMatrix(3);
    }

    Matrix Fincr(3, 3), invFn(3, 3);
    double detFn;
    MathUtils<double>::InvertMatrix(m_F_n, invFn, detFn);

    // compute the Be^trial
    Matrix left_elastic_cauchy_green_tensor_trial(3, 3);

    noalias(Fincr) = prod(F, invFn);
    noalias(left_elastic_cauchy_green_tensor_trial) = prod(Fincr, Matrix(prod(left_elastic_cauchy_green_tensor_n, trans(Fincr))));
    KRATOS_WATCH(left_elastic_cauchy_green_tensor_trial)

    // compute the new Be^trial and the numerical derivatives
    Matrix new_left_elastic_cauchy_green_tensor_trial(3, 3), newF(3, 3), aux(3, 3);

    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newF) = F;
            newF(k, l) += epsilon;

            noalias(Fincr) = prod(newF, invFn);
            noalias(new_left_elastic_cauchy_green_tensor_trial) = prod(Fincr, Matrix(prod(left_elastic_cauchy_green_tensor_n, trans(Fincr))));

            noalias(aux) = (new_left_elastic_cauchy_green_tensor_trial - left_elastic_cauchy_green_tensor_trial) / epsilon;

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    B[i][j](k, l) = aux(i, j);
                }
            }
        }
    }
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeStressDerivatives(Fourth_Order_Tensor& D,
        const Parameters& rValues, const Matrix& F, double epsilon) const
{
    // integrate the stress (1: Cauchy; 2: Kirchhoff)
    Matrix stress_tensor(3, 3), dummy(3, 3);
    this->StressIntegration(rValues, F, stress_tensor, dummy);

    Matrix newF(3, 3), new_stress_tensor(3, 3), aux(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newF) = F;
            newF(k, l) += epsilon;

            // integrate the new stress
            this->StressIntegration(rValues, newF, new_stress_tensor, dummy);

            noalias(aux) = (new_stress_tensor - stress_tensor) / epsilon;

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    D[i][j](k, l) = aux(i, j);
                }
            }
        }
    }
}

//**********************************************************************
template<int TStressType>
int MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    return mpConstitutiveLaw->Check( props, geom, CurrentProcessInfo );
}

//**********************************************************************

template class MultiplicativeFiniteStrainBridgingConstitutiveLaw<1>;
template class MultiplicativeFiniteStrainBridgingConstitutiveLaw<2>;

} // Namespace Kratos

#undef DEBUG_CONSTITUTIVE_LAW
