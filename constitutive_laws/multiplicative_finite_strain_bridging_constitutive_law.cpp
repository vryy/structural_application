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

//#define DEBUG_CONSTITUTIVE_LAW

#ifdef DEBUG_CONSTITUTIVE_LAW
#include <iomanip>
#endif

namespace Kratos
{

//**********************************************************************
MultiplicativeFiniteStrainBridgingConstitutiveLaw::MultiplicativeFiniteStrainBridgingConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//**********************************************************************
MultiplicativeFiniteStrainBridgingConstitutiveLaw::MultiplicativeFiniteStrainBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : ConstitutiveLaw(), mpConstitutiveLaw(pConstitutiveLaw)
{
}

//**********************************************************************
MultiplicativeFiniteStrainBridgingConstitutiveLaw::~MultiplicativeFiniteStrainBridgingConstitutiveLaw()
{
}

//**********************************************************************
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw::Has( const Variable<int>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw::Has( const Variable<double>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw::Has( const Variable<Vector>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
int& MultiplicativeFiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
double& MultiplicativeFiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Vector& MultiplicativeFiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Matrix& MultiplicativeFiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
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
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<ConstitutiveLaw::Pointer>& rThisVariable, ConstitutiveLaw::Pointer rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == CONSTITUTIVE_LAW)
        mpConstitutiveLaw = rValue;
    else
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = dim*(dim+1) / 2;

    mpConstitutiveLaw->InitializeMaterial(props, geom, ShapeFunctionsValues);

    mLastF.resize(3, 3, false);
    mCurrentF.resize(3, 3, false);
    noalias(mLastF) = IdentityMatrix(3);
    noalias(mCurrentF) = mLastF;

    m_left_cauchy_green_tensor_trial.resize(3, 3, false);
    noalias(m_left_cauchy_green_tensor_trial) = IdentityMatrix(3);

    m_stress_n1.resize(3, 3, false);
    noalias(m_stress_n1) = ZeroMatrix(3, 3);

    mPrestressFactor = 1.0;
    mPrestress.resize(strain_size, false);
    noalias(mPrestress) = ZeroVector(strain_size);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->ResetMaterial(props, geom, ShapeFunctionsValues);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);

    // reset the trial left Cauchy-Green tensor
    Matrix elastic_strain_tensor(3, 3), left_cauchy_green_tensor_n(3, 3);

    mpConstitutiveLaw->GetValue(ELASTIC_STRAIN_TENSOR, elastic_strain_tensor);

    EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2, left_cauchy_green_tensor_n, elastic_strain_tensor);

    noalias(m_left_cauchy_green_tensor_trial) = left_cauchy_green_tensor_n;
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);

    noalias(mLastF) = mCurrentF;
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    const auto& geom = rValues.GetElementGeometry();

    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = dim*(dim+1) / 2;

    #ifdef DEBUG_CONSTITUTIVE_LAW
    int ElemId, GaussId;
    mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
    mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, GaussId);
    std::cout << std::setprecision(16);
    #endif

    const Matrix& F = rValues.GetDeformationGradientF();
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
    mCurrentDetF = rValues.GetDeterminantF();

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

    // special treatment when the elastic strain tensor is very small, e.g. at the beginning of the analysis
    const double norm_es = norm_frobenius(elastic_strain_tensor);
    if (norm_es > 1.0e-10)
    {
        EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2,
                left_cauchy_green_tensor_n, elastic_strain_tensor);
    }
    else
    {
        noalias(left_cauchy_green_tensor_n) = IdentityMatrix(3);
    }

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

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(const_params);

    // transform the stress to Cauchy stress
    if (rValues.IsSetStressVector())
    {
        mpConstitutiveLaw->GetValue(CAUCHY_STRESS_TENSOR, m_stress_n1);
        m_stress_n1 /= mCurrentDetF;

        Vector& CauchyStressVector = rValues.GetStressVector();
        SD_MathUtils<double>::StressTensorToVector(m_stress_n1, CauchyStressVector);

        #ifdef DEBUG_CONSTITUTIVE_LAW
        if (ElemId == 200 && GaussId == 0)
        {
            KRATOS_WATCH(CauchyStressVector)
        }
        #endif
    }

    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        this->ComputeTangent(AlgorithmicTangent);
    }

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == 1 && GaussId == 0)
    if (ElemId == 1)
    {
        KRATOS_WATCH(m_stress_n1)
        KRATOS_WATCH("-----------")
    }
    #endif
}

void MultiplicativeFiniteStrainBridgingConstitutiveLaw::ComputeTangent(Matrix& AlgorithmicTangent) const
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

    #ifdef DEBUG_CONSTITUTIVE_LAW
    int ElemId, GaussId;
    mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
    mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, GaussId);
    std::cout << std::setprecision(16);

    if (ElemId == 1)
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

    // Perform extra kinematical operations required by materials of this
    // class at large strains for computation of the spatial modulus 'a'
    // Eq. (14.95), Computational Plasticity, de Souza Neto
    Fourth_Order_Tensor L;
    SD_MathUtils<double>::ComputeDerivativeIsotropicTensorFunction(EigenUtility::logd2, EigenUtility::dlogd2,
        L,
        m_left_cauchy_green_tensor_trial,
        left_cauchy_green_tensor_trial_pri,
        left_cauchy_green_tensor_trial_eigprj,
        1e-10); // tolerance to compare the eigenvalues

    Fourth_Order_Tensor B;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(B);
    const Matrix eye = IdentityMatrix(3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    B[i][j](k, l) = eye(i, k) * m_left_cauchy_green_tensor_trial(j, l)
                                  + eye(j, k) * m_left_cauchy_green_tensor_trial(i, l);

    double J = MathUtils<double>::Det(mCurrentF);
    Fourth_Order_Tensor DL;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(DL);
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0/J, D, L, DL);

    Fourth_Order_Tensor A;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0, DL, B, A);

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        KRATOS_WATCH(m_left_cauchy_green_tensor_trial)

        Matrix Du(4, 4);
        SD_MathUtils<double>::TensorToUnsymmetricMatrix(D, Du);
        KRATOS_WATCH(Du)

        Matrix Ls(3, 3);
        SD_MathUtils<double>::TensorToMatrix(L, Ls);
        KRATOS_WATCH(Ls)

        Matrix Lu(4, 4);
        SD_MathUtils<double>::TensorToUnsymmetricMatrix(L, Lu);
        KRATOS_WATCH(Lu/J)

        Matrix Bu(4, 4);
        SD_MathUtils<double>::TensorToUnsymmetricMatrix(B, Bu);
        KRATOS_WATCH(Bu)

        Matrix DLB(4, 4);
        SD_MathUtils<double>::TensorToUnsymmetricMatrix(A, DLB);
        KRATOS_WATCH(DLB)
    }
    #endif

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    A[i][j](k, l) -= m_stress_n1(i, l) * eye(j, k);

    SD_MathUtils<double>::TensorToUnsymmetricMatrix(A, AlgorithmicTangent);
}

//**********************************************************************
int MultiplicativeFiniteStrainBridgingConstitutiveLaw::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    return mpConstitutiveLaw->Check( props, geom, CurrentProcessInfo );
}

} // Namespace Kratos

#ifdef DEBUG_CONSTITUTIVE_LAW
#undef DEBUG_CONSTITUTIVE_LAW
#endif
