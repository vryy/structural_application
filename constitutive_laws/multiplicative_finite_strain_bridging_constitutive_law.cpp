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
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    BaseType::InitializeMaterial(props, geom, ShapeFunctionsValues);

    m_Be_trial.resize(3, 3, false);
    noalias(m_Be_trial) = IdentityMatrix(3);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    BaseType::InitializeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);

    // reset the trial left Cauchy-Green tensor
    Matrix elastic_strain_tensor(3, 3), Be_n(3, 3);

    if (!mpConstitutiveLaw->Has(ELASTIC_STRAIN_TENSOR))
        KRATOS_ERROR << "Constitutive law is not able to return ELASTIC_STRAIN_TENSOR";
    mpConstitutiveLaw->GetValue(ELASTIC_STRAIN_TENSOR, elastic_strain_tensor);

    EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2, Be_n, elastic_strain_tensor);

    // since the tangent depends on m_Be_trial, resetting it
    // to the last converged value will reset the tangent to the last one
    noalias(m_Be_trial) = Be_n;
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeStrain(Vector& StrainVector,
        Matrix& Be_trial, const Matrix& F) const
{
    const unsigned int strain_size = StrainVector.size();

    /*
     * Integration algorithm for isotropic multiplicative finite strain elastoplasticity
     * Reference: Souza de Neto, Computational Plasticity, Box 14.3
     */

    // compute elastic left Cauchy-Green tensor from previous step
    // It is noted that the constitutive law must be able to provide ELASTIC_STRAIN_TENSOR
    Matrix elastic_strain_tensor(3, 3), elastic_strain_tensor_trial(3, 3), Be_n(3, 3);

    mpConstitutiveLaw->GetValue(ELASTIC_STRAIN_TENSOR, elastic_strain_tensor);
    // EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2, Be_n, elastic_strain_tensor);

    // special treatment when the elastic strain tensor is very small, e.g. at the beginning of the analysis
    const double norm_es = norm_frobenius(elastic_strain_tensor);
    if (norm_es > 1.0e-10)
    {
        EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2,
                Be_n, elastic_strain_tensor);
    }
    else
    {
        noalias(Be_n) = IdentityMatrix(3);
    }

    // compute trial elastic left Cauchy-Green tensor
    Matrix Fincr(3, 3), invFn(3, 3);
    double detFn;
    MathUtils<double>::InvertMatrix(m_F_n, invFn, detFn);
    noalias(Fincr) = prod(F, invFn);

    noalias(Be_trial) = prod(Fincr, Matrix(prod(Be_n, trans(Fincr))));

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == DEBUG_ELEMENT_ID && GaussId == 0)
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(elastic_strain_tensor)
        KRATOS_WATCH(Be_n)
        KRATOS_WATCH(Fincr)
        KRATOS_WATCH(m_J_n1)
        KRATOS_WATCH(Be_trial)
    }
    #endif

    // compute trial elastic logarithmic strain tensor
    std::vector<double> Be_trial_pri(3);
    std::vector<Matrix> Be_trial_eigprj(3);
    EigenUtility::calculate_principle_stresses(
            Be_trial(0, 0),
            Be_trial(1, 1),
            Be_trial(2, 2),
            Be_trial(0, 1),
            Be_trial(1, 2),
            Be_trial(0, 2),
            Be_trial_pri[0],
            Be_trial_pri[1],
            Be_trial_pri[2],
            Be_trial_eigprj[0],
            Be_trial_eigprj[1],
            Be_trial_eigprj[2]);
    SD_MathUtils<double>::ComputeIsotropicTensorFunction(EigenUtility::logd2, elastic_strain_tensor_trial, Be_trial_pri, Be_trial_eigprj);
    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == DEBUG_ELEMENT_ID && GaussId == 0)
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(m_F_n)
        KRATOS_WATCH(m_F_n1)
        KRATOS_WATCH(Fincr)
        KRATOS_WATCH(Be_n)
        KRATOS_WATCH(Be_trial)
        std::cout << "Be_trial_pri:";
        for (unsigned int i = 0; i < Be_trial_pri.size(); ++i)
        {
            std::cout << " " << Be_trial_pri[i];
        }
        std::cout << std::endl;
        KRATOS_WATCH(elastic_strain_tensor_trial)
    }
    #endif

    // create the strain vector as input to the small strain constitutive law
    Vector IncrementalStrainVector(strain_size);
    SD_MathUtils<double>::StrainTensorToVector(elastic_strain_tensor_trial - elastic_strain_tensor, IncrementalStrainVector);
    if (mpConstitutiveLaw->Has(STRAIN_OLD))
        mpConstitutiveLaw->GetValue(STRAIN_OLD, StrainVector);
    else
        KRATOS_ERROR << mpConstitutiveLaw->Info() << " does not provide the STRAIN_OLD variable";
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
        const Matrix& F, Matrix& stress_tensor, Matrix& Be_trial) const
{
    const auto strain_measure = mpConstitutiveLaw->GetStrainMeasure();
    if (strain_measure != ConstitutiveLaw::StrainMeasure_Infinitesimal)
        KRATOS_ERROR << "The strain measure is not StrainMeasure_Infinitesimal";

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

    const auto& geom = rValues.GetElementGeometry();
    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetStrainSize(dim);

    // compute the equivalent (small) strain vector
    Vector StrainVector(strain_size);
    this->ComputeStrain(StrainVector, Be_trial, F);

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

    // for small strain we always integrate the Kirchhoff/Cauchy stress using Cauchy stress integration routine
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(const_params);

    // obtain the stress (1:Cauchy stress; 2: Kirchhoff stress)
    if (!mpConstitutiveLaw->Has(CAUCHY_STRESS_TENSOR))
        KRATOS_ERROR << "Constitutive law is not able to return CAUCHY_STRESS_TENSOR";
    mpConstitutiveLaw->GetValue(CAUCHY_STRESS_TENSOR, stress_tensor);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    // integrate the Kirchhoff stress
    this->UpdateDeformationGradient(m_F_n1, m_J_n1, rValues);
    this->StressIntegration(rValues, m_F_n1, m_stress_n1, m_Be_trial);

    if constexpr (TStressType == 2)
    {
        // transform to Cauchy stress
        m_stress_n1 /= m_J_n1;
    }

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
        BaseType::ComputeTangent(AlgorithmicTangent);
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
template<>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<1>::ComputeTangent(Fourth_Order_Tensor& A) const
{
    Fourth_Order_Tensor D, L, B;
    this->ComputeTangentTerms(D, L, B);

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
                    A[i][j][k][l] += (m_stress_n1(i, j) * eye(k, l) - m_stress_n1(i, l) * eye(j, k));
}

template<>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<2>::ComputeTangent(Fourth_Order_Tensor& A) const
{
    Fourth_Order_Tensor D, L, B;
    this->ComputeTangentTerms(D, L, B);

    Fourth_Order_Tensor DL;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(DL);
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0, D, L, DL);

    // compute tangent tensor A
    // const double J = MathUtils<double>::Det(m_F_n1);
    const double J = m_J_n1; // it is also OK to use this, even for Fbar
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0/J, DL, B, A);

    #ifdef DEBUG_CONSTITUTIVE_LAW
    const unsigned int g_size = 9;
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(g_size)
        KRATOS_WATCH(m_Be_trial)

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
                    A[i][j][k][l] -= m_stress_n1(i, l) * eye(j, k);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeTangentTerms(Fourth_Order_Tensor& D, Fourth_Order_Tensor& L, Fourth_Order_Tensor& B) const
{
    if (!mpConstitutiveLaw->Has(THREED_ALGORITHMIC_TANGENT))
        KRATOS_ERROR << "Constitutive law is not able to return THREED_ALGORITHMIC_TANGENT";

    // obtain the tangent from the small strain constitutive law. It must be from
    // a fourth order tensor to not missing the out-of-plane component in plane strain
    // analysis
    Matrix Dmat(6, 6);
    mpConstitutiveLaw->GetValue(THREED_ALGORITHMIC_TANGENT, Dmat);

    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(D);
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
    std::vector<double> Be_trial_pri(3);
    std::vector<Matrix> Be_trial_eigprj(3);
    EigenUtility::calculate_principle_stresses(
            m_Be_trial(0, 0),
            m_Be_trial(1, 1),
            m_Be_trial(2, 2),
            m_Be_trial(0, 1),
            m_Be_trial(1, 2),
            m_Be_trial(0, 2),
            Be_trial_pri[0],
            Be_trial_pri[1],
            Be_trial_pri[2],
            Be_trial_eigprj[0],
            Be_trial_eigprj[1],
            Be_trial_eigprj[2]);

    // Perform extra kinematical operations required by materials of this
    // class at large strains for computation of the spatial modulus 'a'
    // Eq. (14.95), Computational Plasticity, de Souza Neto
    SD_MathUtils<double>::ComputeDerivativeIsotropicTensorFunction(EigenUtility::logd2, EigenUtility::dlogd2,
        L,
        m_Be_trial,
        Be_trial_pri,
        Be_trial_eigprj,
        1e-10); // tolerance to compare the eigenvalues

    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(B);
    const Matrix eye = IdentityMatrix(3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    B[i][j][k][l] = eye(i, k) * m_Be_trial(j, l)
                                  + eye(j, k) * m_Be_trial(i, l);
}

//**********************************************************************
template<int TStressType>
void MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>::ComputeBetrialDerivatives(Fourth_Order_Tensor& B, const Matrix& F, double epsilon) const
{
    Matrix elastic_strain_tensor(3, 3), Be_n(3, 3);

    mpConstitutiveLaw->GetValue(ELASTIC_STRAIN_TENSOR, elastic_strain_tensor);

    // special treatment when the elastic strain tensor is very small, e.g. at the beginning of the analysis
    const double norm_es = norm_frobenius(elastic_strain_tensor);
    if (norm_es > 1.0e-10)
    {
        EigenUtility::ComputeIsotropicTensorFunction(EigenUtility::exp2,
                Be_n, elastic_strain_tensor);
    }
    else
    {
        noalias(Be_n) = IdentityMatrix(3);
    }

    Matrix Fincr(3, 3), invFn(3, 3);
    double detFn;
    MathUtils<double>::InvertMatrix(m_F_n, invFn, detFn);

    // compute the Be^trial
    Matrix Be_trial(3, 3);

    noalias(Fincr) = prod(F, invFn);
    noalias(Be_trial) = prod(Fincr, Matrix(prod(Be_n, trans(Fincr))));
    KRATOS_WATCH(Be_trial)

    // compute the new Be^trial and the numerical derivatives
    Matrix new_Be_trial(3, 3), newF(3, 3), aux(3, 3);

    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newF) = F;
            newF(k, l) += epsilon;

            noalias(Fincr) = prod(newF, invFn);
            noalias(new_Be_trial) = prod(Fincr, Matrix(prod(Be_n, trans(Fincr))));

            noalias(aux) = (new_Be_trial - Be_trial) / epsilon;

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    B[i][j][k][l] = aux(i, j);
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
                    D[i][j][k][l] = aux(i, j);
                }
            }
        }
    }
}

//**********************************************************************
template<>
std::string MultiplicativeFiniteStrainBridgingConstitutiveLaw<1>::Info() const
{
    return "MultiplicativeFiniteStrainBridgingConstitutiveLaw<Cauchy>";
}

template<>
std::string MultiplicativeFiniteStrainBridgingConstitutiveLaw<2>::Info() const
{
    return "MultiplicativeFiniteStrainBridgingConstitutiveLaw<Kirchhoff>";
}

//**********************************************************************

template class MultiplicativeFiniteStrainBridgingConstitutiveLaw<1>;
template class MultiplicativeFiniteStrainBridgingConstitutiveLaw<2>;

} // Namespace Kratos

#undef DEBUG_CONSTITUTIVE_LAW
