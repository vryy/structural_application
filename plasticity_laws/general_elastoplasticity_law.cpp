/*
LICENSE: see soil_mechanics_application/LICENSE.txt
*/
/* *********************************************************
*
*   Last Modified by:    $Author: Giang Bui-Hoang $
*   Date:                $Date: 16 Feb 2020 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "plasticity_laws/general_elastoplasticity_law.h"

#define CHECK_DERIVATIVES
#define OSTR std::cout

namespace Kratos
{

void GeneralElastoplasticityLaw::ComputeContinuumPlasticTangent(Fourth_Order_Tensor& Cep,
        const Fourth_Order_Tensor& Ce, const Matrix& stress, const Vector& q, const Vector& alpha) const
{
    const unsigned int stress_size = stress.size1();
    const unsigned int q_size = q.size();

    Matrix A(stress_size, stress_size), B(stress_size, stress_size);
    Vector c(q_size), b(q_size), d(q_size);
    this->dFdSigma(A, stress, q);
    this->dGdSigma(B, stress, q);
    this->dFdQ(c, stress, q);
    this->dGdQ(b, stress, q);
    this->dQdPhi(d, stress, q, alpha);
    // #ifdef DEBUG_GENERAL_PLASTICITY_LAW
    // KRATOS_WATCH(A)
    // KRATOS_WATCH(B)
    // KRATOS_WATCH(c)
    // KRATOS_WATCH(b)
    // KRATOS_WATCH(d)
    // #endif

    SD_MathUtils<double>::CopyFourthOrderTensor(Ce, Cep);
    // #ifdef DEBUG_GENERAL_PLASTICITY_LAW
    // Matrix De(6, 6);
    // SD_MathUtils<double>::TensorToMatrix(Cep, De);
    // KRATOS_WATCH(De)
    // #endif
    Matrix CeB(3, 3, 0.0), CeA(3, 3, 0.0);
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0, Ce, B, CeB);
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0, Ce, A, CeA);
    double aux = SD_MathUtils<double>::mat_inner_prod(A, CeB);
    aux -= inner_prod(c, d);
    // #ifdef DEBUG_GENERAL_PLASTICITY_LAW
    // KRATOS_WATCH(aux)
    // #endif
    SD_MathUtils<double>::OuterProductFourthOrderTensor(-1.0/aux, CeB, CeA, Cep);
}

std::vector<double> GeneralElastoplasticityLaw::PlasticIntegration_Substepping(Matrix& stress, Vector& q, Vector& alpha,
        const Matrix& incremental_strain, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters, const unsigned int max_level,
        const int debug_level) const
{
    Matrix stress_trial = stress, stress_last = stress;
    Vector q_last = q, alpha_last = alpha;

    int error_code;
    double factor = 0.0, dfactor = 1.0, factor_last = 0.0, dlambda;
    std::vector<double> cfact;
    int level = 0;

    while (factor_last < 1.0)
    {
        // update the factor
        factor = factor_last + dfactor;

        // compute stress
        noalias(stress) = stress_last;
        SD_MathUtils<double>::ContractFourthOrderTensor(dfactor, Ce, incremental_strain, stress);

        noalias(q) = q_last;
        noalias(alpha) = alpha_last;

        // plastic integration
        noalias(stress_trial) = stress;
        error_code = PlasticIntegration(stress, q, alpha, dlambda, stress_trial, Ce,
                FTOL, max_iters, debug_level);

        // adjust step
        if (error_code == 0)
        {
            #ifdef DEBUG_GENERAL_PLASTICITY_LAW
            if (debug_level > 2)
            {
                std::cout << "PlasticIntegration converges at load step " << factor << std::endl;
                KRATOS_WATCH(stress)
                KRATOS_WATCH(q)
                KRATOS_WATCH(alpha)
            }
            #endif

            // save
            cfact.push_back(factor);
            factor_last = factor;
            noalias(stress_last) = stress;
            noalias(q_last) = q;
            noalias(alpha_last) = alpha;

            // restore
            dfactor *= 2;
            if (dfactor > 1.0 - factor)
                dfactor = 1.0 - factor;
            level = 0;
        }
        else
        {
            dfactor *= 0.5;
            ++level;
            if (level == max_level)
                KRATOS_ERROR << "The internal calculation hits maximum number of allowed levels";
        }
    }

    #ifdef DEBUG_GENERAL_PLASTICITY_LAW
    if (debug_level > 0)
    {
        std::cout << "PlasticIntegration with sub-stepping converges after " << cfact.size() << " steps" << std::endl;
        std::cout << "cfact:";
        for (std::size_t i = 0; i < cfact.size(); ++i)
            std::cout << " " << cfact[i];
        std::cout << std::endl;
    }
    #endif

    return std::move(cfact);
}

int GeneralElastoplasticityLaw::PlasticIntegration_Substepping(Matrix& stress, Vector& q, Vector& alpha,
        const std::vector<double>& loads, const Matrix& incremental_strain, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const int debug_level) const
{
    Matrix stress_trial(3, 3);

    int error_code;
    double factor, dfactor, factor_last = 0.0, dlambda;

    for (std::size_t i = 0; i < loads.size(); ++i)
    {
        // update the factor
        factor = loads[i];
        dfactor = factor - factor_last;
        factor_last = factor;

        // compute stress (trial)
        SD_MathUtils<double>::ContractFourthOrderTensor(dfactor, Ce, incremental_strain, stress);

        // determine if this is elastic or plastic step
        double f = this->F(stress, q);
        if (f < FTOL)
        {
            // elastic step, nothing to do
            #ifdef DEBUG_GENERAL_PLASTICITY_LAW
            if (debug_level > 1)
            {
                std::cout << "PlasticIntegration completes at sub-step " << factor
                          << ", elastic state is detected"
                          << std::endl;
            }
            #endif
        }
        else
        {
            // plastic step, plastic integration
            noalias(stress_trial) = stress;
            error_code = PlasticIntegration(stress, q, alpha, dlambda, stress_trial, Ce,
                    FTOL, max_iters, debug_level);

            #ifdef DEBUG_GENERAL_PLASTICITY_LAW
            if (debug_level > 1)
            {
                std::cout << "PlasticIntegration completes at sub-step " << factor << std::endl;
                KRATOS_WATCH(error_code)
                KRATOS_WATCH(stress)
                KRATOS_WATCH(q)
                KRATOS_WATCH(alpha)
            }
            #endif

            // check convergence
            if (error_code != 0)
            {
                KRATOS_WATCH_STD_CON(loads)
                // KRATOS_ERROR << "The plastic integration does not converge at load step " << factor;
                std::cout << "The plastic integration does not converge at load step " << factor << std::endl;
                return error_code;
            }
        }
    }

    return 0;
}

int GeneralElastoplasticityLaw::PlasticIntegration(Matrix& stress, Vector& q, Vector& alpha, double& dlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const int debug_level) const
{
    const unsigned int nvars = this->NumberOfInternalVariables();

    Matrix n(3, 3), m(3, 3), b(3, 3), r(3, 3), dm_dlambda(3, 3), ninvA(3, 3), Ceb(3, 3);
    Vector dfdq(nvars), dqdphi(nvars), dalphadphi(nvars);
    double ddlambda, l, f;
    const Matrix eye = IdentityMatrix(3);

    #if defined(DEBUG_GENERAL_PLASTICITY_LAW) && defined(CHECK_DERIVATIVES)
    Matrix num_n(3, 3);
    Matrix num_m(3, 3);
    Fourth_Order_Tensor num_dm_dsigma;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(num_dm_dsigma);
    #endif

    Third_Order_Tensor dmdq;
    SD_MathUtils<double>::InitializeThirdOrderTensor(dmdq, 3, 3, nvars, false);

    Fourth_Order_Tensor A, invA, dm_dsigma;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(invA);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dm_dsigma);

    #ifdef DEBUG_GENERAL_PLASTICITY_LAW
    Fourth_Order_Tensor AinvA;
    Matrix ninvAA(3, 3);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(AinvA);
    #endif

    bool converged = false;
    int it = 0;
    dlambda = ddlambda = 0.0;
    double norm_r;
    const Vector alphan = alpha;
    double norm_stress_trial = norm_frobenius(stress_trial);
    if (norm_stress_trial < 1.0) norm_stress_trial = 1.0;
    do
    {
        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(it)
        }
        #endif

        // compute f
        f = this->F(stress, q);

        // compute gradient m
        this->dGdSigma(m, stress, q);

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
        }
        #endif

        // compute r
        noalias(r) = stress - stress_trial;
        SD_MathUtils<double>::ContractFourthOrderTensor(dlambda, Ce, m, r);
        norm_r = norm_frobenius(r);

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(f)
            KRATOS_WATCH(norm_r)
            KRATOS_WATCH(FTOL)
        }
        #endif

        // check convergence
        if (std::abs(f) + norm_r/norm_stress_trial < FTOL)
        {
            converged = true;
            break;
        }

        // compute A
        this->d2GdSigma2(dm_dsigma, stress, q);
        SD_MathUtils<double>::ZeroFourthOrderTensor(A); // A = 0
        SD_MathUtils<double>::SpecialProduct1FourthOrderTensor(1.0, eye, eye, A); // A += delta_ik delta_jl
        SD_MathUtils<double>::ProductFourthOrderTensor(dlambda, Ce, dm_dsigma, A); // A += lambda*Ce:dm_dsigma

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            // KRATOS_WATCH(A)
        }
        #endif

        // compute inverse of A
        SD_MathUtils<double>::InvertFourthOrderTensor(A, invA);

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            // SD_MathUtils<double>::ZeroFourthOrderTensor(AinvA);
            // SD_MathUtils<double>::ProductFourthOrderTensor(1.0, A, invA, AinvA);
            // KRATOS_WATCH(AinvA)
        }
        #endif

        // compute gradient n
        this->dFdSigma(n, stress, q);

        #if defined(DEBUG_GENERAL_PLASTICITY_LAW) && defined(CHECK_DERIVATIVES)
        if (debug_level > 2)
        {
            this->Num_dFdSigma(num_n, stress, q, 1e-6);
            KRATOS_WATCH(n)
            KRATOS_WATCH(num_n)
            std::cout << "diff_n: " << norm_frobenius(n-num_n) << std::endl;
            this->Num_dGdSigma(num_m, stress, q, 1e-6);
            // KRATOS_WATCH(m)
            // KRATOS_WATCH(num_m)
            std::cout << "diff_m: " << norm_frobenius(m-num_m) << std::endl;
            this->Num_d2GdSigma2(num_dm_dsigma, stress, q, 1e-6);
            // KRATOS_WATCH(dm_dsigma)
            // KRATOS_WATCH(num_dm_dsigma)
            SD_MathUtils<double>::AddFourthOrderTensor(-1.0, dm_dsigma, num_dm_dsigma);
            const double norm_diff_dm_dsigma = SD_MathUtils<double>::normTensor(num_dm_dsigma);
            std::cout << "diff_dm_dsigma: " << norm_diff_dm_dsigma << std::endl;
        }
        #endif

        // compute b
        this->dQdPhi(dqdphi, stress, q, alpha);
        this->d2GdSigmadQ(dmdq, stress, q);
        dm_dlambda.clear();
        SD_MathUtils<double>::ContractThirdOrderTensor(1.0, dmdq, dqdphi, dm_dlambda);
        noalias(b) = m + dlambda*dm_dlambda;

        // compute l
        this->dFdQ(dfdq, stress, q);
        l = inner_prod(dfdq, dqdphi);

        #if defined(DEBUG_GENERAL_PLASTICITY_LAW) && defined(CHECK_DERIVATIVES)
        if (debug_level > 2)
        {
            KRATOS_WATCH(dfdq)
            KRATOS_WATCH(dqdphi)
            KRATOS_WATCH(l)
            KRATOS_WATCH(dmdq)
            KRATOS_WATCH(dm_dlambda)
        }
        #endif

        // compute left and right hand side
        ninvA.clear();
        SD_MathUtils<double>::ContractFourthOrderTensor(1.0, n, invA, ninvA);
        Ceb.clear();
        SD_MathUtils<double>::ContractFourthOrderTensor(1.0, Ce, b, Ceb);
        double lhs = SD_MathUtils<double>::mat_inner_prod(ninvA, Ceb) - l;
        double rhs = f - SD_MathUtils<double>::mat_inner_prod(ninvA, r);

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            const double num_lhs = PlasticIntegration_ComputeNumLHS(stress, q, alpha,
                    dlambda, 1e-8, stress_trial, Ce, FTOL, max_iters);

            KRATOS_WATCH(dqdphi)
            KRATOS_WATCH(l)
            KRATOS_WATCH(lhs)
            KRATOS_WATCH(num_lhs)
            std::cout << "diff_lhs: " << (lhs-num_lhs) << std::endl;
            std::cout << "diff_lhs/lhs: " << (lhs-num_lhs)/lhs << std::endl;
            // lhs = num_lhs; // try
            KRATOS_WATCH(rhs)
            // ninvAA.clear();
            // SD_MathUtils<double>::ContractFourthOrderTensor(1.0, ninvA, A, ninvAA);
            // KRATOS_WATCH(n)
            // KRATOS_WATCH(ninvAA)
        }
        #endif

        // compute ddlambda
        ddlambda = rhs / lhs;

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(ddlambda)
        }
        #endif

        // update lambda
        dlambda += ddlambda;
        this->dAlphadPhi(dalphadphi, alpha);
        noalias(alpha) += dalphadphi*ddlambda;

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(dlambda)
        }
        #endif

        // update stress
        SD_MathUtils<double>::ContractFourthOrderTensor(ddlambda, Ce, b, r);
        SD_MathUtils<double>::ContractFourthOrderTensor(-1.0, invA, r, stress);

        // update thermodynamics forces
        noalias(q) += dqdphi*ddlambda;

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(stress)
            std::cout << "------------" << std::endl;
        }
        #endif

        ++it;
    }
    while (!converged && it < max_iters);

    #ifdef DEBUG_GENERAL_PLASTICITY_LAW
    if (debug_level > 1)
    {
        std::cout << "PlasticIntegration completed" << std::endl;
        KRATOS_WATCH(it)
        KRATOS_WATCH(f)
        KRATOS_WATCH(stress_trial)
        KRATOS_WATCH(stress)
        KRATOS_WATCH(q)
        KRATOS_WATCH(alpha)
        KRATOS_WATCH(dlambda)
        std::cout << "-------------------" << std::endl;
    }
    #endif

    if (it >= max_iters && !converged)
        return 1;

    return 0;
}

int GeneralElastoplasticityLaw::PlasticIntegration_CuttingPlane(Matrix& stress, Vector& q, Vector& alpha, double& dlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const int debug_level) const
{
    const unsigned int nvars = this->NumberOfInternalVariables();

    Matrix n(3, 3), m(3, 3), r(3, 3), Cem(3, 3);
    Vector dfdq(nvars), dqdphi(nvars), dalphadphi(nvars);
    double ddlambda, l, f;
    const Matrix eye = IdentityMatrix(3);

    #if defined(DEBUG_GENERAL_PLASTICITY_LAW) && defined(CHECK_DERIVATIVES)
    Matrix num_n(3, 3);
    Matrix num_m(3, 3);
    #endif

    bool converged = false;
    int it = 0;
    dlambda = ddlambda = 0.0;
    double norm_r;
    const Vector alphan = alpha;
    do
    {
        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(it)
        }
        #endif

        // compute f
        f = this->F(stress, q);

        // compute gradient m
        this->dGdSigma(m, stress, q);

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
        }
        #endif

        // compute r
        noalias(r) = stress - stress_trial;
        SD_MathUtils<double>::ContractFourthOrderTensor(dlambda, Ce, m, r);
        norm_r = norm_frobenius(r);

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(f)
            KRATOS_WATCH(norm_r)
            KRATOS_WATCH(FTOL)
        }
        #endif

        // check convergence
        if (std::abs(f) + norm_r < FTOL)
        {
            converged = true;
            break;
        }

        // compute gradient n
        this->dFdSigma(n, stress, q);

        #if defined(DEBUG_GENERAL_PLASTICITY_LAW) && defined(CHECK_DERIVATIVES)
        if (debug_level > 2)
        {
            this->Num_dFdSigma(num_n, stress, q, 1e-6);
            // KRATOS_WATCH(n)
            // KRATOS_WATCH(num_n)
            std::cout << "diff_n: " << norm_frobenius(n-num_n) << std::endl;
            this->Num_dGdSigma(num_m, stress, q, 1e-6);
            // KRATOS_WATCH(m)
            // KRATOS_WATCH(num_m)
            std::cout << "diff_m: " << norm_frobenius(m-num_m) << std::endl;
        }
        #endif

        // compute l
        this->dQdPhi(dqdphi, stress, q, alpha);
        this->dFdQ(dfdq, stress, q);
        l = inner_prod(dfdq, dqdphi);

        #if defined(DEBUG_GENERAL_PLASTICITY_LAW) && defined(CHECK_DERIVATIVES)
        if (debug_level > 2)
        {
            KRATOS_WATCH(dfdq)
            KRATOS_WATCH(dqdphi)
            KRATOS_WATCH(l)
        }
        #endif

        // compute left and right hand side
        Cem.clear();
        SD_MathUtils<double>::ContractFourthOrderTensor(1.0, Ce, m, Cem);
        double lhs = SD_MathUtils<double>::mat_inner_prod(n, Cem) - l;
        double rhs = f;

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(dqdphi)
            KRATOS_WATCH(l)
            KRATOS_WATCH(lhs)
            KRATOS_WATCH(rhs)
        }
        #endif

        // compute ddlambda
        ddlambda = rhs / lhs;

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(ddlambda)
        }
        #endif

        // update lambda
        dlambda += ddlambda;
        this->dAlphadPhi(dalphadphi, alpha);
        noalias(alpha) += dalphadphi*ddlambda;

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(dlambda)
        }
        #endif

        // update stress
        noalias(stress) = stress_trial;
        SD_MathUtils<double>::ContractFourthOrderTensor(-dlambda, Ce, m, stress);

        // update thermodynamics forces
        noalias(q) += dqdphi*ddlambda;

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 2)
        {
            KRATOS_WATCH(stress)
            std::cout << "------------" << std::endl;
        }
        #endif

        ++it;
    }
    while (!converged && it < max_iters);

    #ifdef DEBUG_GENERAL_PLASTICITY_LAW
    if (debug_level > 1)
    {
        std::cout << "PlasticIntegration completed" << std::endl;
        KRATOS_WATCH(it)
        KRATOS_WATCH(f)
        KRATOS_WATCH(stress_trial)
        KRATOS_WATCH(stress)
        KRATOS_WATCH(q)
        KRATOS_WATCH(alpha)
        KRATOS_WATCH(dlambda)
        std::cout << "-------------------" << std::endl;
    }
    #endif

    if (it >= max_iters && !converged)
        return 1;

    return 0;
}

int GeneralElastoplasticityLaw::PlasticIntegration_ComputeStress(Matrix& stress, const Vector& q, const Vector& alpha,
        const double dlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters) const
{
    Matrix n(3, 3), m(3, 3), r(3, 3), ninvA(3, 3);
    const Matrix eye = IdentityMatrix(3);

    Fourth_Order_Tensor A, invA, dm_dsigma;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(invA);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dm_dsigma);

    double ddlambda = 0.0;
    unsigned int it = 0;

    do
    {
        // compute gradient m
        this->dGdSigma(m, stress, q);

        // compute r
        noalias(r) = stress - stress_trial;
        SD_MathUtils<double>::ContractFourthOrderTensor(dlambda, Ce, m, r);

        if (norm_frobenius(r) < FTOL)
            break;

        // compute A
        this->d2GdSigma2(dm_dsigma, stress, q);
        SD_MathUtils<double>::ZeroFourthOrderTensor(A); // A = 0
        SD_MathUtils<double>::SpecialProduct1FourthOrderTensor(1.0, eye, eye, A); // A += delta_ik delta_jl
        SD_MathUtils<double>::ProductFourthOrderTensor(dlambda, Ce, dm_dsigma, A); // A += lambda*Ce:dm_dsigma

        // compute inverse of A
        SD_MathUtils<double>::InvertFourthOrderTensor(A, invA);

        // update stress
        SD_MathUtils<double>::ContractFourthOrderTensor(-1.0, invA, r, stress);

        ++it;
    }
    while (it < max_iters);

    if (it >= max_iters)
        return -1;

    return 0;
}

// double GeneralElastoplasticityLaw::PlasticIntegration_ComputeRHS(const Matrix& stress, const Vector& q, const Vector& alpha,
//         const double dlambda,
//         const Matrix& stress_trial, const Fourth_Order_Tensor& Ce) const
// {
//     double f;

//     // compute f
//     f = this->F(stress, q);

//     return f;
// }

double GeneralElastoplasticityLaw::PlasticIntegration_ComputeRHS(const Matrix& stress, const Vector& q, const Vector& alpha,
        const double dlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce) const
{
    double f;
    Matrix n(3, 3), m(3, 3), r(3, 3), ninvA(3, 3);
    const Matrix eye = IdentityMatrix(3);

    Fourth_Order_Tensor A, invA, dm_dsigma;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(invA);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dm_dsigma);

    // compute f
    f = this->F(stress, q);

    // compute gradient m
    this->dGdSigma(m, stress, q);

    // compute r
    noalias(r) = stress - stress_trial;
    SD_MathUtils<double>::ContractFourthOrderTensor(dlambda, Ce, m, r);

    // compute A
    this->d2GdSigma2(dm_dsigma, stress, q);
    SD_MathUtils<double>::ZeroFourthOrderTensor(A); // A = 0
    SD_MathUtils<double>::SpecialProduct1FourthOrderTensor(1.0, eye, eye, A); // A += delta_ik delta_jl
    SD_MathUtils<double>::ProductFourthOrderTensor(dlambda, Ce, dm_dsigma, A); // A += lambda*Ce:dm_dsigma

    // compute inverse of A
    SD_MathUtils<double>::InvertFourthOrderTensor(A, invA);

    // compute gradient n
    this->dFdSigma(n, stress, q);

    // compute right hand side
    ninvA.clear();
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0, n, invA, ninvA);
    double rhs = f - SD_MathUtils<double>::mat_inner_prod(ninvA, r);

    // std::cout << "Within PlasticIntegration_ComputeRHS:" << std::endl;
    // KRATOS_WATCH(stress)
    // KRATOS_WATCH(q)
    // KRATOS_WATCH(alpha)
    // KRATOS_WATCH(f)
    // KRATOS_WATCH(n)
    // KRATOS_WATCH(m)
    // KRATOS_WATCH(r)
    // KRATOS_WATCH(invA)
    // KRATOS_WATCH(ninvA)
    // std::cout << "------------------------" << std::endl;

    return rhs;
}

double GeneralElastoplasticityLaw::PlasticIntegration_ComputeNumLHS(const Matrix& stress, const Vector& q, const Vector& alpha,
        const double dlambda, const double ddlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters) const
{
    const unsigned int nvars = this->NumberOfInternalVariables();

    double rref = PlasticIntegration_ComputeRHS(stress, q, alpha, dlambda, stress_trial, Ce);

    Matrix m(3, 3);
    this->dGdSigma(m, stress, q);
    Matrix new_stress = stress;
    SD_MathUtils<double>::ContractFourthOrderTensor(-ddlambda, Ce, m, new_stress);

    // Matrix new_stress = stress;
    // int error_code = PlasticIntegration_ComputeStress(new_stress, q, alpha,
    //     dlambda + ddlambda, stress_trial, Ce,
    //     FTOL, max_iters);
    // if (error_code != 0)
    //     KRATOS_ERROR << "PlasticIntegration_ComputeStress does not converge";

    Vector new_q = q, new_alpha = alpha;

    Vector dqdphi(nvars), dalphadphi(nvars);
    this->dQdPhi(dqdphi, stress, q, alpha);
    noalias(new_q) += ddlambda*dqdphi;

    this->dAlphadPhi(dalphadphi, alpha);
    noalias(new_alpha) += ddlambda*dalphadphi;

    // Matrix new_stress = stress;
    // int error_code = PlasticIntegration_ComputeStress(new_stress, new_q, new_alpha,
    //     dlambda + ddlambda, stress_trial, Ce,
    //     FTOL, max_iters);
    // if (error_code != 0)
    //     KRATOS_ERROR << "PlasticIntegration_ComputeStress does not converge";

    double rnew = PlasticIntegration_ComputeRHS(new_stress, new_q, new_alpha, dlambda + ddlambda, stress_trial, Ce);

    return (rref - rnew) / ddlambda;
}

/// Because the tangent is evaluated always with reference to the yield state, the substepping scheme
/// does not fit with the DC constitutive law
void GeneralElastoplasticityLaw::ComputeConsistentPlasticTangent_Substepping(Fourth_Order_Tensor& Cep,
        const Fourth_Order_Tensor& Ce, const Matrix& stressn, const Vector& qn, const Vector& alphan,
        const std::vector<double>& loads, const Matrix& incremental_strain,
        const double FTOL, const int max_iters) const
{
    Matrix stress = stressn;
    Matrix stress_trial = stress;
    Vector alpha = alphan, q = qn;

    int error_code;
    const int debug_level = 0;
    double factor, dfactor, factor_last = 0.0, dlambda;

    const unsigned int nvars = this->NumberOfInternalVariables();
    Matrix n(3, 3), m(3, 3), b(3, 3), r(3, 3), dm_dlambda(3, 3), ninvA(3, 3), Ceb(3, 3), invACeb(3, 3), dlambda_dsigmatrial(3, 3);
    Vector dfdq(nvars), dqdphi(nvars), dalphadphi(nvars);
    const Matrix eye = IdentityMatrix(3);

    Third_Order_Tensor dmdq;
    SD_MathUtils<double>::InitializeThirdOrderTensor(dmdq, 3, 3, nvars, false);

    Fourth_Order_Tensor A, invA, dm_dsigma, dsigma_dsigmatrial, Cep_temp;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(invA);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dm_dsigma);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dsigma_dsigmatrial);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Cep_temp);

    SD_MathUtils<double>::ZeroFourthOrderTensor(Cep);
    for (std::size_t i = 0; i < loads.size(); ++i)
    {
        // update the factor
        factor = loads[i];
        dfactor = factor - factor_last;
        factor_last = factor;

        // compute stress (trial)
        SD_MathUtils<double>::ContractFourthOrderTensor(dfactor, Ce, incremental_strain, stress);

        // determine if this is elastic or plastic step
        double f = this->F(stress, q);
        if (f < FTOL)
        {
            // elastic step
            SD_MathUtils<double>::AddFourthOrderTensor(dfactor, Ce, Cep);
        }
        else
        {
            // plastic integration
            noalias(stress_trial) = stress;
            error_code = PlasticIntegration(stress, q, alpha, dlambda, stress_trial, Ce,
                    FTOL, max_iters, debug_level);

            // check convergence
            if (error_code != 0)
                KRATOS_ERROR << "The plastic integration does not converge at load step " << factor;

            if (i == 0)
            {
                /* compute the first tangent */
                ComputeConsistentPlasticTangent(Cep, Ce, stress, q, alpha, dlambda);
                SD_MathUtils<double>::ScaleFourthOrderTensor(Cep, dfactor);
            }
            else
            {
                /* update the tangent */

                // compute gradient m
                this->dGdSigma(m, stress, q);

                // compute gradient n
                this->dFdSigma(n, stress, q);

                // compute A
                this->d2GdSigma2(dm_dsigma, stress, q);
                SD_MathUtils<double>::ZeroFourthOrderTensor(A); // A = 0
                SD_MathUtils<double>::SpecialProduct1FourthOrderTensor(1.0, eye, eye, A); // A += delta_ik delta_jl
                SD_MathUtils<double>::ProductFourthOrderTensor(dlambda, Ce, dm_dsigma, A); // A += lambda*Ce:dm_dsigma

                // compute inverse of A
                SD_MathUtils<double>::InvertFourthOrderTensor(A, invA);

                // compute b
                this->dQdPhi(dqdphi, stress, q, alpha);
                this->d2GdSigmadQ(dmdq, stress, q);
                dm_dlambda.clear();
                SD_MathUtils<double>::ContractThirdOrderTensor(1.0, dmdq, dqdphi, dm_dlambda);
                noalias(b) = m + dlambda*dm_dlambda;

                // compute l
                this->dFdQ(dfdq, stress, q);
                double l = inner_prod(dfdq, dqdphi);

                // compute dlambda_dsigmatrial
                ninvA.clear();
                SD_MathUtils<double>::ContractFourthOrderTensor(1.0, n, invA, ninvA);
                Ceb.clear();
                SD_MathUtils<double>::ContractFourthOrderTensor(1.0, Ce, b, Ceb);
                double lhs = SD_MathUtils<double>::mat_inner_prod(ninvA, Ceb) - l;
                noalias(dlambda_dsigmatrial) = (1.0/lhs) * ninvA;

                // compute dsigma_dsigmatrial
                invACeb.clear();
                SD_MathUtils<double>::ContractFourthOrderTensor(1.0, invA, Ceb, invACeb);
                SD_MathUtils<double>::CopyFourthOrderTensor(invA, dsigma_dsigmatrial);
                SD_MathUtils<double>::OuterProductFourthOrderTensor(-1.0, invACeb, dlambda_dsigmatrial, dsigma_dsigmatrial);

                // update the tangent
                SD_MathUtils<double>::AddFourthOrderTensor(dfactor, Ce, Cep);
                SD_MathUtils<double>::ZeroFourthOrderTensor(Cep_temp);
                SD_MathUtils<double>::ProductFourthOrderTensor(1.0, dsigma_dsigmatrial, Cep, Cep_temp);
                SD_MathUtils<double>::CopyFourthOrderTensor(Cep_temp, Cep);
            }
        }
    }
}

void GeneralElastoplasticityLaw::ComputeNumericalPlasticTangent_Substepping(Fourth_Order_Tensor& Cep,
        const Fourth_Order_Tensor& Ce, const Matrix& stressn, const Vector& qn, const Vector& alphan,
        const std::vector<double>& loads, const Matrix& incremental_strain, const double epsilon,
        const double FTOL, const int max_iters) const
{
    Matrix stress = stressn, new_stress = stressn;
    Vector q = qn, alpha = alphan, new_q = qn, new_alpha = alphan;

    PlasticIntegration_Substepping(stress, q, alpha, loads, incremental_strain,
            Ce, FTOL, max_iters);

    Matrix new_incremental_strain(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(new_incremental_strain) = incremental_strain;
            new_incremental_strain(k, l) += 0.5*epsilon;
            new_incremental_strain(l, k) += 0.5*epsilon;

            noalias(new_stress) = stressn;
            noalias(new_q) = qn;
            noalias(new_alpha) = alphan;

            PlasticIntegration_Substepping(new_stress, new_q, new_alpha, loads, new_incremental_strain,
                    Ce, FTOL, max_iters);

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    Cep[i][j][k][l] = (new_stress(i, j) - stress(i, j)) / epsilon;
                }
            }
        }
    }
}

void GeneralElastoplasticityLaw::ComputeConsistentPlasticTangent(Fourth_Order_Tensor& Cep,
        const Fourth_Order_Tensor& Ce, const Matrix& stress, const Vector& q, const Vector& alpha,
        const double dlambda) const
{
    const unsigned int nvars = this->NumberOfInternalVariables();

    // initialize some variables
    Matrix n(3, 3), m(3, 3), b(3, 3), dm_dlambda(3, 3), ninvA(3, 3), Ceb(3, 3), ddlambla_de(3, 3);
    Vector dfdq(nvars), dqdphi(nvars);
    double l;
    const Matrix eye = IdentityMatrix(3);

    Third_Order_Tensor dmdq;
    SD_MathUtils<double>::InitializeThirdOrderTensor(dmdq, 3, 3, nvars, false);

    Fourth_Order_Tensor A, invA, dm_dsigma, Aux;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(invA);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dm_dsigma);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Aux);

    // compute gradients n, m and dq_dphi
    this->dQdPhi(dqdphi, stress, q, alpha);
    this->dGdSigma(m, stress, q);
    this->dFdSigma(n, stress, q);

    // compute A
    this->d2GdSigma2(dm_dsigma, stress, q);
    SD_MathUtils<double>::ZeroFourthOrderTensor(A); // A = 0
    SD_MathUtils<double>::SpecialProduct1FourthOrderTensor(1.0, eye, eye, A); // A += delta_ik delta_jl
    SD_MathUtils<double>::ProductFourthOrderTensor(dlambda, Ce, dm_dsigma, A); // A += lambda*Ce:dm_dsigma

    // compute inverse of A
    SD_MathUtils<double>::InvertFourthOrderTensor(A, invA);
    // KRATOS_WATCH(invA)

    // compute b
    this->d2GdSigmadQ(dmdq, stress, q);
    dm_dlambda.clear();
    SD_MathUtils<double>::ContractThirdOrderTensor(1.0, dmdq, dqdphi, dm_dlambda);
    noalias(b) = m + dlambda*dm_dlambda;

    // compute l
    this->dFdQ(dfdq, stress, q);
    l = inner_prod(dfdq, dqdphi);
    ninvA.clear();
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0, n, invA, ninvA);
    Ceb.clear();
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0, Ce, b, Ceb);
    double lhs = SD_MathUtils<double>::mat_inner_prod(ninvA, Ceb) - l;
    // KRATOS_WATCH(ninvA)
    // KRATOS_WATCH(Ceb)
    // KRATOS_WATCH(l)
    // KRATOS_WATCH(lhs)

    // compute ddlamda_de
    ddlambla_de.clear();
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0/lhs, ninvA, Ce, ddlambla_de);
    // KRATOS_WATCH(dlambda)
    // KRATOS_WATCH(ddlambla_de)

    // compute dsigma_de
    SD_MathUtils<double>::ZeroFourthOrderTensor(Cep);
    SD_MathUtils<double>::CopyFourthOrderTensor(Ce, Aux);
    SD_MathUtils<double>::OuterProductFourthOrderTensor(-1.0, Ceb, ddlambla_de, Aux);
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0, invA, Aux, Cep);
}

} // Namespace Kratos

#undef CHECK_DERIVATIVES
#undef DEBUG_GENERAL_PLASTICITY_LAW
#undef OSTR
