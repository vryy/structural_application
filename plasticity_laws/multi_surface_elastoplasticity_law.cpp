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
#include "plasticity_laws/multi_surface_elastoplasticity_law.h"

// #define DEBUG_MULTI_SURFACE_PLASTICITY_LAW
#define OSTR std::cout

namespace Kratos
{

int MultiSurfaceElastoplasticityLaw::PlasticIntegration(const std::vector<int>& active_surfaces,
        Matrix& stress, std::vector<Vector>& q, std::vector<Vector>& alpha, std::vector<double>& dlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props,
        const int debug_level) const
{
    const unsigned int nsurfaces = mpPlasticityLaws.size();
    const unsigned int nactive_surfaces = active_surfaces.size();

    Matrix r(3, 3), dm_dlambda(3, 3), ninvA(3, 3), Ceb(3, 3);
    std::vector<Matrix> m(nsurfaces), n(nsurfaces), b(nsurfaces);
    std::vector<Vector> dqdphi(nsurfaces), dfdq(nsurfaces), dalphadphi(nsurfaces);
    Vector ddlambda(nactive_surfaces), l(nsurfaces), f(nsurfaces);
    Matrix lhs(nactive_surfaces, nactive_surfaces), inv_lhs(nactive_surfaces, nactive_surfaces);
    Vector rhs(nactive_surfaces);
    const Matrix eye = IdentityMatrix(3);

    if (dlambda.size() != nsurfaces)
        dlambda.resize(nsurfaces);

    std::vector<Third_Order_Tensor> dmdq(nsurfaces);
    for (unsigned int i = 0; i < nsurfaces; ++i)
    {
        const unsigned int nvars = mpPlasticityLaws[i]->NumberOfInternalVariables();
        SD_MathUtils<double>::InitializeThirdOrderTensor(dmdq[i], 3, 3, nvars, false);
    }

    Fourth_Order_Tensor A, invA, dm_dsigma;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(invA);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dm_dsigma);

    for (unsigned int i = 0; i < nsurfaces; ++i)
    {
        const unsigned int nvars = mpPlasticityLaws[i]->NumberOfInternalVariables();
        dfdq[i].resize(nvars, false);
        dfdq[i].clear();
        dqdphi[i].resize(nvars, false);
        dqdphi[i].clear();
        dalphadphi[i].resize(nvars, false);
        dalphadphi[i].clear();
        n[i].resize(3, 3, false);
        n[i].clear();
        m[i].resize(3, 3, false);
        m[i].clear();
        b[i].resize(3, 3, false);
        b[i].clear();
        dlambda[i] = 0.0;
    }

    bool converged = false;
    int it = 0;
    ddlambda.clear();
    double norm_r, f_sum, f_sum_ref = 1.0;
    do
    {
        // compute f
        f_sum = 0.0;
        for (int ia : active_surfaces)
        {
            f[ia] = mpPlasticityLaws[ia]->F(stress, q[ia], CurrentProcessInfo, props);
            f_sum += std::abs(f[ia]);
        }

        if (it == 0)
        {
            if (f_sum > 1.0) f_sum_ref = f_sum;
        }

        // compute r
        noalias(r) = stress - stress_trial;
        for (int ia : active_surfaces)
        {
            // compute gradient m
            mpPlasticityLaws[ia]->dGdSigma(m[ia], stress, q[ia], CurrentProcessInfo, props);

            // update r
            SD_MathUtils<double>::ContractFourthOrderTensor(dlambda[ia], Ce, m[ia], r);
        }
        norm_r = norm_frobenius(r);

        if (debug_level > 1)
        {
            std::cout << "At step " << it << ":" << std::endl;
            std::cout << " "; KRATOS_WATCH(norm_r)
            std::cout << " "; KRATOS_WATCH(f_sum)
            std::cout << " "; KRATOS_WATCH(f_sum_ref)
        }

        // check convergence
        if (f_sum/f_sum_ref + norm_r < FTOL)
        {
            converged = true;
            break;
        }

        // compute A
        SD_MathUtils<double>::ZeroFourthOrderTensor(A); // A = 0
        SD_MathUtils<double>::SpecialProduct1FourthOrderTensor(1.0, eye, eye, A); // A += delta_ik delta_jl
        for (int ia : active_surfaces)
        {
            mpPlasticityLaws[ia]->d2GdSigma2(dm_dsigma, stress, q[ia], CurrentProcessInfo, props);
            SD_MathUtils<double>::ProductFourthOrderTensor(dlambda[ia], Ce, dm_dsigma, A); // A += lambda*Ce:dm_dsigma
        }

        // compute inverse of A
        SD_MathUtils<double>::InvertFourthOrderTensor(A, invA);

        // compute gradient n
        for (int ia : active_surfaces)
        {
            mpPlasticityLaws[ia]->dFdSigma(n[ia], stress, q[ia], CurrentProcessInfo, props);
        }

        // compute b
        for (int ia : active_surfaces)
        {
            // compute dqdphi
            mpPlasticityLaws[ia]->dQdPhi(dqdphi[ia], stress, q[ia], alpha[ia], CurrentProcessInfo, props);

            // update b
            mpPlasticityLaws[ia]->d2GdSigmadQ(dmdq[ia], stress, q[ia], CurrentProcessInfo, props);
            dm_dlambda.clear();
            SD_MathUtils<double>::ContractThirdOrderTensor(1.0, dmdq[ia], dqdphi[ia], dm_dlambda);
            noalias(b[ia]) = m[ia] + dlambda[ia]*dm_dlambda;
        }

        // compute l
        for (int ia : active_surfaces)
        {
            mpPlasticityLaws[ia]->dFdQ(dfdq[ia], stress, q[ia], CurrentProcessInfo, props);
            l[ia] = inner_prod(dfdq[ia], dqdphi[ia]);
        }

        // compute the local linear system
        for (unsigned int i = 0; i < nactive_surfaces; ++i)
        {
            const int ia = active_surfaces[i];
            ninvA.clear();
            SD_MathUtils<double>::ContractFourthOrderTensor(1.0, n[ia], invA, ninvA);
            for (unsigned int j = 0; j < nactive_surfaces; ++j)
            {
                const int ja = active_surfaces[j];

                Ceb.clear();
                SD_MathUtils<double>::ContractFourthOrderTensor(1.0, Ce, b[ja], Ceb);

                lhs(i, j) = SD_MathUtils<double>::mat_inner_prod(ninvA, Ceb);
            }

            lhs(i, i) -= l[ia];

            rhs(i) = f[ia] - SD_MathUtils<double>::mat_inner_prod(ninvA, r);
        }

        // compute ddlambda
        if (nactive_surfaces > 3)
        {
            SD_MathUtils<double>::Solve(lhs, ddlambda, rhs);
        }
        else if (nactive_surfaces > 1)
        {
            Matrix inv_lhs(lhs.size1(), lhs.size2());
            double det_lhs;
            MathUtils<double>::InvertMatrix(lhs, inv_lhs, det_lhs);
            noalias(ddlambda) = prod(inv_lhs, rhs);
        }
        else
        {
            ddlambda[0] = rhs(0) / lhs(0, 0);
        }

        // update lambda
        for (unsigned int i = 0; i < nactive_surfaces; ++i)
        {
            const int ia = active_surfaces[i];
            dlambda[ia] += ddlambda[i];
        }

        // update alpha
        for (int ia : active_surfaces)
        {
            mpPlasticityLaws[ia]->dAlphadPhi(dalphadphi[ia], alpha[ia], CurrentProcessInfo, props);
            noalias(alpha[ia]) += dalphadphi[ia]*dlambda[ia];
        }

        // update stress
        for (unsigned int i = 0; i < nactive_surfaces; ++i)
        {
            const int ia = active_surfaces[i];
            SD_MathUtils<double>::ContractFourthOrderTensor(ddlambda[i], Ce, b[ia], r);
        }
        SD_MathUtils<double>::ContractFourthOrderTensor(-1.0, invA, r, stress);

        // update thermodynamics forces
        for (unsigned int i = 0; i < nactive_surfaces; ++i)
        {
            const int ia = active_surfaces[i];
            noalias(q[ia]) += ddlambda[i]*dqdphi[ia];
        }

        ++it;
    }
    while (!converged && it < max_iters);

    if (it >= max_iters && !converged)
        return 1;

    return 0;
}

void MultiSurfaceElastoplasticityLaw::ComputeConsistentPlasticTangent(const std::vector<int>& active_surfaces,
        Fourth_Order_Tensor& Cep, const Fourth_Order_Tensor& Ce,
        const Matrix& stress, const std::vector<Vector>& q,
        const std::vector<Vector>& alpha, const std::vector<double>& dlambda,
        const ProcessInfo& CurrentProcessInfo, const Properties& props,
        const int debug_level) const
{
    const unsigned int nsurfaces = mpPlasticityLaws.size();
    const unsigned int nactive_surfaces = active_surfaces.size();

    Matrix m(3, 3), dm_dlambda(3, 3), ninvA(3, 3), Ceb(3, 3);
    std::vector<Matrix> n(nsurfaces), b(nsurfaces), aux(nactive_surfaces), ddlambla_de(nsurfaces);
    std::vector<Vector> dqdphi(nsurfaces), dfdq(nsurfaces), dalphadphi(nsurfaces);
    Vector l(nsurfaces), f(nsurfaces);
    Matrix lhs(nactive_surfaces, nactive_surfaces);
    Vector rhs(nactive_surfaces), sol(nactive_surfaces);
    const Matrix eye = IdentityMatrix(3);

    for (unsigned int i = 0; i < nsurfaces; ++i)
    {
        const unsigned int nvars = mpPlasticityLaws[i]->NumberOfInternalVariables();
        dfdq[i].resize(nvars, false);
        dqdphi[i].resize(nvars, false);
        dalphadphi[i].resize(nvars, false);
        n[i].resize(3, 3, false);
        b[i].resize(3, 3, false);
        ddlambla_de[i].resize(3, 3, false);
    }

    for (unsigned int i = 0; i < nactive_surfaces; ++i)
    {
        aux[i].resize(3, 3, false);
    }

    std::vector<Third_Order_Tensor> dmdq(nsurfaces);
    for (unsigned int i = 0; i < nsurfaces; ++i)
    {
        const unsigned int nvars = mpPlasticityLaws[i]->NumberOfInternalVariables();
        SD_MathUtils<double>::InitializeThirdOrderTensor(dmdq[i], 3, 3, nvars, false);
    }

    Fourth_Order_Tensor A, invA, dm_dsigma, Aux;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(invA);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dm_dsigma);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Aux);

    // compute A
    SD_MathUtils<double>::ZeroFourthOrderTensor(A); // A = 0
    SD_MathUtils<double>::SpecialProduct1FourthOrderTensor(1.0, eye, eye, A); // A += delta_ik delta_jl
    for (int ia : active_surfaces)
    {
        mpPlasticityLaws[ia]->d2GdSigma2(dm_dsigma, stress, q[ia], CurrentProcessInfo, props);
        SD_MathUtils<double>::ProductFourthOrderTensor(dlambda[ia], Ce, dm_dsigma, A); // A += lambda*Ce:dm_dsigma
    }

    // compute inverse of A
    SD_MathUtils<double>::InvertFourthOrderTensor(A, invA);

    // compute dq_dphi
    for (int ia : active_surfaces)
    {
        mpPlasticityLaws[ia]->dQdPhi(dqdphi[ia], stress, q[ia], alpha[ia], CurrentProcessInfo, props);
    }

    // compute gradient n
    for (int ia : active_surfaces)
    {
        mpPlasticityLaws[ia]->dFdSigma(n[ia], stress, q[ia], CurrentProcessInfo, props);
    }

    // compute b
    for (int ia : active_surfaces)
    {
        // compute gradient m
        mpPlasticityLaws[ia]->dGdSigma(m, stress, q[ia], CurrentProcessInfo, props);

        // update b
        mpPlasticityLaws[ia]->d2GdSigmadQ(dmdq[ia], stress, q[ia], CurrentProcessInfo, props);
        dm_dlambda.clear();
        SD_MathUtils<double>::ContractThirdOrderTensor(1.0, dmdq[ia], dqdphi[ia], dm_dlambda);
        noalias(b[ia]) = m + dlambda[ia]*dm_dlambda;
    }

    // compute l
    for (int ia : active_surfaces)
    {
        mpPlasticityLaws[ia]->dFdQ(dfdq[ia], stress, q[ia], CurrentProcessInfo, props);
        l[ia] = inner_prod(dfdq[ia], dqdphi[ia]);
    }

    // compute the local linear system
    for (unsigned int i = 0; i < nactive_surfaces; ++i)
    {
        const int ia = active_surfaces[i];
        ninvA.clear();
        SD_MathUtils<double>::ContractFourthOrderTensor(1.0, n[ia], invA, ninvA);

        for (unsigned int j = 0; j < nactive_surfaces; ++j)
        {
            const int ja = active_surfaces[j];

            Ceb.clear();
            SD_MathUtils<double>::ContractFourthOrderTensor(1.0, Ce, b[ja], Ceb);

            lhs(i, j) = SD_MathUtils<double>::mat_inner_prod(ninvA, Ceb);
        }

        lhs(i, i) -= l[ia];

        aux[i].clear();
        SD_MathUtils<double>::ContractFourthOrderTensor(1.0, ninvA, Ce, aux[i]);
    }

    // compute ddlamda_de
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            for (unsigned int k = 0; k < nactive_surfaces; ++k)
            {
                rhs[k] = aux[k](i, j);
            }

            if (nactive_surfaces > 3)
            {
                SD_MathUtils<double>::Solve(lhs, sol, rhs);
            }
            else if (nactive_surfaces > 1)
            {
                Matrix inv_lhs(lhs.size1(), lhs.size2());
                double det_lhs;
                MathUtils<double>::InvertMatrix(lhs, inv_lhs, det_lhs);
                noalias(sol) = prod(inv_lhs, rhs);
            }
            else
            {
                sol[0] = rhs(0) / lhs(0, 0);
            }

            for (unsigned int k = 0; k < nactive_surfaces; ++k)
            {
                const int ia = active_surfaces[k];
                ddlambla_de[ia](i, j) = sol[k];
            }
        }
    }

    // compute dsigma_de
    SD_MathUtils<double>::ZeroFourthOrderTensor(Cep);
    SD_MathUtils<double>::CopyFourthOrderTensor(Ce, Aux);
    for (int ia : active_surfaces)
    {
        Ceb.clear();
        SD_MathUtils<double>::ContractFourthOrderTensor(1.0, Ce, b[ia], Ceb);
        SD_MathUtils<double>::OuterProductFourthOrderTensor(-1.0, Ceb, ddlambla_de[ia], Aux);
    }
    SD_MathUtils<double>::ProductFourthOrderTensorTN(1.0, invA, Aux, Cep);
}

void MultiSurfaceElastoplasticityLaw::PlasticIntegration_ComputeRHS(Vector& rhs, const std::vector<int>& active_surfaces,
        const Matrix& stress, const std::vector<Vector>& q, const std::vector<Vector>& alpha,
        const std::vector<double>& dlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const
{
    const unsigned int nsurfaces = mpPlasticityLaws.size();
    const unsigned int nactive_surfaces = active_surfaces.size();

    if (rhs.size() != nactive_surfaces)
        rhs.resize(nactive_surfaces, false);

    Matrix m(3, 3), r(3, 3), ninvA(3, 3);
    std::vector<Matrix> n(nsurfaces);
    Vector f(nsurfaces);
    const Matrix eye = IdentityMatrix(3);

    Fourth_Order_Tensor A, invA, dm_dsigma;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(invA);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dm_dsigma);

    for (unsigned int i = 0; i < nsurfaces; ++i)
    {
        n[i].resize(3, 3, false);
    }

    // compute f
    for (int ia : active_surfaces)
    {
        f[ia] = mpPlasticityLaws[ia]->F(stress, q[ia], CurrentProcessInfo, props);
    }

    // compute r
    noalias(r) = stress - stress_trial;
    for (int ia : active_surfaces)
    {
        // compute gradient m
        mpPlasticityLaws[ia]->dGdSigma(m, stress, q[ia], CurrentProcessInfo, props);

        // update r
        SD_MathUtils<double>::ContractFourthOrderTensor(dlambda[ia], Ce, m, r);
    }

    // compute A
    SD_MathUtils<double>::ZeroFourthOrderTensor(A); // A = 0
    SD_MathUtils<double>::SpecialProduct1FourthOrderTensor(1.0, eye, eye, A); // A += delta_ik delta_jl
    for (int ia : active_surfaces)
    {
        mpPlasticityLaws[ia]->d2GdSigma2(dm_dsigma, stress, q[ia], CurrentProcessInfo, props);
        SD_MathUtils<double>::ProductFourthOrderTensor(dlambda[ia], Ce, dm_dsigma, A); // A += lambda*Ce:dm_dsigma
    }

    // compute inverse of A
    SD_MathUtils<double>::InvertFourthOrderTensor(A, invA);
    /////
    // Matrix T(6, 6), InvT(6, 6);
    // SD_MathUtils<double>::TensorToMatrix(A, T);
    // SD_MathUtils<double>::InvertMatrix(T, InvT);
    // SD_MathUtils<double>::MatrixToTensor(InvT, invA);
    /////

    // compute gradient n
    for (int ia : active_surfaces)
    {
        mpPlasticityLaws[ia]->dFdSigma(n[ia], stress, q[ia], CurrentProcessInfo, props);
    }

    // compute the local linear system
    for (unsigned int i = 0; i < nactive_surfaces; ++i)
    {
        const int ia = active_surfaces[i];
        ninvA.clear();
        SD_MathUtils<double>::ContractFourthOrderTensor(1.0, n[ia], invA, ninvA);
        rhs(i) = f[ia] - SD_MathUtils<double>::mat_inner_prod(ninvA, r);
    }
}

void MultiSurfaceElastoplasticityLaw::PlasticIntegration_ComputeNumLHS(Matrix& lhs, const std::vector<int>& active_surfaces,
        const Matrix& stress, const std::vector<Vector>& q, const std::vector<Vector>& alpha,
        const std::vector<double> dlambda, const double ddlambda,
        const Matrix& stress_trial, const Fourth_Order_Tensor& Ce,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const
{
    const unsigned int nsurfaces = mpPlasticityLaws.size();
    const unsigned int nactive_surfaces = active_surfaces.size();

    if (lhs.size1() != nactive_surfaces || lhs.size2() != nactive_surfaces)
        lhs.resize(nactive_surfaces, nactive_surfaces, false);

    Vector ref_rhs(nactive_surfaces), new_rhs(nactive_surfaces);
    PlasticIntegration_ComputeRHS(ref_rhs, active_surfaces, stress, q, alpha, dlambda, stress_trial, Ce, CurrentProcessInfo, props);

    Matrix m(3, 3), new_stress(3, 3);
    std::vector<Vector> new_q(nsurfaces), new_alpha(nsurfaces);
    std::vector<double> new_dlambda(nsurfaces);
    for (unsigned int i = 0; i < nactive_surfaces; ++i)
    {
        int ia = active_surfaces[i];
        const unsigned int nvars = mpPlasticityLaws[ia]->NumberOfInternalVariables();

        new_q[ia].resize(nvars, false);
        new_alpha[ia].resize(nvars, false);
    }

    for (unsigned int i = 0; i < nactive_surfaces; ++i)
    {
        int ia = active_surfaces[i];

        noalias(new_stress) = stress;
        for (unsigned int j = 0; j < nsurfaces; ++j)
        {
            noalias(new_q[j]) = q[j];
            noalias(new_alpha[j]) = alpha[j];
            new_dlambda[j] = dlambda[j];
        }
        new_dlambda[ia] += ddlambda;

        for (unsigned int j = 0; j < nactive_surfaces; ++j)
        {
            int ja = active_surfaces[j];

            const unsigned int nvars = mpPlasticityLaws[ja]->NumberOfInternalVariables();

            mpPlasticityLaws[ja]->dGdSigma(m, stress, q[ja], CurrentProcessInfo, props);
            SD_MathUtils<double>::ContractFourthOrderTensor(dlambda[ja] - new_dlambda[ja], Ce, m, new_stress);

            Vector dqdphi(nvars);
            mpPlasticityLaws[ja]->dQdPhi(dqdphi, stress, q[ja], alpha[ja], CurrentProcessInfo, props);
            noalias(new_q[ja]) += (new_dlambda[ja] - dlambda[ja])*dqdphi;

            Vector dalphadphi(nvars);
            mpPlasticityLaws[ja]->dAlphadPhi(dalphadphi, alpha[ja], CurrentProcessInfo, props);
            noalias(new_alpha[ja]) += (new_dlambda[ja] - dlambda[ja])*dalphadphi;
        }

        PlasticIntegration_ComputeRHS(new_rhs, active_surfaces, new_stress, new_q, new_alpha, new_dlambda, stress_trial, Ce, CurrentProcessInfo, props);

        noalias(column(lhs, i)) = (ref_rhs - new_rhs) / ddlambda;
    }
}

} // Namespace Kratos

#undef DEBUG_MULTI_SURFACE_PLASTICITY_LAW
#undef OSTR
