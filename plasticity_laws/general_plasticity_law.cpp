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
#include "plasticity_laws/general_plasticity_law.h"

// #define DEBUG_GENERAL_PLASTICITY_LAW
#define CHECK_DERIVATIVES
#define OSTR std::cout

namespace Kratos
{

double GeneralPlasticityLaw::ComputeFirstYieldPoint(const Matrix& stress0, const Matrix& delta_stress,
        const Vector& q, const int ndiv, const double NTOL, const double FTOL) const
{
    double dx = 1.0/(double)ndiv;

    int i;
    double xleft, xright, fleft, fright;
    bool found = false;
    for (i = 0; i < ndiv; ++i)
    {
        xleft = i*dx;
        xright = (i+1)*dx;
        fleft = this->F(stress0 + xleft*delta_stress, q);
        fright = this->F(stress0 + xright*delta_stress, q);

        if (fleft*fright < 0.0)
        {
            found = true;
            break;
        }
    }

    if (!found)
    {
        double f0 = this->F(stress0, q);
        if (fabs(f0) < FTOL) return 0.0;

        KRATOS_WATCH("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        // KRATOS_WATCH(mElemId)
        // KRATOS_WATCH(mGaussId)
        KRATOS_WATCH(f0)
        KRATOS_WATCH(FTOL)
        KRATOS_WATCH(this->F(stress0, q))
        KRATOS_WATCH(this->F(stress0 + delta_stress, q))
        KRATOS_ERROR << "The subincrement to identify yield cut is not identifiable";
    }

    double xmid, fmid;
    while(fabs(xleft - xright) > NTOL)
    {
        xmid = 0.5*(xleft + xright);
        fmid = this->F(stress0 + xmid*delta_stress, q);
        if (fleft*fmid > 0.0) xleft = xmid;
        else xright = xmid;
    }

    return xmid;
}

int GeneralPlasticityLaw::DriftCorrection(Matrix& stress, Vector& q, Vector& alpha, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const int debug_level) const
{
    const unsigned int nvars = this->NumberOfInternalVariables();

    Matrix stress_old = stress;
    Vector q_old = q;
    Vector alpha_old = alpha;

    double f = this->F(stress, q);
    #ifdef DEBUG_GENERAL_PLASTICITY_LAW
    if (debug_level > 0)
    {
        OSTR << "  f = " << this->F(stress, q) << ", FTOL = " << FTOL << std::endl;
        OSTR << "  "; this->ReportStressState(OSTR, stress) << std::endl;
    }
    #endif
    int it = 0;
    double dgamma, fnew;
    Matrix A(3, 3), B(3, 3), CexB(3, 3);
    Vector b(nvars), c(nvars), d(nvars);
    while((fabs(f) > FTOL) && (it++ < max_iters))
    {
        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        OSTR << "  drift correction step " << it << ":" << std::endl;
        #endif

        this->dFdSigma(A, stress, q);
        this->dGdSigma(B, stress, q);
        this->dFdQ(c, stress, q);
        this->dGdQ(b, stress, q);
        this->dQdPhi(d, stress, q, alpha);
        CexB.clear();
        SD_MathUtils<double>::ContractFourthOrderTensor(1.0, Ce, B, CexB);
        dgamma = f / (SD_MathUtils<double>::mat_inner_prod(A, CexB) - inner_prod(c, d));

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 1)
        {
            OSTR << "   dgamma: " << dgamma << std::endl;
        }
        #endif

        noalias(stress) -= dgamma*CexB;
        noalias(q) += dgamma*d;
        noalias(alpha) -= dgamma*b;
        fnew = this->F(stress, q);

        if (fabs(fnew) > fabs(f))
        {
            dgamma = f/SD_MathUtils<double>::mat_inner_prod(A, A);
            noalias(stress) = stress_old - dgamma*A;
            noalias(q) = q_old;
            noalias(alpha) = alpha_old;
            fnew = this->F(stress, q);
        }

        f = fnew;
        noalias(stress_old) = stress;
        noalias(q_old) = q;
        noalias(alpha_old) = alpha;

        #ifdef DEBUG_GENERAL_PLASTICITY_LAW
        if (debug_level > 1)
        {
            OSTR << "   f: " << f << ", FTOL = " << FTOL << std::endl;
        }
        #endif
    }

    if ((it >= max_iters) && (fabs(f) > FTOL))
        return 1;

    return 0;
}

{

void GeneralPlasticityLaw::Num_dFdSigma(Matrix& n, const Matrix& stress, const Vector& q, const double epsilon) const
{
    const double F = this->F(stress, q);

    Matrix new_stress(3, 3);
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = i; j < 3; ++j)
        {
            noalias(new_stress) = stress;
            new_stress(i, j) += 0.5*epsilon;
            new_stress(j, i) += 0.5*epsilon;
            double new_F = this->F(new_stress, q);

            n(i, j) = (new_F - F) / epsilon;
            n(j, i) = n(i, j);
        }
    }
}

void GeneralPlasticityLaw::Num_d2FdSigma2(Fourth_Order_Tensor& dn_dsigma, const Matrix& stress, const Vector& q, const double epsilon) const
{
    Matrix n(3, 3), new_n(3, 3);
    this->dFdSigma(n, stress, q);

    Matrix new_stress(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(new_stress) = stress;
            new_stress(k, l) += 0.5*epsilon;
            new_stress(l, k) += 0.5*epsilon;
            this->dFdSigma(new_n, new_stress, q);

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    dn_dsigma[i][j][k][l] = (new_n(i, j) - n(i, j)) / epsilon;
                }
            }
        }
    }
}

void GeneralPlasticityLaw::Num_dGdSigma(Matrix& m, const Matrix& stress, const Vector& q, const double epsilon) const
{
    const double G = this->G(stress, q);

    Matrix new_stress(3, 3);
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = i; j < 3; ++j)
        {
            noalias(new_stress) = stress;
            new_stress(i, j) += 0.5*epsilon;
            new_stress(j, i) += 0.5*epsilon;
            double new_G = this->G(new_stress, q);

            m(i, j) = (new_G - G) / epsilon;
            m(j, i) = m(i, j);
        }
    }
}

void GeneralPlasticityLaw::Num_d2GdSigma2(Fourth_Order_Tensor& dm_dsigma, const Matrix& stress, const Vector& q, const double epsilon) const
{
    Matrix m(3, 3), new_m(3, 3);
    this->dGdSigma(m, stress, q);

    Matrix new_stress(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(new_stress) = stress;
            new_stress(k, l) += 0.5*epsilon;
            new_stress(l, k) += 0.5*epsilon;
            this->dGdSigma(new_m, new_stress, q);

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    dm_dsigma[i][j][k][l] = (new_m(i, j) - m(i, j)) / epsilon;
                }
            }
        }
    }
}

void GeneralPlasticityLaw::Num_dFdQ(Vector& dfdq, const Matrix& stress, const Vector& q, const double epsilon) const
{
    const double F = this->F(stress, q);

    Vector new_q(q.size());
    for (unsigned int i = 0; i < q.size(); ++i)
    {
        noalias(new_q) = q;
        new_q(i) += epsilon;
        double new_F = this->F(stress, new_q);

        dfdq(i) = (new_F - F) / epsilon;
    }
}

void GeneralPlasticityLaw::Num_d2GdSigmadQ(Third_Order_Tensor& dmdq, const Matrix& stress, const Vector& q, const double epsilon) const
{
    Matrix m(3, 3), new_m(3, 3);
    this->dGdSigma(m, stress, q);

    Vector new_q(q.size());
    for (unsigned int k = 0; k < q.size(); ++k)
    {
        noalias(new_q) = q;
        new_q(k) += epsilon;
        this->dGdSigma(new_m, stress, new_q);

        for (unsigned int i = 0; i < 3; ++i)
        {
            for (unsigned int j = 0; j < 3; ++j)
            {
                dmdq[i][j][k] = (new_m(i, j) - m(i, j)) / epsilon;
            }
        }
    }
}

} // Namespace Kratos

#undef CHECK_DERIVATIVES
#undef DEBUG_GENERAL_PLASTICITY_LAW
#undef OSTR
