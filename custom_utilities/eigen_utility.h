//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Jul 2017 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_EIGEN_UTILITY_INCLUDED )
#define  KRATOS_EIGEN_UTILITY_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "utilities/math_utils.h"
#include "custom_utilities/eig/eig.h"
#include "custom_utilities/eig/eig3.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{

class EigenUtility
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( EigenUtility );

    typedef double (*unitary_func_t)(double);

    typedef typename SD_MathUtils<double>::Fourth_Order_Tensor Fourth_Order_Tensor;

    EigenUtility() {}

    virtual ~EigenUtility() {}

    /// Unitary function
    static double exp2(double x) {return std::exp(2.0*x);}
    static double log(double x) {return std::log(x);}
    static double dlog(double x) {return 1.0/x;}
    static double logd2(double x) {return 0.5*std::log(x);}
    static double dlogd2(double x) {return 0.5/x;}
    static double sqrt(double x) {return std::sqrt(x);}

    /*
    * Compute principal stresses and direction using eig3
    * http://barnesc.blogspot.de/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
    * Remarks: sigma_1, sigma_2, sigma_3 is sorted; the eigen vector is sorted accordingly
    */
    static void calculate_principle_stresses(
        double sigma_xx,
        double sigma_yy,
        double sigma_zz,
        double sigma_xy,
        double sigma_yz,
        double sigma_zx,
        double& sigma_1, double& sigma_2, double& sigma_3,
        array_1d<double, 3>& dir1, array_1d<double, 3>& dir2, array_1d<double, 3>& dir3
    )
    {
        double A[3][3];
        double V[3][3];
        double d[3];

        A[0][0] = sigma_xx;
        A[0][1] = sigma_xy;
        A[0][2] = sigma_zx;

        A[1][0] = sigma_xy;
        A[1][1] = sigma_yy;
        A[1][2] = sigma_yz;

        A[2][0] = sigma_zx;
        A[2][1] = sigma_yz;
        A[2][2] = sigma_zz;

        eig3::eigen_decomposition(A, V, d);

        sigma_1 = d[2];
        sigma_2 = d[1];
        sigma_3 = d[0];

        dir1[0] = V[0][2];
        dir1[1] = V[1][2];
        dir1[2] = V[2][2];

        dir2[0] = V[0][1];
        dir2[1] = V[1][1];
        dir2[2] = V[2][1];

        dir3[0] = V[0][0];
        dir3[1] = V[1][0];
        dir3[2] = V[2][0];
    }

    /*
    * Compute principal stresses and direction using eig3
    * http://barnesc.blogspot.de/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
    * Remarks: sigma_1, sigma_2, sigma_3 is sorted; the eigen projection tensor is sorted accordingly
    */
    static void calculate_principle_stresses(
        double sigma_xx,
        double sigma_yy,
        double sigma_zz,
        double sigma_xy,
        double sigma_yz,
        double sigma_zx,
        double& sigma_1, double& sigma_2, double& sigma_3,
        Matrix& eigprj1, Matrix& eigprj2, Matrix& eigprj3
    )
    {
        double A[3][3];
        double V[3][3];
        double d[3];

        A[0][0] = sigma_xx;
        A[0][1] = sigma_xy;
        A[0][2] = sigma_zx;

        A[1][0] = sigma_xy;
        A[1][1] = sigma_yy;
        A[1][2] = sigma_yz;

        A[2][0] = sigma_zx;
        A[2][1] = sigma_yz;
        A[2][2] = sigma_zz;

        eig3::eigen_decomposition(A, V, d);

        sigma_1 = d[2];
        sigma_2 = d[1];
        sigma_3 = d[0];

        if(eigprj1.size1() != 3 || eigprj1.size2() != 3)
            eigprj1.resize(3, 3, false);

        if(eigprj2.size1() != 3 || eigprj2.size2() != 3)
            eigprj2.resize(3, 3, false);

        if(eigprj3.size1() != 3 || eigprj3.size2() != 3)
            eigprj3.resize(3, 3, false);

        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                eigprj1(i, j) = V[i][2] * V[j][2]; // dir1[i] * dir1[j];
                eigprj2(i, j) = V[i][1] * V[j][1]; // dir2[i] * dir2[j];
                eigprj3(i, j) = V[i][0] * V[j][0]; // dir3[i] * dir3[j];
            }
        }
    }

    /*
    * Compute principal stresses and direction using eig3
    * http://barnesc.blogspot.de/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
    * Remarks: sigma_1, sigma_2, sigma_3 is sorted descending except when the reversed flag is true
    */
    static void calculate_principle_stresses(
        double sigma_xx,
        double sigma_yy,
        double sigma_zz,
        double sigma_xy,
        double sigma_yz,
        double sigma_zx,
        double& sigma_1, double& sigma_2, double& sigma_3,
        bool reversed = false
    )
    {
        double A[3][3];
        double V[3][3];
        double d[3];

        A[0][0] = sigma_xx;
        A[0][1] = sigma_xy;
        A[0][2] = sigma_zx;

        A[1][0] = sigma_xy;
        A[1][1] = sigma_yy;
        A[1][2] = sigma_yz;

        A[2][0] = sigma_zx;
        A[2][1] = sigma_yz;
        A[2][2] = sigma_zz;

        eig3::eigen_decomposition(A, V, d);

        if (!reversed)
        {
            sigma_1 = d[2];
            sigma_2 = d[1];
            sigma_3 = d[0];
        }
        else
        {
            sigma_1 = d[0];
            sigma_2 = d[1];
            sigma_3 = d[2];
        }
    }

    /*
    * Compute principal strain and direction using for plane strain case
    * REF: Souza de Neto, Computational Plasticity, box A.2
    * Remarks: sigma_1, sigma_2, sigma_3 is sorted
    */
    static void calculate_principle_strain(
        double epsilon_xx,
        double epsilon_yy,
        double epsilon_xy,
        double& epsilon_1, double& epsilon_2,
        Matrix& eigprj1, Matrix& eigprj2
    )
    {
        // compute the eigenvalues
        double I1 = epsilon_xx + epsilon_yy;
        double I2 = epsilon_xx*epsilon_yy - epsilon_xy*epsilon_xy;
        double x1 = 0.5 * (I1 + sqrt(I1*I1 - 4.0*I2));
        double x2 = 0.5 * (I1 - sqrt(I1*I1 - 4.0*I2));

        // compute the eigenprojection tensor
        if(eigprj1.size1() != 3 || eigprj1.size2() != 3)
            eigprj1.resize(3, 3, false);

        if(eigprj2.size1() != 3 || eigprj2.size2() != 3)
            eigprj2.resize(3, 3, false);

        double diff = fabs(x1 - x2);
        double maxeig = std::max(fabs(x1), fabs(x2));
        if(maxeig != 0.0) diff /= maxeig;
        const double tol = 1.0e-5;
        if(diff < tol)
        {
            // repeated eigenvalues
            Matrix aux(2, 2);
            aux(0, 0) = epsilon_xx;
            aux(1, 1) = epsilon_yy;
            aux(0, 1) = epsilon_xy;
            aux(1, 0) = epsilon_xy;

            Vector eig(2);
            Matrix eigvec(2, 2);
            unsigned int nrot;
            jacobi(aux, eig, eigvec, nrot);
//            KRATOS_WATCH(norm_2(column(eigvec, 0)))

            unsigned int ii, jj;
            if(eig(0) >= eig(1))
            {
                ii = 0;
                jj = 1;
            }
            else
            {
                ii = 1;
                jj = 0;
            }

            x1 = eig(ii);
            x2 = eig(jj);
            noalias(eigprj1) = ZeroMatrix(3, 3);
            noalias(eigprj2) = ZeroMatrix(3, 3);
            for(unsigned int i = 0; i < 2; ++i)
            {
                for(unsigned int j = 0; j < 2; ++j)
                {
                    eigprj1(i, j) = eigvec(i, ii) * eigvec(j, ii);
                    eigprj2(i, j) = eigvec(i, jj) * eigvec(j, jj);
                }
            }
        }
        else
        {
            // non-repeated eigenvalues
            double b = x1 - I1;
            double c = 1.0 / (x1 + b);
            noalias(eigprj1) = ZeroMatrix(3, 3);
            eigprj1(0, 0) = c * (epsilon_xx + b);
            eigprj1(1, 1) = c * (epsilon_yy + b);
            eigprj1(0, 1) = c * epsilon_xy;
            eigprj1(1, 0) = eigprj1(0, 1);
//            KRATOS_WATCH(eigprj1)

            b = x2 - I1;
            c = 1.0 / (x2 + b);
            noalias(eigprj2) = ZeroMatrix(3, 3);
            eigprj2(0, 0) = c * (epsilon_xx + b);
            eigprj2(1, 1) = c * (epsilon_yy + b);
            eigprj2(0, 1) = c * epsilon_xy;
            eigprj2(1, 0) = eigprj2(0, 1);
//            KRATOS_WATCH(eigprj2)
        }
    }

    /*
    * Compute principal stresses and direction using for plane strain case
    * REF: Souza de Neto, Computational Plasticity, box A.2
    * Remarks: sigma_1, sigma_2, sigma_3 is sorted
    */
    static void calculate_principle_stresses(
        double sigma_xx,
        double sigma_yy,
        double sigma_zz,
        double sigma_xy,
        double& sigma_1, double& sigma_2, double& sigma_3,
        Matrix& eigprj1, Matrix& eigprj2, Matrix& eigprj3
    )
    {
        // compute the eigenvalues
        double I1 = sigma_xx + sigma_yy;
        double I2 = sigma_xx*sigma_yy - sigma_xy*sigma_xy;
        double x1 = 0.5 * (I1 + sqrt(I1*I1 - 4.0*I2));
        double x2 = 0.5 * (I1 - sqrt(I1*I1 - 4.0*I2));
        double x3 = sigma_zz;

        // compute the eigenprojection tensor
        if(eigprj1.size1() != 3 || eigprj1.size2() != 3)
            eigprj1.resize(3, 3, false);

        if(eigprj2.size1() != 3 || eigprj2.size2() != 3)
            eigprj2.resize(3, 3, false);

        if(eigprj3.size1() != 3 || eigprj3.size2() != 3)
            eigprj3.resize(3, 3, false);
        noalias(eigprj3) = ZeroMatrix(3, 3);
        eigprj3(2, 2) = 1.0;

        double diff = fabs(x1 - x2);
        double maxeig = std::max(fabs(x1), fabs(x2));
        if(maxeig != 0.0) diff /= maxeig;
        const double tol = 1.0e-5;
        if(diff < tol)
        {
            // repeated eigenvalues
            Matrix aux(2, 2);
            aux(0, 0) = sigma_xx;
            aux(1, 1) = sigma_yy;
            aux(0, 1) = sigma_xy;
            aux(1, 0) = sigma_xy;

            Vector eig(2);
            Matrix eigvec(2, 2);
            unsigned int nrot;
            jacobi(aux, eig, eigvec, nrot);
//            KRATOS_WATCH(norm_2(column(eigvec, 0)))

            unsigned int ii, jj;
            if(eig(0) >= eig(1))
            {
                ii = 0;
                jj = 1;
            }
            else
            {
                ii = 1;
                jj = 0;
            }

            x1 = eig(ii);
            x2 = eig(jj);
            noalias(eigprj1) = ZeroMatrix(3, 3);
            noalias(eigprj2) = ZeroMatrix(3, 3);
            for(unsigned int i = 0; i < 2; ++i)
            {
                for(unsigned int j = 0; j < 2; ++j)
                {
                    eigprj1(i, j) = eigvec(i, ii) * eigvec(j, ii);
                    eigprj2(i, j) = eigvec(i, jj) * eigvec(j, jj);
                }
            }
        }
        else
        {
            // non-repeated eigenvalues
            double b = x1 - I1;
            double c = 1.0 / (x1 + b);
            noalias(eigprj1) = ZeroMatrix(3, 3);
            eigprj1(0, 0) = c * (sigma_xx + b);
            eigprj1(1, 1) = c * (sigma_yy + b);
            eigprj1(0, 1) = c * sigma_xy;
            eigprj1(1, 0) = eigprj1(0, 1);
//            KRATOS_WATCH(eigprj1)

            b = x2 - I1;
            c = 1.0 / (x2 + b);
            noalias(eigprj2) = ZeroMatrix(3, 3);
            eigprj2(0, 0) = c * (sigma_xx + b);
            eigprj2(1, 1) = c * (sigma_yy + b);
            eigprj2(0, 1) = c * sigma_xy;
            eigprj2(1, 0) = eigprj2(0, 1);
//            KRATOS_WATCH(eigprj2)
        }

        // sort the eigenvalues
        if(x1 >= x3)
        {
            sigma_1 = x1;
            if(x2 >= x3)
            {
                sigma_2 = x2;
                sigma_3 = x3;
            }
            else
            {
                sigma_2 = x3;
                sigma_3 = x2;
                Matrix tmp = eigprj2;
                eigprj2 = eigprj3;
                eigprj3 = tmp;
            }
        }
        else
        {
            sigma_1 = x3;
            sigma_2 = x1;
            sigma_3 = x2;

            Matrix tmp1 = eigprj1;
            Matrix tmp2 = eigprj2;
            Matrix tmp3 = eigprj3;
            eigprj1 = tmp3;
            eigprj2 = tmp1;
            eigprj3 = tmp2;
        }
    }

    /*
     * Compute the isotropic function of the type
     *       Y(X) = sum{ y(x_i) E_i }
     * WHERE Y AND X ARE SYMMETRIC TENSORS, x_i AND E_i ARE, RESPECTIVELY
     * THE EIGENVALUES AND EIGENPROJECTIONS OF X, AND y(.) IS A SCALAR
     * FUNCTION.
     * X must be 3 x 3 matrix
     * Y will be resized accordingly
     * Rererence: Section A.5.2, Computational Plasticity, de Souza Neto.
     */
    template<typename TMatrixType>
    static void ComputeIsotropicTensorFunction(
        unitary_func_t func,
        TMatrixType& Y,
        const TMatrixType& X
    )
    {
        TMatrixType V(3, 3);
        std::vector<double> e(3);

        eig::eigen_decomposition<3>(X, V, e);

        std::vector<Matrix> eigprj(3);

        for (unsigned int d = 0; d < 3; ++d)
        {
            eigprj[d].resize(3, 3, false);

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    eigprj[d](i, j) = V(i, d) * V(j, d);
                }
            }
        }

        if ((Y.size1() != 3) || (Y.size2() != 3))
            Y.resize(3, 3, false);

        Y.clear();
        for (unsigned int d = 0; d < 3; ++d)
            noalias(Y) += func(e[d]) * eigprj[d];
    }

    /*
     * Compute the isotropic function of the type
     *       Y(X) = sum{ y(x_i) E_i }
     * WHERE Y AND X ARE SYMMETRIC TENSORS, x_i AND E_i ARE, RESPECTIVELY
     * THE EIGENVALUES AND EIGENPROJECTIONS OF X, AND y(.) IS A SCALAR
     * FUNCTION.
     * X must be 3 x 3 matrix
     * Y will be resized accordingly
     * Rererence: Section A.5.2, Computational Plasticity, de Souza Neto.
     */
    template<typename TMatrixType>
    static void ComputeIsotropicTensorFunction(
        unitary_func_t func,
        TMatrixType& Y,
        const TMatrixType& X,
        double TOL // tolerance to compare the eigenvalues
    )
    {
        TMatrixType V(3, 3);
        std::vector<double> e(3);

        eig::eigen_decomposition<3>(X, V, e);

        std::vector<Matrix> eigprj(3);

        for (unsigned int d = 0; d < 3; ++d)
        {
            eigprj[d].resize(3, 3, false);

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    eigprj[d](i, j) = V(i, d) * V(j, d);
                }
            }
        }

        if ((Y.size1() != 3) || (Y.size2() != 3))
            Y.resize(3, 3, false);

        Y.clear();
        // It is noted that e is sorted ascending
        if (fabs(e[0] - e[1]) > TOL && fabs(e[1] - e[2]) > TOL)
        {
            for (unsigned int d = 0; d < 3; ++d)
                noalias(Y) += func(e[d]) * eigprj[d];
        }
        else
        {
            const Matrix eye = IdentityMatrix(3);

            if (fabs(e[0] - e[1]) < TOL && fabs(e[1] - e[2]) < TOL)
            {
                noalias(Y) = func(e[0]) * eye;
            }
            else
            {
                if (fabs(e[0] - e[1]) < TOL)
                {
                    noalias(Y) = func(e[2]) * eigprj[2]
                               + func(e[0]) * (eye - eigprj[0]);
                }
                else // (fabs(e[1] - e[2]) < TOL)
                {
                    noalias(Y) = func(e[0]) * eigprj[0]
                               + func(e[2]) * (eye - eigprj[2]);
                }
            }
        }
    }

    /***********************************************************************
    * JACOBI ITERATIVE PROCEDURE FOR SPECTRAL DECOMPOSITION OF A
    * N-DIMENSIONAL SYMMETRIC MATRIX
    *
    * REFERENCE: WH Press, SA Teukolsky, WT Vetting & BP Flannery. Numerical
    *            recipes in FORTRAN: The art of scientific computing. 2nd
    *            Edn., Cambridge University Press, 1992.
    ***********************************************************************/
//    static void jacobi(const Matrix& A,
//            std::vector<double>& eigenvalues, std::vector<Vector>& eigenvectors)
//    {
//        const std::size_t n = A.size1();

//        if(eigenvectors.size() != n)
//            eigenvectors.resize(n);

//        Matrix _A = A;
//        Vector B(n);
//        Vector Z(n);

//        for(unsigned int i = 0; i < n; ++i)
//        {
//            eigenvectors[i] = ZeroVector(n);
//            eigenvectors[i](i) = 1.0;

//            B(i) = A(i, i);
//            D(i) = B(i);
//            Z(i) = 0.0;
//        }

//        // TODO
//
//    }

    /***********************************************************************
    * JACOBI ITERATIVE PROCEDURE FOR SPECTRAL DECOMPOSITION OF A
    * N-DIMENSIONAL SYMMETRIC MATRIX
    *
    * REFERENCE: WH Press, SA Teukolsky, WT Vetting & BP Flannery. Numerical
    *            recipes in C++, Cambridge University Press, 1992.
    ***********************************************************************/
    static inline void rot(Matrix& a, double s, double tau, unsigned int i,
        unsigned int j, unsigned int k, unsigned int l)
    {
        double  g,h;

        g=a(i,j);
        h=a(k,l);
        a(i,j)=g-s*(h+g*tau);
        a(k,l)=h+s*(g-h*tau);
    }

    // a: input matrix
    // d: eigenvalues
    // v: eigenvectors (column-wise)
    static void jacobi(Matrix& a, Vector& d, Matrix& v, unsigned int& nrot)
    {
        unsigned int i,j,ip,iq;
        double tresh,theta,tau,t,sm,s,h,g,c;

        unsigned int n = a.size1();
        Vector b(n),z(n);
        for (ip=0;ip<n;ip++) {
            for (iq=0;iq<n;iq++) v(ip,iq)=0.0;
            v(ip,ip)=1.0;
        }
        for (ip=0;ip<n;ip++) {
            b[ip]=d[ip]=a(ip,ip);
            z[ip]=0.0;
        }
        nrot=0;
        for (i=1;i<=50;i++) {
            sm=0.0;
            for (ip=0;ip<n-1;ip++) {
                for (iq=ip+1;iq<n;iq++)
                    sm += fabs(a(ip,iq));
            }
            if (sm == 0.0)
                return;
            if (i < 4)
                tresh=0.2*sm/(n*n);
            else
                tresh=0.0;
            for (ip=0;ip<n-1;ip++) {
                for (iq=ip+1;iq<n;iq++) {
                    g=100.0*fabs(a(ip,iq));
                    if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
                        && (fabs(d[iq])+g) == fabs(d[iq]))
                            a(ip,iq)=0.0;
                    else if (fabs(a(ip,iq)) > tresh) {
                        h=d[iq]-d[ip];
                        if ((fabs(h)+g) == fabs(h))
                            t=(a(ip,iq))/h;
                        else {
                            theta=0.5*h/(a(ip,iq));
                            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                            if (theta < 0.0) t = -t;
                        }
                        c=1.0/sqrt(1+t*t);
                        s=t*c;
                        tau=s/(1.0+c);
                        h=t*a(ip,iq);
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        a(ip,iq)=0.0;
                        for (j=0;j<ip;j++)
                            rot(a,s,tau,j,ip,j,iq);
                        for (j=ip+1;j<iq;j++)
                            rot(a,s,tau,ip,j,j,iq);
                        for (j=iq+1;j<n;j++)
                            rot(a,s,tau,ip,j,iq,j);
                        for (j=0;j<n;j++)
                            rot(v,s,tau,j,ip,j,iq);
                        ++nrot;
                    }
                }
            }
            for (ip=0;ip<n;ip++) {
                b[ip] += z[ip];
                d[ip]=b[ip];
                z[ip]=0.0;
            }
        }
        KRATOS_THROW_ERROR(std::logic_error, "Too many iterations in routine jacobi", "");
    }

};

}

#endif // KRATOS_EIGEN_UTILITY_INCLUDED
