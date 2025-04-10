
/* Eigen-decomposition for symmetric 3x3 real matrices.
   Public domain, adapted from the public domain Java library JAMA. */

#ifndef _eig3_h
#define _eig3_h

#include "includes/kratos_export_api.h"

#include "eig.h"

namespace Kratos
{

namespace eig3
{

/* Symmetric matrix A => eigenvectors in columns of V, corresponding
   eigenvalues in d. */
void KRATOS_API(STRUCTURAL_APPLICATION) eigen_decomposition(double A[3][3], double V[3][3], double d[3]);

/*
 * Compute principal stresses and direction using eig3
 * http://barnesc.blogspot.de/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
 * Remarks: sigma_1, sigma_2, sigma_3 is sorted
 */
template<typename TVectorType> // = array_1d<double, 3>
void spectral_decomposition(
                           double sigma_xx,
                           double sigma_yy,
                           double sigma_zz,
                           double sigma_xy,
                           double sigma_yz,
                           double sigma_zx,
                           double& sigma_1, double& sigma_2, double& sigma_3,
                           TVectorType& dir1, TVectorType& dir2, TVectorType& dir3)
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

   eigen_decomposition(A, V, d);

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

void KRATOS_API(STRUCTURAL_APPLICATION) spectral_decomposition(
                           double sigma_xx,
                           double sigma_yy,
                           double sigma_zz,
                           double sigma_xy,
                           double sigma_yz,
                           double sigma_zx,
                           double& sigma_1, double& sigma_2, double& sigma_3);

} // end namespace eig3

} // end namespace Kratos

#endif
