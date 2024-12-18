
/* Eigen decomposition code for symmetric 3x3 matrices, adapted from the public
   domain Java Matrix library JAMA. */

#include "eig3.h"

namespace Kratos
{

namespace eig3
{

void eigen_decomposition(double A[3][3], double V[3][3], double d[3]) {
  double e[3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      V[i][j] = A[i][j];
    }
  }
  eig::tred2<3>(V, d, e);
  eig::tql2<3>(V, d, e);
}

void spectral_decomposition(
        double sigma_xx,
        double sigma_yy,
        double sigma_zz,
        double sigma_xy,
        double sigma_yz,
        double sigma_zx,
        double& sigma_1, double& sigma_2, double& sigma_3
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

    eigen_decomposition(A, V, d);

    sigma_1 = d[2];
    sigma_2 = d[1];
    sigma_3 = d[0];
}

} // end namespace eig3

} // end namespace Kratos
