
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

} // end namespace eig3

} // end namespace Kratos
