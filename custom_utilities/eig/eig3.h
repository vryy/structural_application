
/* Eigen-decomposition for symmetric 3x3 real matrices.
   Public domain, adapted from the public domain Java library JAMA. */

#ifndef _eig3_h
#define _eig3_h

#include <cmath>
#include "eig.h"

namespace Kratos
{

namespace eig3
{

/* Symmetric matrix A => eigenvectors in columns of V, corresponding
   eigenvalues in d. */
void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);

} // end namespace eig3

} // end namespace Kratos

#endif
