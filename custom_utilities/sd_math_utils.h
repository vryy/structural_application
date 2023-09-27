/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

#if !defined(SD_MATH_UTILS)
#define SD_MATH_UTILS

#include <cmath>
#include "utilities/math_utils.h"
#include "geometries/point.h"

namespace Kratos
{
template<class TDataType> class SD_MathUtils
{
public:
    /**
     * @name type definitions
     * @{
     */
    typedef Matrix MatrixType;

    typedef SymmetricMatrix SymmetricMatrixType;

    typedef Vector VectorType;

    typedef unsigned int IndexType;

    typedef unsigned int SizeType;

    typedef MathUtils<TDataType> MathUtilsType;

    typedef boost::numeric::ublas::vector<VectorType> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector

    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;

    typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<MatrixType> > Fourth_Order_Tensor;

    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matriz.

    typedef double (*unitary_func_t)(double);

    #ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Point PointType;
    #else
    typedef Point<3> PointType;
    #endif

    /// Return constant Pi
    static constexpr TDataType Pi()
    {
        #if defined(__clang__)
        return 3.1415926535897932384626433832795028841971693;
        #elif defined(__GNUC__) || defined(__GNUG__)
        return std::atan(1.0)*4;
        #else
        return std::atan(1.0)*4;
        #endif
    }

    /**
     * @}
     */
    /**
     * calculates the solutions for a given cubic polynomial equation
     * 0= a*x^3+b*x^2+c*x+d
     * @param a coefficient
     * @param b coefficient
     * @param c coefficient
     * @param d coefficient
     * @param ZeroTol number treated as zero
     * @return Vector of solutions
     * WARNING only valid cubic (not quadratic, not linear, not constant) equations with
     * three real (not complex) solutions
     */

    static inline bool CardanoFormula(double a, double b, double c, double d, VectorType& solution)
    {
        solution.resize(3,false);
        noalias(solution)= ZeroVector(3);

        if(a==0)
        {
            std::cout<<"This is not a cubic equation: CardanoFormula"<<std::endl;

            return false;
        }

        double p= (3.0*a*c-b*b)/(3.0*a*a);

        double q= 2.0*b*b*b/(27.0*a*a*a)-b*c/(3.0*a*a)+d/a;

        double discriminante= p*p*p/27.0+q*q/4.0;

        if(discriminante>0)
        {
            return false;
        }

        if(discriminante==0)
        {
            if( a == 0 )
                return false;

            solution(0)= pow(q/2.0, 1.0/3.0)-b/(3*a);
            solution(1)= pow(q/2.0, 1.0/3.0)-b/(3*a);
            solution(2)= pow(-4.0*q, 1.0/3.0)-b/(3*a);

            return true;
        }

        solution(0)=
            -sqrt(-4.0/3.0*p)*cos(1.0/3.0*acos(-q/2.0*sqrt(-27.0/(p*p*p)))+Pi()/3.0)
            -b/(3*a);
        solution(1)=
            sqrt(-4.0/3.0*p)*cos(1.0/3.0*acos(-q/2.0*sqrt(-27.0/(p*p*p))))-b/(3*a)
            ;
        solution(2)=
            -sqrt(-4.0/3.0*p)*cos(1.0/3.0*acos(-q/2.0*sqrt(-27.0/(p*p*p)))-Pi()/3.0)
            -b/(3*a);

#ifdef _DEBUG
        if(std::isnan<double>(solution(0)) || std::isnan<double>(solution(1))|| std::isnan<double>(solution(2)))
        {
            return false;
        }
#endif

        return true;
    }
    /**
     * @}
      */
    /**
     * calculates Eigenvalues of given square matrix A.
     * The QR Algorithm with shifts is used
     * @param A the given square matrix the eigenvalues are to be calculated.
     * @param crit convergence criteria
     * @param zero number treated as zero
     * @return Vector of eigenvalues
     * WARNING only valid for 2*2 and 3*3 Matrices yet
     */

    static inline VectorType EigenValues(const MatrixType& A, double crit, double zero)
    {
        int dim= A.size1();


        MatrixType Convergence(2,dim);

        double delta;

        double abs;

        VectorType Result=ZeroVector(dim);

        MatrixType HelpA= ZeroMatrix(dim, dim);

        MatrixType HelpQ= ZeroMatrix(dim, dim);

        MatrixType HelpR= ZeroMatrix(dim, dim);

        HelpA=A;

        bool is_converged=false;

        while(!(is_converged))
        {
            double shift= HelpA((dim-1),(dim-1));
            //
            for(int i=0; i<dim; i++)
            {
                HelpA(i,i) = HelpA(i,i)- shift;
            }

            QRFactorization(HelpA, HelpQ, HelpR);

            HelpA= ZeroMatrix(dim, dim);

            for(int i=0; i<dim; i++)
            {
                HelpA(i,i) += shift;
                for(int j=0; j< dim; j++)
                {
                    for(int k=0; k< dim; k++)
                    {
                        HelpA(i,j) += HelpR(i,k)*HelpQ(k,j);
                    }
                }
            }

            delta= 0.0;

            abs = 0.0;

            for(int i=0; i<dim; i++)
            {
                Convergence(0,i)=Convergence(1,i);
                Convergence(1,i)=HelpA(i,i);
                delta+= (Convergence(1,i)-Convergence(0,i))*(Convergence(1,i)-Convergence(0,i));
                abs+=(Convergence(1,i))*(Convergence(1,i));
            }

            delta= sqrt(delta);

            abs=sqrt(abs);

            if(abs< zero)
                abs=1.0;

            if(delta < zero || (delta/abs) < crit)
                is_converged=true;

        }

        for(int i=0; i<dim; i++)
        {
            Result(i)= HelpA(i,i);

            if(fabs(Result(i)) <zero)
                Result(i)=0.0;
        }

        return Result;
    }

    // Given a real symmetric 3x3 matrix A, compute the eigenvalues
    // Note that acos and cos operate on angles in radian
    // REF: https://en.wikipedia.org/wiki/Eigenvalue_algorithm
    static inline VectorType EigenValuesSym3x3(const MatrixType& A)
    {
        VectorType Result(3);

        double p1 = pow(A(0, 1), 2) + pow(A(1, 2), 2) + pow(A(0, 2), 2);
        if (p1 == 0)
        {
            // A is diagonal.
            Result[0] = A(0, 0);
            Result[1] = A(1, 1);
            Result[2] = A(2, 2);
        }
        else
        {
            double q = (A(0, 0) + A(1, 1) + A(2, 2))/3; // trace(A) is the sum of all diagonal values
            double p2 = pow(A(0, 0) - q, 2) + pow(A(1, 1) - q, 2) + pow(A(2, 2) - q, 2) + 2 * p1;
            double p = sqrt(p2 / 6);
            MatrixType B = (1 / p) * (A - q * IdentityMatrix(3));
            double r = MathUtils<TDataType>::Det3(B) / 2;

            // In exact arithmetic for a symmetric matrix  -1 <= r <= 1
            // but computation error can leave it slightly outside this range.
            double phi;
            if (r <= -1)
            {
                phi = Pi() / 3;
            }
            else if (r >= 1)
            {
                phi = 0;
            }
            else
                phi = acos(r) / 3;

           // the eigenvalues satisfy eig3 <= eig2 <= eig1
           Result[0] = q + 2 * p * cos(phi);
           Result[2] = q + 2 * p * cos(phi + (2*Pi()/3));
           Result[1] = 3 * q - Result[0] - Result[2];     // since trace(A) = eig1 + eig2 + eig3
        }

        return Result;
    }

    /**
    Organize a 3D eigenvalues vector to the sequence 1 > 2 > 3
     */
    template<class TValuesContainerType>
    static inline VectorType OrganizeEigenvalues(const TValuesContainerType& rEigenvalues)
    {
        VectorType Result(3);

        if(rEigenvalues[0] >= rEigenvalues[1] && rEigenvalues[0] >= rEigenvalues[2])
        {
            Result[0] = rEigenvalues[0];
            if(rEigenvalues[1] >= rEigenvalues[2])
            {
                Result[1] = rEigenvalues[1];
                Result[2] = rEigenvalues[2];
            }
            else
            {
                Result[1] = rEigenvalues[2];
                Result[2] = rEigenvalues[1];
            }
        }
        else if(rEigenvalues[1] >= rEigenvalues[0] && rEigenvalues[1] >= rEigenvalues[2])
        {
            Result[0] = rEigenvalues[1];
            if(rEigenvalues[0] >= rEigenvalues[2])
            {
                Result[1] = rEigenvalues[0];
                Result[2] = rEigenvalues[2];
            }
            else
            {
                Result[1] = rEigenvalues[2];
                Result[2] = rEigenvalues[0];
            }
        }
        else if(rEigenvalues[2] >= rEigenvalues[0] && rEigenvalues[2] >= rEigenvalues[1])
        {
            Result[0] = rEigenvalues[2];
            if(rEigenvalues[0] >= rEigenvalues[1])
            {
                Result[1] = rEigenvalues[0];
                Result[2] = rEigenvalues[1];
            }
            else
            {
                Result[1] = rEigenvalues[1];
                Result[2] = rEigenvalues[0];
            }
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Something must be wrong. This case can't happe. At", __FUNCTION__)

        return Result;
    }


    /**
     * @}
     */
    /**
     * calculates the QR Factorization of given square matrix A=QR.
     * The Factorization is performed using the householder algorithm
     * @param A the given square matrix the factorization is to be calculated.
     * @param Q the result matrix Q
     * @param R the result matrix R
     */

    static inline void QRFactorization(const MatrixType& A, MatrixType& Q, MatrixType& R)
    {

        //QR Factorization with Householder-Algo
        int dim= A.size1();

        VectorType y(dim);

        VectorType w(dim);

        R.resize(dim,dim,false);

        R=ZeroMatrix(dim,dim);

        Q.resize(dim,dim,false);

        Q=ZeroMatrix(dim,dim);

        MatrixType Help= A;

        MatrixType unity= ZeroMatrix(dim,dim);

        for(int j=0; j<dim; j++)
            unity(j,j)=1.0;

        std::vector<MatrixType> HelpQ(dim-1);

        std::vector<MatrixType> HelpR(dim-1);

        for(int i=0; i< dim-1; i++)
        {
            HelpQ[i].resize(dim,dim,false);
            HelpR[i].resize(dim,dim,false);
            noalias(HelpQ[i])= unity;
            noalias(HelpR[i])= ZeroMatrix(dim,dim);
        }

        for(int iteration=0; iteration< dim-1; iteration++)
        {
            //Vector y
            for(int i=iteration; i<dim; i++)
                y(i)= Help(i,iteration);


            //Helpvalue l
            double normy=0.0;

            for(int i=iteration; i<dim; i++)
                normy += y(i)*y(i);

            normy= sqrt(normy);

            double l= sqrt((normy*(normy+fabs(y(iteration))))/2);

            double k=0.0;

            if(y[iteration] !=0)
                k= - y(iteration)/fabs(y(iteration))*normy;
            else
                k= -normy;

            for(int i=iteration; i<dim; i++)
            {
                double e=0;

                if(i==iteration)
                    e=1;

                w(i)= 1/(2*l)*(y(i)-k*e);
            }

            for(int i=iteration; i<dim; i++)
                for(int j=iteration; j<dim; j++)
                    HelpQ[iteration](i,j)= unity(i,j)- 2*w(i)*w(j);


            for(int i=iteration; i<dim; i++)
                for(int j=iteration; j<dim; j++)
                    for(int k=iteration; k<dim; k++)
                        HelpR[iteration](i,j)+= HelpQ[iteration](i,k)*Help(k,j);

            Help= HelpR[iteration];

        }

        //Assembling R
        for(int k=0; k<dim-1; k++)
        {
            for(int i=k; i<dim; i++)
                for(int j=k; j<dim; j++)
                    R(i,j) =HelpR[k](i,j);

        }


        for(int k=1; k<dim-1; k++)
        {
            for(int i=0; i<dim; i++)
                for(int j=0; j<dim; j++)
                    for(int l=0; l<dim; l++)
                        Q(i,j)+= HelpQ[(k-1)](i,l)*HelpQ[k](l,j);
            noalias(HelpQ[k])=Q;
        }
        if(dim-1==1)
            noalias(Q)=HelpQ[0];

    }

    /**
     * @}
     */
    /**
     * calculates the eigenvectors and eigenvalues of given symmetric matrix A.
     * The eigenvectors and eigenvalues are calculated using the iterative
     * Gauss-Seidel-method
     * @param A the given symmetric matrix the eigenvectors are to be calculated.
     * :WARNING: Matrix A will be overwritten and has to be symmetric
     * @param V the result matrix (will be overwritten with the eigenvectors)
     * @param zero_tolerance the largest value considered to be zero
     */


    static inline void EigenVectors(const MatrixType& A, MatrixType& vectors, VectorType& lambda, double zero_tolerance =1e-9, int max_iterations = 10)
    {
        MatrixType Help= A;

        for(int i=0; i<3; i++)
            for(int j=0; j<3; j++)
                Help(i,j)= Help(i,j);


        vectors.resize(Help.size1(),Help.size2(),false);

        lambda.resize(Help.size1(),false);

        MatrixType HelpDummy(Help.size1(),Help.size2());

        bool is_converged = false;

        MatrixType unity=ZeroMatrix(Help.size1(),Help.size2());

        for(unsigned int i=0; i< Help.size1(); i++)
            unity(i,i)= 1.0;

        MatrixType V= unity;

        MatrixType VDummy(Help.size1(),Help.size2());

        MatrixType Rotation(Help.size1(),Help.size2());


        for(int iterations=0; iterations<max_iterations; iterations++)
        {

            is_converged= true;

            double a= 0.0;

            unsigned int index1= 0;

            unsigned int index2= 1;

            for(unsigned int i=0; i< Help.size1(); i++)
            {
                for(unsigned int j=(i+1); j< Help.size2(); j++)
                {
                    if((fabs(Help(i,j)) > a ) && (fabs(Help(i,j)) > zero_tolerance))
                    {
                        a= fabs(Help(i,j));

                        index1= i;
                        index2= j;

                        is_converged= false;
                    }
                }
            }

//                 KRATOS_WATCH(Help);

            if(is_converged)
                break;

            //Calculation of Rotationangle

            double gamma= (Help(index2,index2)-Help(index1,index1))/(2*Help(index1,index2));

            double u=1.0;

            if(fabs(gamma) > zero_tolerance && fabs(gamma)< (1/zero_tolerance))
            {
                u= gamma/fabs(gamma)*1.0/(fabs(gamma)+sqrt(1.0+gamma*gamma));
            }
            else
            {
                if  (fabs(gamma)>= (1.0/zero_tolerance))
                    u= 0.5/gamma;
            }

            double c= 1.0/(sqrt(1.0+u*u));

            double s= c*u;

            double teta= s/(1.0+c);

            //Ratotion of the Matrix
            HelpDummy= Help;

            HelpDummy(index2,index2)= Help(index2,index2)+u*Help(index1,index2);
            HelpDummy(index1,index1)= Help(index1,index1)-u*Help(index1,index2);
            HelpDummy(index1,index2)= 0.0;
            HelpDummy(index2,index1)= 0.0;

            for(unsigned int i=0; i<Help.size1(); i++)
            {
                if((i!= index1) && (i!= index2))
                {
                    HelpDummy(index2,i)=Help(index2,i)+s*(Help(index1,i)- teta*Help(index2,i));
                    HelpDummy(i,index2)=Help(index2,i)+s*(Help(index1,i)- teta*Help(index2,i));

                    HelpDummy(index1,i)=Help(index1,i)-s*(Help(index2,i)+ teta*Help(index1,i));
                    HelpDummy(i,index1)=Help(index1,i)-s*(Help(index2,i)+ teta*Help(index1,i));
                }
            }


            Help= HelpDummy;

            //Calculation of the eigenvectors V
            Rotation =unity;
            Rotation(index2,index1)=-s;
            Rotation(index1,index2)=s;
            Rotation(index1,index1)=c;
            Rotation(index2,index2)=c;

//                 Help=ZeroMatrix(A.size1(),A.size1());

            VDummy = ZeroMatrix(Help.size1(), Help.size2());

            for(unsigned int i=0; i< Help.size1(); i++)
            {
                for(unsigned int j=0; j< Help.size1(); j++)
                {
                    for(unsigned int k=0; k< Help.size1(); k++)
                    {
                        VDummy(i,j) += V(i,k)*Rotation(k,j);
                    }
                }
            }
            V= VDummy;

//                 MatrixType VTA= ZeroMatrix(3,3);
//                 for(int i=0; i< Help.size1(); i++)
//                 {
//                     for(int j=0; j< Help.size1(); j++)
//                     {
//                         for(int k=0; k< Help.size1(); k++)
//                         {
//                             VTA(i,j) += V(k,i)*A(k,j);
//                         }
//                     }
//                 }
//
//                 for(int i=0; i< Help.size1(); i++)
//                 {
//                     for(int j=0; j< Help.size1(); j++)
//                     {
//                         for(int k=0; k< Help.size1(); k++)
//                         {
//                             Help(i,j) += VTA(i,k)*V(k,j);
//                         }
//                     }
//                 }

        }

        if(!(is_converged))
        {
            std::cout<<"########################################################"<<std::endl;
            std::cout<<"Max_Iterations exceed in Jacobi-Seidel-Iteration (eigenvectors)"<<std::endl;
            std::cout<<"########################################################"<<std::endl;
        }

        for(unsigned int i=0; i< Help.size1(); i++)
        {
            for(unsigned int j=0; j< Help.size1(); j++)
            {
                vectors(i,j)= V(j,i);
            }
        }

        for(unsigned int i=0; i<Help.size1(); i++)
            lambda(i)= Help(i,i);

        return;
    }

    /**
     * @}
     */
    /**
     * calculates the eigenvectors and eigenvalues of given matrix A.
     * The eigenvectors and eigenvalues are calculated using the iterative
     * JACOBI-method
     * @param A the given matrix the eigenvectors are to be calculated.
     * :WARNING: Matrix A will be overwritten
     * @param V the result matrix (will be overwritten with the eigenvectors)
     * @param error_tolerance the desired accuracy for the convergence check
     * @param zero_tolerance the largest value considered to be zero
     */
    static inline void EigenVectors( MatrixType& A,
                                     MatrixType& V,
                                     TDataType& error_tolerance,
                                     TDataType zero_tolerance)
    {
        //initial error
        TDataType error = 1.0;
        int n = A.size2();
        //setting V to identity matrix
        V = IdentityMatrix( V.size1() );
        //calculation loop (as long as there is no convergence)
        //WARNING: iteration never exceeds
        while( error > error_tolerance )
        {
            for( int i=0; i<n; i++ )
            {
                for( int j=i+1; j<n; j++ )
                {
                    double theta = 0.0;
                    if( MathUtilsType::Abs( A(i,j) ) >= zero_tolerance )
                    {
                        if( MathUtilsType::Abs( A(i,i)-A(j,j) ) > 0.0 )
                        {
                            theta = 0.5*atan(2*A(i,j)/(A(i,i)-A(j,j)));
                        }
                        else theta = 0.25*Pi();
                    }
                    MatrixType T = IdentityMatrix( n );

                    T(i,i) = cos(theta);
                    T(i,j) = -sin(theta);
                    T(j,i) = -T(i,j);
                    T(j,j) = T(i,i);

                    A = Mult( A, T );
                    MatrixType TT = Transpose(T);
                    A = Mult( TT, A );
                    V = Mult( V, T );
                }
            }
            double sTot = 0.0;
            double sDiag = 0.0;
            for( unsigned int i=0; i<A.size1(); i++ )
            {
                for( unsigned int j=0; j<A.size2(); j++ )
                {
                    sTot += MathUtilsType::Abs(A(i,j));
                }
                sDiag+= MathUtilsType::Abs(A(i,i));
            }
            error=(sTot-sDiag)/sDiag;
        }
        //sorting eigenvalues
        int maxIndex = 0;
        TDataType maxEv = A(0,0);
        for( unsigned int i=0; i<A.size1(); i++ )
        {
            for( unsigned int j=i; j<A.size1(); j++ )
            {
                //searching current maximum
                if( A(j,j) > maxEv )
                {
                    maxIndex = j;
                    maxEv = A(j,j);
                }
                //swapping eigenvalue matrix
                TDataType dummy = A(i,i);
                A(i,i) = A(maxIndex,maxIndex);
                A(maxIndex,maxIndex) = dummy;
                //swapping eigenvector matrix
                for( unsigned int k=0; k<A.size2(); k++ )
                {
                    dummy = V(k,i);
                    V(k,i) = V(k,maxIndex);
                    V(k,maxIndex) = dummy;
                }
            }

        }
    }//EigenVectors

    /**
     * creates identity matrix.
     * Given matrix will be overwritten
     * @param given matrix to be overwritten by identity matrix
     */
    static inline MatrixType IdentityMatrix( SizeType size )
    {
        MatrixType A = ZeroMatrix( size );
        for( unsigned int i=0; i<size ; i++ )
        {
            A(i,i) = 1.0;
        }
        return A;
    }//IdentityMatrix

    /**
     * Adds two matrices. first argument is overwritten by sum of both
     * Matrices are assumed to be of same dimension (no check on boundaries is made!)
     * @param A first matrix argument (overwritten by solution)
     * @param B second matrix argument
     */
    static inline void Add( MatrixType& A, MatrixType& B )
    {
        for( unsigned int i=0; i<A.size1(); i++ )
        {
            for( unsigned int j=0; j<A.size2(); j++ )
            {
                A(i,j) += B(i,j);
            }
        }
    }

    /**
     * multiplies two matrices. Performs operation \f$ C = A B \f$
     * @param A matrix A
     * @param B matrix B
     * @return matrix \f$ C = A B \f$
     */
    static inline MatrixType Mult( MatrixType& A,
                                   MatrixType& B)
    {
        MatrixType C(A.size1(),B.size2());
        for( unsigned int i=0; i<A.size1(); i++ )
        {
            for( unsigned int j=0; j<B.size2(); j++ )
            {
                for( unsigned int k=0; k<A.size2(); k++ )
                {
                    C(i,j) += B(k,j)*A(i,k);
                }
            }
        }
        return C;
    }//Mult

    /**
     * multiplies a matrix by a scalar
     */
    static inline void Mult( MatrixType& M, TDataType a )
    {
        for( unsigned int i=0; i<M.size1(); i++ )
        {
            for( unsigned int j=0; j<M.size2(); j++ )
            {
                M(i,j) = M(i,j)*a;
            }
        }
    }

    /**
     * multiplies a vector by a scalar
     */
    static inline void Mult( VectorType& v, TDataType a )
    {
        for( unsigned int i=0; i<v.size(); i++ )
        {
            v(i) = v(i)*a;
        }
    }//Mult

    /**
     * transposes matrix A. Matrix A is not overwritten!
     * @param A the given Matrix
     * @return the transposed matrix \f$ A^T \f$
     */
    static inline MatrixType Transpose( MatrixType& A )
    {
        MatrixType AT = ZeroMatrix(A.size2(),A.size1());
        for( unsigned int i=0; i<A.size1(); i++ )
        {
            for( unsigned int j=0; j<A.size2(); j++ )
            {
                AT(j,i) = A(i,j);
            }
        }
        return AT;
    }//Transpose

    /**
     * normalises a vector. Vector is scaled by \f$ V_{norm} = \frac{V}{|V|} \f$
     */
    static inline void Normalize( VectorType& v )
    {
        Mult( v, 1.0/(MathUtilsType::Norm( v )) );
    }

    /**
     * converts a strain vector into a matrix. Strains are assumed to be stored
     * in the following way:
     * \f$ [ e11, e22, e33, 2*e12, 2*e23, 2*e13 ] \f$ for 3D case and
     * \f$ [ e11, e22, 2*e12 ] \f$ fir 2D case.
     * Hence the deviatoric components of the strain vector are divided by 2
     * while they are stored into the matrix
     * @param Strains the given strain vector
     * @return the corresponding strain tensor in matrix form
     */
    static inline MatrixType StrainVectorToTensor( const VectorType& Strains )
    {
        KRATOS_TRY

        MatrixType StrainTensor;
        //KRATOS_WATCH(Strains)
        if (Strains.size() == 3)
        {
            StrainTensor.resize(2, 2, false);
            //KRATOS_WATCH(StrainTensor)
            StrainTensor(0,0) = Strains[0];
            StrainTensor(0,1) = 0.5*Strains[2];
            StrainTensor(1,0) = 0.5*Strains[2];
            StrainTensor(1,1) = Strains[1];
        }
        else if (Strains.size() == 4)
        {
            StrainTensor.resize(3, 3, false);
            //KRATOS_WATCH(StrainTensor)
            StrainTensor(0,0) = Strains[0];
            StrainTensor(0,1) = 0.5*Strains[2];
            StrainTensor(1,0) = 0.5*Strains[2];
            StrainTensor(1,1) = Strains[1];
            StrainTensor(2,2) = Strains[3];
        }
        else if (Strains.size() == 6)
        {
            StrainTensor.resize(3, 3, false);
            StrainTensor(0,0) = Strains[0];
            StrainTensor(0,1) = 0.5*Strains[3];
            StrainTensor(0,2) = 0.5*Strains[5];
            StrainTensor(1,0) = 0.5*Strains[3];
            StrainTensor(1,1) = Strains[1];
            StrainTensor(1,2) = 0.5*Strains[4];
            StrainTensor(2,0) = 0.5*Strains[5];
            StrainTensor(2,1) = 0.5*Strains[4];
            StrainTensor(2,2) = Strains[2];
        }

        //KRATOS_WATCH(StrainTensor)
        return StrainTensor;
        KRATOS_CATCH("")
    }

    template<typename TVectorType, typename TMatrixType>
    static inline void StrainVectorToTensor(const TVectorType& StrainVector, TMatrixType& StrainTensor)
    {
        if(StrainVector.size() == 3) // plane strain
        {
            noalias(StrainTensor) = ZeroMatrix(3, 3);
            StrainTensor(0, 0) = StrainVector(0);
            StrainTensor(0, 1) = 0.5 * StrainVector(2);
            StrainTensor(1, 0) = 0.5 * StrainVector(2);
            StrainTensor(1, 1) = StrainVector(1);
        }
        else if(StrainVector.size() == 4) // axisymmetric
        {
            noalias(StrainTensor) = ZeroMatrix(3, 3);
            StrainTensor(0, 0) = StrainVector(0);
            StrainTensor(0, 1) = 0.5 * StrainVector(2);
            StrainTensor(1, 0) = 0.5 * StrainVector(2);
            StrainTensor(1, 1) = StrainVector(1);
            StrainTensor(2, 2) = StrainVector(3);
        }
        else if(StrainVector.size() == 6) // 3d
        {
            StrainTensor(0, 0) = StrainVector(0);
            StrainTensor(0, 1) = 0.5 * StrainVector(3);
            StrainTensor(0, 2) = 0.5 * StrainVector(5);
            StrainTensor(1, 0) = 0.5 * StrainVector(3);
            StrainTensor(1, 1) = StrainVector(1);
            StrainTensor(1, 2) = 0.5 * StrainVector(4);
            StrainTensor(2, 0) = 0.5 * StrainVector(5);
            StrainTensor(2, 1) = 0.5 * StrainVector(4);
            StrainTensor(2, 2) = StrainVector(2);
        }
    }

    static inline VectorType TensorToStrainVector( const MatrixType& Tensor )
    {
        KRATOS_TRY
        VectorType StrainVector;

        if (Tensor.size1()==2)
        {
            StrainVector.resize(3);
            noalias(StrainVector) = ZeroVector(3);
            StrainVector(0) = Tensor(0,0);
            StrainVector(1) = Tensor(1,1);
            StrainVector(2) = 2.00*Tensor(0,1);
        }
        else if (Tensor.size1()==3)
        {
            StrainVector.resize(6);
            noalias(StrainVector) = ZeroVector(6);
            StrainVector(0) = Tensor(0,0);
            StrainVector(1) = Tensor(1,1);
            StrainVector(2) = Tensor(2,2);
            StrainVector(3) = 2.0*Tensor(0,1);
            StrainVector(4) = 2.0*Tensor(1,2);
            StrainVector(5) = 2.0*Tensor(0,2);
        }

        return StrainVector;
        KRATOS_CATCH("")
    }

    template<typename TVectorType, typename TMatrixType>
    static inline void StrainTensorToVector(const TMatrixType& StrainTensor, TVectorType& StrainVector)
    {
        if(StrainVector.size() == 3)
        {
            StrainVector[0] = StrainTensor(0, 0);
            StrainVector[1] = StrainTensor(1, 1);
            StrainVector[2] = 2.0*StrainTensor(0, 1);
        }
        else if(StrainVector.size() == 4)
        {
            StrainVector[0] = StrainTensor(0, 0);
            StrainVector[1] = StrainTensor(1, 1);
            StrainVector[2] = 2.0*StrainTensor(0, 1);
            StrainVector[3] = StrainTensor(2, 2);
        }
        else if(StrainVector.size() == 6)
        {
            StrainVector[0] = StrainTensor(0, 0);
            StrainVector[1] = StrainTensor(1, 1);
            StrainVector[2] = StrainTensor(2, 2);
            StrainVector[3] = 2.0*StrainTensor(0, 1);
            StrainVector[4] = 2.0*StrainTensor(1, 2);
            StrainVector[5] = 2.0*StrainTensor(0, 2);
        }
    }

    template<typename TVectorType, typename TMatrixType>
    static inline void StressVectorToTensor(const TVectorType& StressVector, TMatrixType& StressTensor)
    {
        if(StressVector.size() == 3) // plane strain
        {
            noalias(StressTensor) = ZeroMatrix(3, 3);
            StressTensor(0, 0) = StressVector(0);
            StressTensor(0, 1) = StressVector(2);
            StressTensor(1, 0) = StressVector(2);
            StressTensor(1, 1) = StressVector(1);
        }
        else if(StressVector.size() == 4) // axisymmetric
        {
            noalias(StressTensor) = ZeroMatrix(3, 3);
            StressTensor(0, 0) = StressVector(0);
            StressTensor(0, 1) = StressVector(2);
            StressTensor(1, 0) = StressVector(2);
            StressTensor(1, 1) = StressVector(1);
            StressTensor(2, 2) = StressVector(3);
        }
        else if(StressVector.size() == 6) // 3d
        {
            StressTensor(0, 0) = StressVector(0);
            StressTensor(0, 1) = StressVector(3);
            StressTensor(0, 2) = StressVector(5);
            StressTensor(1, 0) = StressVector(3);
            StressTensor(1, 1) = StressVector(1);
            StressTensor(1, 2) = StressVector(4);
            StressTensor(2, 0) = StressVector(5);
            StressTensor(2, 1) = StressVector(4);
            StressTensor(2, 2) = StressVector(2);
        }
    }

    template<typename TVectorType, typename TMatrixType>
    static inline void StressTensorToVector(const TMatrixType& StressTensor, TVectorType& StressVector)
    {
        if(StressVector.size() == 3) // plane strain
        {
            StressVector[0] = StressTensor(0, 0);
            StressVector[1] = StressTensor(1, 1);
            StressVector[2] = StressTensor(0, 1);
        }
        if(StressVector.size() == 4) // axisymmetric
        {
            StressVector[0] = StressTensor(0, 0);
            StressVector[1] = StressTensor(1, 1);
            StressVector[2] = StressTensor(0, 1);
            StressVector[3] = StressTensor(2, 2);
        }
        else if(StressVector.size() == 6) // 3d
        {
            StressVector[0] = StressTensor(0, 0);
            StressVector[1] = StressTensor(1, 1);
            StressVector[2] = StressTensor(2, 2);
            StressVector[3] = StressTensor(0, 1);
            StressVector[4] = StressTensor(1, 2);
            StressVector[5] = StressTensor(0, 2);
        }
    }

    template<typename TVectorType, typename TMatrixType>
    static inline void AddStressTensorToVector(double c, const TMatrixType& StressTensor, TVectorType& StressVector)
    {
        if(StressVector.size() == 3)
        {
            StressVector[0] += c*StressTensor(0, 0);
            StressVector[1] += c*StressTensor(1, 1);
            StressVector[2] += c*StressTensor(0, 1);
        }
        else if(StressVector.size() == 4)
        {
            StressVector[0] += c*StressTensor(0, 0);
            StressVector[1] += c*StressTensor(1, 1);
            StressVector[2] += c*StressTensor(0, 1);
            StressVector[3] += c*StressTensor(2, 2);
        }
        else if(StressVector.size() == 6)
        {
            StressVector[0] += c*StressTensor(0, 0);
            StressVector[1] += c*StressTensor(1, 1);
            StressVector[2] += c*StressTensor(2, 2);
            StressVector[3] += c*StressTensor(0, 1);
            StressVector[4] += c*StressTensor(1, 2);
            StressVector[5] += c*StressTensor(0, 2);
        }
    }

    /**
    * Builds the Inverse of Matrix input
    * @param input the given Matrix
    * @param inverse inverse of the given Matrix
    */
    static int InvertMatrix( const MatrixType& input, MatrixType& inverse )
    {
        int singular = 0;
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<std::size_t> pmatrix;
        MatrixType A(input);
        pmatrix pm(A.size1());
        singular = lu_factorize(A,pm);
        inverse.assign( IdentityMatrix(A.size1()));
        lu_substitute(A, pm, inverse);
        return singular;
    }

    /**
    * Compute the Frobenius norm of a given second order tensor
    * @param Tensor the given second order tensor
    * @return the norm of the given tensor
    */
    template<typename TMatrixType>
    static double normTensor(const TMatrixType& Tensor)
    {
        double result=0.0;
        for(unsigned int i=0; i<Tensor.size1(); i++)
            for(unsigned int j=0; j<Tensor.size2(); j++)
                result += Tensor(i,j)*Tensor(i,j);
        return sqrt(result);
    }

    /**
    * Compute the Frobenius norm of a given fourth order tensor
    * @param Tensor the given fourth order tensor
    * @return the norm of the given tensor
    */
    static double normTensor(const Fourth_Order_Tensor& T)
    {
        double v = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                v += pow(norm_frobenius(T[i][j]), 2);
        return sqrt(v);
    }

    /**
    * Transforms a given 6*1 Vector to a corresponing symmetric Tensor of second order (3*3)
    * @param Vector the given vector
    * @param Tensor the symmetric second order tensor
    */
    template<typename TVectorType, typename TMatrixType>
    static inline void VectorToTensor(const TVectorType& Stress, TMatrixType& Tensor)
    {
        if(Stress.size()==6)
        {
            Tensor.resize(3,3);
            Tensor(0,0)= Stress(0);
            Tensor(0,1)= Stress(3);
            Tensor(0,2)= Stress(5);
            Tensor(1,0)= Stress(3);
            Tensor(1,1)= Stress(1);
            Tensor(1,2)= Stress(4);
            Tensor(2,0)= Stress(5);
            Tensor(2,1)= Stress(4);
            Tensor(2,2)= Stress(2);
        }
        if(Stress.size()==3)
        {
            Tensor.resize(2,2);
            Tensor(0,0)= Stress(0);
            Tensor(0,1)= Stress(2);
            Tensor(1,0)= Stress(2);
            Tensor(1,1)= Stress(1);
        }
        return;
    }

    /**
    * Transforms a given symmetric Tensor of second order (3*3) to a corresponing 6*1 Vector
    * @param Tensor the given symmetric second order tensor
    * @param Vector the vector
    */
    template<typename TVectorType, typename TMatrixType>
    static void TensorToVector( const TMatrixType& Tensor, TVectorType& Vector)
    {
        //if(Vector.size()!= 6)
        unsigned int  dim  =  Tensor.size1();
        if (dim==3)
        {
            Vector.resize(6,false);
            Vector(0)= Tensor(0,0);
            Vector(1)= Tensor(1,1);
            Vector(2)= Tensor(2,2);
            Vector(3)= Tensor(0,1);
            Vector(4)= Tensor(1,2);
            Vector(5)= Tensor(2,0);
        }
        else if(dim==2)
        {
            Vector.resize(3,false);
            Vector(0)= Tensor(0,0);
            Vector(1)= Tensor(1,1);
            Vector(2)= Tensor(0,1);
        }
        return;
    }

    /// Transformation from a fourth order tensor to symmetric matrix. The matrix
    /// can then be used to multiply with vector form of strain tensor
    /// This transformation uses the notation [o_xx o_yy o_zz o_xy o_yz o_xz]
    /// Note that, this already accounts for factor 2 in [e_xx e_yy e_zz 2e_xy 2e_yz 2e_xz] (see https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Constitutive_relations)
    template<typename TMatrixType>
    static inline void TensorToMatrix(const Fourth_Order_Tensor& Tensor, TMatrixType& Matrix)
    {
        if (Matrix.size1() == 6)
        {
            Matrix(0, 0) = Tensor[0][0](0, 0); // xx-xx
            Matrix(0, 1) = Tensor[0][0](1, 1); // xx-yy
            Matrix(0, 2) = Tensor[0][0](2, 2); // xx-zz
            Matrix(0, 3) = Tensor[0][0](0, 1); // xx-xy
            Matrix(0, 4) = Tensor[0][0](1, 2); // xx-yz
            Matrix(0, 5) = Tensor[0][0](0, 2); // xx-xz

            Matrix(1, 0) = Tensor[1][1](0, 0);
            Matrix(1, 1) = Tensor[1][1](1, 1);
            Matrix(1, 2) = Tensor[1][1](2, 2);
            Matrix(1, 3) = Tensor[1][1](0, 1);
            Matrix(1, 4) = Tensor[1][1](1, 2);
            Matrix(1, 5) = Tensor[1][1](0, 2);

            Matrix(2, 0) = Tensor[2][2](0, 0);
            Matrix(2, 1) = Tensor[2][2](1, 1);
            Matrix(2, 2) = Tensor[2][2](2, 2);
            Matrix(2, 3) = Tensor[2][2](0, 1);
            Matrix(2, 4) = Tensor[2][2](1, 2);
            Matrix(2, 5) = Tensor[2][2](0, 2);

            Matrix(3, 0) = Tensor[0][1](0, 0);
            Matrix(3, 1) = Tensor[0][1](1, 1);
            Matrix(3, 2) = Tensor[0][1](2, 2);
            Matrix(3, 3) = Tensor[0][1](0, 1);
            Matrix(3, 4) = Tensor[0][1](1, 2);
            Matrix(3, 5) = Tensor[0][1](0, 2);

            Matrix(4, 0) = Tensor[1][2](0, 0);
            Matrix(4, 1) = Tensor[1][2](1, 1);
            Matrix(4, 2) = Tensor[1][2](2, 2);
            Matrix(4, 3) = Tensor[1][2](0, 1);
            Matrix(4, 4) = Tensor[1][2](1, 2);
            Matrix(4, 5) = Tensor[1][2](0, 2);

            Matrix(5, 0) = Tensor[0][2](0, 0);
            Matrix(5, 1) = Tensor[0][2](1, 1);
            Matrix(5, 2) = Tensor[0][2](2, 2);
            Matrix(5, 3) = Tensor[0][2](0, 1);
            Matrix(5, 4) = Tensor[0][2](1, 2);
            Matrix(5, 5) = Tensor[0][2](0, 2);
        }
        else if(Matrix.size1() == 4)
        {
            Matrix(0, 0) = Tensor[0][0](0, 0); // xx-xx
            Matrix(0, 1) = Tensor[0][0](1, 1); // xx-yy
            Matrix(0, 2) = Tensor[0][0](0, 1); // xx-xy
            Matrix(0, 3) = Tensor[0][0](2, 2); // xx-zz

            Matrix(1, 0) = Tensor[1][1](0, 0); // yy-xx
            Matrix(1, 1) = Tensor[1][1](1, 1); // yy-yy
            Matrix(1, 2) = Tensor[1][1](0, 1); // yy-xy
            Matrix(1, 3) = Tensor[1][1](2, 2); // yy-zz

            Matrix(2, 0) = Tensor[0][1](0, 0); // xy-xx
            Matrix(2, 1) = Tensor[0][1](1, 1); // xy-yy
            Matrix(2, 2) = Tensor[0][1](0, 1); // xy-xy
            Matrix(2, 3) = Tensor[0][1](2, 2); // xy-zz

            Matrix(3, 0) = Tensor[2][2](0, 0); // zz-xx
            Matrix(3, 1) = Tensor[2][2](1, 1); // zz-yy
            Matrix(3, 2) = Tensor[2][2](0, 1); // zz-xy
            Matrix(3, 3) = Tensor[2][2](2, 2); // zz-zz
        }
        else if(Matrix.size1() == 3)
        {
            Matrix(0, 0) = Tensor[0][0](0, 0); // xx-xx
            Matrix(0, 1) = Tensor[0][0](1, 1); // xx-yy
            Matrix(0, 2) = Tensor[0][0](0, 1); // xx-xy
            Matrix(1, 0) = Tensor[1][1](0, 0); // yy-xx
            Matrix(1, 1) = Tensor[1][1](1, 1); // yy-yy
            Matrix(1, 2) = Tensor[1][1](0, 1); // yy-xy
            Matrix(2, 0) = Tensor[0][1](0, 0); // xy-xx
            Matrix(2, 1) = Tensor[0][1](1, 1); // xy-yy
            Matrix(2, 2) = Tensor[0][1](0, 1); // xy-xy
        }
        else
            KRATOS_ERROR << "Invalid matrix size (" << Matrix.size1() << ", " << Matrix.size2() << ")";
    }

    /// Transformation from a fourth order tensor to unsymmetric matrix. The matrix
    /// can then be used to multiply with vector form of an unsymmetric tensor
    /// This transformation uses the notation [o_xx o_yx o_zx o_xy o_yy o_zy o_xz o_yz o_zz] in 3D and
    /// [o_xx o_yx o_xy o_yy] in 2D
    /// For axisymmetric, it is [o_xx o_yx o_xy o_yy o_zz]
    /// Reference: Souza de Neto, Computational Plasticity, Appendix D.2.1
    template<typename TMatrixType>
    static inline void TensorToUnsymmetricMatrix(const Fourth_Order_Tensor& Tensor, TMatrixType& Matrix)
    {
        if (Matrix.size1() == 9)
        {
            /// ((DO NOT DELETE)) (KEEP AS REFERENCE)
            // Matrix(0, 0) = Tensor[0][0](0, 0); // xx-xx
            // Matrix(0, 1) = Tensor[0][0](1, 0); // xx-yx
            // Matrix(0, 2) = Tensor[0][0](2, 0); // xx-zx
            // Matrix(0, 3) = Tensor[0][0](0, 1); // xx-xy
            // Matrix(0, 4) = Tensor[0][0](1, 1); // xx-yy
            // Matrix(0, 5) = Tensor[0][0](2, 1); // xx-zy
            // Matrix(0, 6) = Tensor[0][0](0, 2); // xx-xz
            // Matrix(0, 7) = Tensor[0][0](1, 2); // xx-yz
            // Matrix(0, 8) = Tensor[0][0](2, 2); // xx-zz
            //
            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    for (unsigned int k = 0; k < 3; ++k)
                    {
                        for (unsigned int l = 0; l < 3; ++l)
                        {
                            Matrix(3*j+i, 3*l+k) = Tensor[i][j](k, l);
                        }
                    }
                }
            }
        }
        else if(Matrix.size1() == 5)
        {
            Matrix(0, 0) = Tensor[0][0](0, 0); // xx-xx
            Matrix(0, 1) = Tensor[0][0](1, 0); // xx-yx
            Matrix(0, 2) = Tensor[0][0](0, 1); // xx-xy
            Matrix(0, 3) = Tensor[0][0](1, 1); // xx-yy
            Matrix(0, 4) = Tensor[0][0](2, 2); // xx-zz
            //
            Matrix(1, 0) = Tensor[1][0](0, 0); // yx-xx
            Matrix(1, 1) = Tensor[1][0](1, 0); // yx-yx
            Matrix(1, 2) = Tensor[1][0](0, 1); // yx-xy
            Matrix(1, 3) = Tensor[1][0](1, 1); // yx-yy
            Matrix(1, 4) = Tensor[1][0](2, 2); // yx-zz
            //
            Matrix(2, 0) = Tensor[0][1](0, 0); // xy-xx
            Matrix(2, 1) = Tensor[0][1](1, 0); // xy-yx
            Matrix(2, 2) = Tensor[0][1](0, 1); // xy-xy
            Matrix(2, 3) = Tensor[0][1](1, 1); // xy-yy
            Matrix(2, 4) = Tensor[0][1](2, 2); // xy-zz
            //
            Matrix(3, 0) = Tensor[1][1](0, 0); // yy-xx
            Matrix(3, 1) = Tensor[1][1](1, 0); // yy-yx
            Matrix(3, 2) = Tensor[1][1](0, 1); // yy-xy
            Matrix(3, 3) = Tensor[1][1](1, 1); // yy-yy
            Matrix(3, 4) = Tensor[1][1](2, 2); // yy-zz
            //
            Matrix(4, 0) = Tensor[2][2](0, 0); // zz-xx
            Matrix(4, 1) = Tensor[2][2](1, 0); // zz-yx
            Matrix(4, 2) = Tensor[2][2](0, 1); // zz-xy
            Matrix(4, 3) = Tensor[2][2](1, 1); // zz-yy
            Matrix(4, 4) = Tensor[2][2](2, 2); // zz-zz
        }
        else if(Matrix.size1() == 4)
        {
            /// ((DO NOT DELETE)) (KEEP AS REFERENCE)
            // Matrix(0, 0) = Tensor[0][0](0, 0); // xx-xx
            // Matrix(0, 1) = Tensor[0][0](1, 0); // xx-yx
            // Matrix(0, 2) = Tensor[0][0](0, 1); // xx-xy
            // Matrix(0, 3) = Tensor[0][0](1, 1); // xx-yy
            // //
            // Matrix(1, 0) = Tensor[1][0](0, 0); // yx-xx
            // Matrix(1, 1) = Tensor[1][0](1, 0); // yx-yx
            // Matrix(1, 2) = Tensor[1][0](0, 1); // yx-xy
            // Matrix(1, 3) = Tensor[1][0](1, 1); // yx-yy
            // //
            // Matrix(2, 0) = Tensor[0][1](0, 0); // xy-xx
            // Matrix(2, 1) = Tensor[0][1](1, 0); // xy-yx
            // Matrix(2, 2) = Tensor[0][1](0, 1); // xy-xy
            // Matrix(2, 3) = Tensor[0][1](1, 1); // xy-yy
            // //
            // Matrix(3, 0) = Tensor[1][1](0, 0); // yy-xx
            // Matrix(3, 1) = Tensor[1][1](1, 0); // yy-yx
            // Matrix(3, 2) = Tensor[1][1](0, 1); // yy-xy
            // Matrix(3, 3) = Tensor[1][1](1, 1); // yy-yy
            /// Or
            for (unsigned int i = 0; i < 2; ++i)
            {
                for (unsigned int j = 0; j < 2; ++j)
                {
                    for (unsigned int k = 0; k < 2; ++k)
                    {
                        for (unsigned int l = 0; l < 2; ++l)
                        {
                            Matrix(2*j+i, 2*l+k) = Tensor[i][j](k, l);
                        }
                    }
                }
            }
        }
        else
            KRATOS_ERROR << "Invalid matrix size (" << Matrix.size1() << ", " << Matrix.size2() << ")";
    }

    template<typename TMatrixType>
    static inline void UnsymmetricMatrixToTensor(const TMatrixType& Matrix, Fourth_Order_Tensor& Tensor)
    {
        if (Matrix.size1() == 9)
        {
            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    for (unsigned int k = 0; k < 3; ++k)
                    {
                        for (unsigned int l = 0; l < 3; ++l)
                        {
                            Tensor[i][j](k, l) = Matrix(3*j+i, 3*l+k);
                        }
                    }
                }
            }
        }
        else
            KRATOS_ERROR << "If matrix size is not 9, the 4th order tensor can't be filled since information is not sufficient";
    }

    /// Inversed operation of TensorToMatrix
    template<typename TMatrixType>
    static inline void MatrixToTensor(const TMatrixType& Matrix, Fourth_Order_Tensor& Tensor)
    {
        CalculateFourthOrderZeroTensor(Tensor);

        if (Matrix.size1() == 6)
        {
            Tensor[0][0](0, 0) = Matrix(0, 0); // xx-xx
            Tensor[0][0](1, 1) = Matrix(0, 1); // xx-yy
            Tensor[0][0](2, 2) = Matrix(0, 2); // xx-zz
            Tensor[0][0](0, 1) = Matrix(0, 3); // xx-xy
            Tensor[0][0](1, 2) = Matrix(0, 4); // xx-yz
            Tensor[0][0](0, 2) = Matrix(0, 5); // xx-xz

            Tensor[1][1](0, 0) = Matrix(1, 0);
            Tensor[1][1](1, 1) = Matrix(1, 1);
            Tensor[1][1](2, 2) = Matrix(1, 2);
            Tensor[1][1](0, 1) = Matrix(1, 3);
            Tensor[1][1](1, 2) = Matrix(1, 4);
            Tensor[1][1](0, 2) = Matrix(1, 5);

            Tensor[2][2](0, 0) = Matrix(2, 0);
            Tensor[2][2](1, 1) = Matrix(2, 1);
            Tensor[2][2](2, 2) = Matrix(2, 2);
            Tensor[2][2](0, 1) = Matrix(2, 3);
            Tensor[2][2](1, 2) = Matrix(2, 4);
            Tensor[2][2](0, 2) = Matrix(2, 5);

            Tensor[0][1](0, 0) = Matrix(3, 0);
            Tensor[0][1](1, 1) = Matrix(3, 1);
            Tensor[0][1](2, 2) = Matrix(3, 2);
            Tensor[0][1](0, 1) = Matrix(3, 3);
            Tensor[0][1](1, 2) = Matrix(3, 4);
            Tensor[0][1](0, 2) = Matrix(3, 5);

            Tensor[1][2](0, 0) = Matrix(4, 0);
            Tensor[1][2](1, 1) = Matrix(4, 1);
            Tensor[1][2](2, 2) = Matrix(4, 2);
            Tensor[1][2](0, 1) = Matrix(4, 3);
            Tensor[1][2](1, 2) = Matrix(4, 4);
            Tensor[1][2](0, 2) = Matrix(4, 5);

            Tensor[0][2](0, 0) = Matrix(5, 0);
            Tensor[0][2](1, 1) = Matrix(5, 1);
            Tensor[0][2](2, 2) = Matrix(5, 2);
            Tensor[0][2](0, 1) = Matrix(5, 3);
            Tensor[0][2](1, 2) = Matrix(5, 4);
            Tensor[0][2](0, 2) = Matrix(5, 5);

            for (unsigned int j = 0; j < 3; ++j)
            {
                for (unsigned int i = 0; i <= j; ++i)
                {
                    for (unsigned int k = 0; k < 3; ++k)
                    {
                        for (unsigned int l = 0; l < k; ++l)
                        {
                            Tensor[i][j](k, l) = Tensor[i][j](l, k);
                        }
                    }
                }
            }

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < i; ++j)
                {
                    for (unsigned int k = 0; k < 3; ++k)
                    {
                        for (unsigned int l = 0; l < 3; ++l)
                        {
                            Tensor[i][j](k, l) = Tensor[j][i](k, l);
                        }
                    }
                }
            }
        }
        else
            KRATOS_ERROR << "If matrix size is not 6, the 4th order tensor can't be filled since information is not sufficient";
    }

    // THis uses the notation [o_xx o_yy o_zz o_xy o_xz o_yz]
    // Note that, this already accounts for factor 2 in [e_xx e_yy e_zz 2e_xy 2e_xz 2e_yz] (see https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Constitutive_relations)
    template<typename TMatrixType>
    static inline void TensorToMatrix2(const Fourth_Order_Tensor& Tensor, TMatrixType& Matrix)
    {
        // Simetrias seguras
        //  Cijkl = Cjilk;
        //  Cijkl = Cklji;
        if (Tensor[0].size()== 3)
        {
            // Tensor de cuarto orden cuyos componentes correspondes a una matriz de 3x3
            if(Matrix.size1()!=6 || Matrix.size2()!=6)
                Matrix.resize(6,6,false);
            Matrix(0,0) = Tensor[0][0](0,0);
            Matrix(0,1) = Tensor[0][0](1,1);
            Matrix(0,2) = Tensor[0][0](2,2);
            Matrix(0,3) = Tensor[0][0](0,1);
            Matrix(0,4) = Tensor[0][0](0,2);
            Matrix(0,5) = Tensor[0][0](1,2);

            Matrix(1,0) = Tensor[1][1](0,0);
            Matrix(1,1) = Tensor[1][1](1,1);
            Matrix(1,2) = Tensor[1][1](2,2);
            Matrix(1,3) = Tensor[1][1](0,1);
            Matrix(1,4) = Tensor[1][1](0,2);
            Matrix(1,5) = Tensor[1][1](1,2);

            Matrix(2,0) = Tensor[2][2](0,0);
            Matrix(2,1) = Tensor[2][2](1,1);
            Matrix(2,2) = Tensor[2][2](2,2);
            Matrix(2,3) = Tensor[2][2](0,1);
            Matrix(2,4) = Tensor[2][2](0,2);
            Matrix(2,5) = Tensor[2][2](1,2);

            Matrix(3,0) = Tensor[0][1](0,0);
            Matrix(3,1) = Tensor[0][1](1,1);
            Matrix(3,2) = Tensor[0][1](2,2);
            Matrix(3,3) = Tensor[0][1](0,1);
            Matrix(3,4) = Tensor[0][1](0,2);
            Matrix(3,5) = Tensor[0][1](1,2);

            Matrix(4,0) = Tensor[0][2](0,0);
            Matrix(4,1) = Tensor[0][2](1,1);
            Matrix(4,2) = Tensor[0][2](2,2);
            Matrix(4,3) = Tensor[0][2](0,1);
            Matrix(4,4) = Tensor[0][2](0,2);
            Matrix(4,5) = Tensor[0][2](1,2);

            Matrix(5,0) = Tensor[1][2](0,0);
            Matrix(5,1) = Tensor[1][2](1,1);
            Matrix(5,2) = Tensor[1][2](2,2);
            Matrix(5,3) = Tensor[1][2](0,1);
            Matrix(5,4) = Tensor[1][2](0,2);
            Matrix(5,5) = Tensor[1][2](1,2);
        }
        else
        {
            // Tensor de cuarto orden cuyos componentes correspondes a una matriz de 2x2
            if(Matrix.size1()!=3 || Matrix.size2()!=3)
                Matrix.resize(3,3,false);
            Matrix(0,0) = Tensor[0][0](0,0);
            Matrix(0,1) = Tensor[0][0](1,1);
            Matrix(0,2) = Tensor[0][0](0,1);
            Matrix(1,0) = Tensor[1][1](0,0);
            Matrix(1,1) = Tensor[1][1](1,1);
            Matrix(1,2) = Tensor[1][1](0,1);
            Matrix(2,0) = Tensor[0][1](0,0);
            Matrix(2,1) = Tensor[0][1](1,1);
            Matrix(2,2) = Tensor[0][1](0,1);
        }
        return;
    }

    /**
    * Transforms a given 6*6 Matrix to a corresponding 4th order tensor
    * @param Tensor the given Matrix
    * @param Vector the Tensor
    */
    template<typename TMatrixType, typename Fourth_Order_Tensor_Type>
    static void MatrixToTensor2(const TMatrixType& A, Fourth_Order_Tensor_Type& Tensor)
    {
        int help1 = 0;
        int help2 = 0;
        double coeff = 1.0;

        Tensor.resize(3);

        for(unsigned int i=0; i<3; i++)
        {
            Tensor[i].resize(3);
            for(unsigned int j=0; j<3; j++)
            {
                Tensor[i][j].resize(3,3,false);
                noalias(Tensor[i][j])= ZeroMatrix(3,3);
                for(unsigned int k=0; k<3; k++)
                    for(unsigned int l=0; l<3; l++)
                    {
                        if(i==j) help1= i;
                        else
                        {
                            if((i==0 && j==1) || (i==1 && j==0)) help1= 3;
                            if((i==1 && j==2) || (i==2 && j==1)) help1= 4;
                            if((i==2 && j==0) || (i==0 && j==2)) help1= 5;
                        }
                        if(k==l)
                        {
                            help2= k;
                            coeff=1.0;
                        }
                        else
                        {
                            coeff=0.5;
                            if((k==0 && l==1) || (k==1 && l==0)) help2= 3;
                            if((k==1 && l==2) || (k==2 && l==1)) help2= 4;
                            if((k==2 && l==0) || (k==0 && l==2)) help2= 5;
                        }

                        Tensor[i][j](k,l)= A(help1,help2)*coeff;
                    }
            }
        }

        return;
    }

    /**
    * Transforms a given 6*6 Matrix to a corresponing 4th order tensor
    * @param Tensor the given Matrix
    * @param Vector the Tensor
    */
    template<typename TMatrixType>
    static void MatrixToTensor(const TMatrixType& A, array_1d<double, 81>& Tensor)
    {
        int help1 = 0;
        int help2 = 0;
        double coeff = 1.0;
        std::fill(Tensor.begin(), Tensor.end(), 0.0);
        for(unsigned int i=0; i<3; i++)
        {
            for(unsigned int j=0; j<3; j++)
            {
                for(unsigned int k=0; k<3; k++)
                    for(unsigned int l=0; l<3; l++)
                    {
                        if(i==j) help1= i;
                        else
                        {
                            if((i==0 && j==1) || (i==1 && j==0)) help1= 3;
                            if((i==1 && j==2) || (i==2 && j==1)) help1= 4;
                            if((i==2 && j==0) || (i==0 && j==2)) help1= 5;
                        }
                        if(k==l)
                        {
                            help2= k;
                            coeff=1.0;
                        }
                        else
                        {
                            coeff=0.5;
                            if((k==0 && l==1) || (k==1 && l==0)) help2= 3;
                            if((k==1 && l==2) || (k==2 && l==1)) help2= 4;
                            if((k==2 && l==0) || (k==0 && l==2)) help2= 5;
                        }

                        Tensor[i*27+j*9+k*3+l]= A(help1,help2)*coeff;
                    }
            }
        }

        return;
    }
    /**
    * Transforms a given 4th order tensor to a corresponing 6*6 Matrix
    * @param Tensor the given Tensor
    * @param Vector the Matrix
    */
    template<typename TMatrixType>
    static void TensorToMatrix(const std::vector<std::vector<MatrixType> >& Tensor, TMatrixType& Matrix)
    {
        int help1 = 0;
        int help2 = 0;
        int help3 = 0;
        int help4 = 0;
        double coeff = 1.0;

        if(Matrix.size1()!=6 || Matrix.size2()!=6)
            Matrix.resize(6,6,false);

        for(unsigned int i=0; i<6; i++)
            for(unsigned int j=0; j<6; j++)
            {
                if(i<3)
                {
                    help1= i;
                    help2= i;
                }
                else
                {
                    if(i==3)
                    {
                        help1= 0;
                        help2= 1;
                    }
                    if(i==4)
                    {
                        help1= 1;
                        help2= 2;
                    }
                    if(i==5)
                    {
                        help1= 2;
                        help2= 0;
                    }
                }

                if(j<3)
                {
                    help3= j;
                    help4= j;
                    coeff= 1.0;
                }
                else
                {
                    if(j==3)
                    {
                        help3= 0;
                        help4= 1;
                    }
                    if(j==4)
                    {
                        help3= 1;
                        help4= 2;
                    }
                    if(j==5)
                    {
                        help3= 2;
                        help4= 0;
                    }
                    coeff= 2.0;
                }

                Matrix(i,j)= Tensor[help1][help2](help3,help4)*coeff;
            }

        return;
    }

    /**
     * Transforms a given 4th order tensor to a corresponing 6*6 Matrix
     * @param Tensor the given Tensor
     * @param Vector the Matrix
     */
    template<typename TMatrixType>
    static void TensorToMatrix( const array_1d<double, 81>& Tensor, TMatrixType& Matrix )
    {
        if(Matrix.size1()!=6 || Matrix.size2()!=6)
            Matrix.resize(6,6,false);

        Matrix(0,0) = Tensor[0];
        Matrix(0,1) = Tensor[4];
        Matrix(0,2) = Tensor[8];
        Matrix(0,3) = 2.0*Tensor[1];
        Matrix(0,4) = 2.0*Tensor[5];
        Matrix(0,5) = 2.0*Tensor[6];

        Matrix(1,0) = Tensor[36];
        Matrix(1,1) = Tensor[40];
        Matrix(1,2) = Tensor[44];
        Matrix(1,3) = 2.0*Tensor[37];
        Matrix(1,4) = 0.0*Tensor[41];
        Matrix(1,5) = 0.0*Tensor[42];

        Matrix(2,0) = Tensor[72];
        Matrix(2,1) = Tensor[76];
        Matrix(2,2) = Tensor[80];
        Matrix(2,3) = 2.0*Tensor[73];
        Matrix(2,4) = 2.0*Tensor[77];
        Matrix(2,5) = 2.0*Tensor[78];

        Matrix(3,0) = Tensor[9];
        Matrix(3,1) = Tensor[13];
        Matrix(3,2) = Tensor[18];
        Matrix(3,3) = 2.0*Tensor[10];
        Matrix(3,4) = 2.0*Tensor[14];
        Matrix(3,5) = 2.0*Tensor[15];

        Matrix(4,0) = Tensor[45];
        Matrix(4,1) = Tensor[49];
        Matrix(4,2) = Tensor[53];
        Matrix(4,3) = 2.0*Tensor[46];
        Matrix(4,4) = 0.0*Tensor[50];
        Matrix(4,5) = 0.0*Tensor[51];

        Matrix(5,0) = Tensor[54];
        Matrix(5,1) = Tensor[58];
        Matrix(5,2) = Tensor[62];
        Matrix(5,3) = 2.0*Tensor[55];
        Matrix(5,4) = 2.0*Tensor[59];
        Matrix(5,5) = 2.0*Tensor[60];

        return;
    }

    template<typename TMatrixType1, typename TMatrixType2>
    static void ExtractVolumetricDeviatoricTensor( const TMatrixType1& C, TMatrixType2& dev, double& vol )
    {
        vol = C(0, 0) + C(1, 1) + C(2, 2);
        noalias(dev) = C;
        dev(0, 0) -= vol / 3;
        dev(1, 1) -= vol / 3;
        dev(2, 2) -= vol / 3;
    }

    static inline void CalculateFourthOrderDeviatoricTensor( Fourth_Order_Tensor& C )
    {
        MatrixType kronecker(3, 3);
        noalias(kronecker) = ZeroMatrix(3, 3);
        for(unsigned int i = 0; i < 3; ++i)
            kronecker(i, i) = 1;

        C.resize(3);
        for(unsigned int i = 0; i < 3; ++i)
        {
            C[i].resize(3);
            for(unsigned int j = 0; j < 3; ++j)
            {
                C[i][j].resize(3, 3);
                noalias(C[i][j]) = ZeroMatrix(3, 3);
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                        C[i][j](k,l) = 0.5 * kronecker(i, k) * kronecker(j, l)
                                     + 0.5 * kronecker(i, l) * kronecker(j, k)
                                     - 1.0 / 3 * kronecker(i, j) * kronecker(k, l);
                }
            }
        }
    }

    static inline void CalculateFourthOrderZeroTensor( Fourth_Order_Tensor& C )
    {
        C.resize(3, false);
        for(unsigned int i = 0; i < 3; ++i)
        {
            C[i].resize(3, false);
            for(unsigned int j = 0; j < 3; ++j)
            {
                C[i][j].resize(3, 3, false);
                noalias(C[i][j]) = ZeroMatrix(3, 3);
            }
        }
    }

    static inline void CalculateFourthOrderSymmetricTensor( Fourth_Order_Tensor& C )
    {
        MatrixType kronecker = IdentityMatrix(3);
        C.resize(3);
        for(unsigned int i = 0; i < 3; ++i)
        {
            C[i].resize(3);
            for(unsigned int j = 0; j < 3; ++j)
            {
                C[i][j].resize(3, 3);
                noalias(C[i][j]) = ZeroMatrix(3, 3);
                for(unsigned int k = 0; k < 3; ++k)
                    for(unsigned int l = 0; l < 3; ++l)
                        C[i][j](k,l) = 0.5 * (kronecker(i, k) * kronecker(j, l) + kronecker(i, l) * kronecker(j, k));
            }
        }
    }

    /// REF: Eq. (2.100) Souze de Neto, Computational Plasticity
    /// Note that this tensor is not symmetric
    static inline void CalculateFourthOrderUnitTensor( Fourth_Order_Tensor& C )
    {
        MatrixType kronecker = IdentityMatrix(3);

        C.resize(3);
        for(unsigned int i = 0; i < 3; ++i)
        {
            C[i].resize(3);
            for(unsigned int j = 0; j < 3; ++j)
            {
                C[i][j].resize(3, 3);
                noalias(C[i][j]) = ZeroMatrix(3, 3);
                for(unsigned int k = 0; k < 3; ++k)
                    for(unsigned int l = 0; l < 3; ++l)
                        C[i][j](k,l) = kronecker(i, k) * kronecker(j, l);
            }
        }
    }

    /**
     * Zero out a given 4th order tensor
     * @param C the given Tensor
     */
    static void ZeroFourthOrderTensor( Fourth_Order_Tensor& C )
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                        C[i][j](k,l) = 0.0;
                }
            }
        }
    }

    /**
     * Scales a given 4th order tensor by a scalar (C = C*a)
     * @param C the given Tensor
     * @param alpha
     */
    static void ScaleFourthOrderTensor( Fourth_Order_Tensor& C, double alpha )
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                        C[i][j](k,l) *= alpha;
                }
            }
        }
    }

    /**
     * Copy the fourth order tensor A -> B
     * @param A the source Tensor
     * @param B the target tensor
     */
    static void CopyFourthOrderTensor( const Fourth_Order_Tensor& A, Fourth_Order_Tensor& B )
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                        B[i][j](k,l) = A[i][j](k,l);
                }
            }
        }
    }

    /**
     * Computes contraction of a fourth order tensor and matrix and add to a second order tensor (Result += alpha * (AA : B))
     * @param alpha
     * @param AA the fourth order Tensor
     * @param B the second order Tensor
     */
    static void ContractFourthOrderTensor(double alpha, const Fourth_Order_Tensor& A, const MatrixType& B, MatrixType& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                    {
                        Result(i, j) += alpha * A[i][j](k, l) * B(k, l);
                    }
                }
            }
        }
    }

    /**
     * Computes contraction of a fourth order tensor and matrix and add to a second order tensor (Result += alpha * (AA : B))
     * @param alpha
     * @param AA the fourth order Tensor
     * @param B the second order Tensor
     */
    static void ContractFourthOrderTensor(double alpha, const Fourth_Order_Tensor& A, const SymmetricMatrixType& B, SymmetricMatrixType& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = i; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                    {
                        Result(i, j) += alpha * A[i][j](k, l) * B(k, l);
                    }
                }
            }
        }
    }

    /**
     * Computes contraction of a fourth order tensor and matrix and add to a second order tensor (Result += alpha * (A : B))
     * @param alpha
     * @param A the second order Tensor
     * @param BB the fourth order Tensor
     */
    static void ContractFourthOrderTensor(double alpha, const MatrixType& A, const Fourth_Order_Tensor& BB, MatrixType& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                    {
                        Result(k, l) += alpha * A(i, j) * BB[i][j](k, l);
                    }
                }
            }
        }
    }

    /**
     * Computes contraction of a fourth order tensor and symmetric matrix and add to a second order tensor (Result += alpha * (A : B))
     * @param alpha
     * @param A the second order Tensor
     * @param BB the fourth order Tensor
     */
    static void ContractFourthOrderTensor(double alpha, const SymmetricMatrixType& A, const Fourth_Order_Tensor& BB, SymmetricMatrixType& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = k; l < 3; ++l)
                    {
                        Result(k, l) += alpha * A(i, j) * BB[i][j](k, l);
                    }
                }
            }
        }
    }

    /**
     * Computes outer product of two 2nd order tensors (Matrix) and add to a given 4th order tensor (Result += alpha * (A \odot B))
     * In the indices notation: Result(i,j,k,l) += alpha * A(i, j) * B(k, l)
     * @param C the given Tensor
     * @param alpha
     */
    template<typename TMatrixType1, typename TMatrixType2>
    static void OuterProductFourthOrderTensor(double alpha, const TMatrixType1& A, const TMatrixType2& B, Fourth_Order_Tensor& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                    {
                        Result[i][j](k, l) += alpha * A(i, j) * B(k, l);
                    }
                }
            }
        }
    }

    /**
     * In the indices notation: Result(i,j,k,l) += alpha * A(i, k) * B(j, l)
     * @param C the given Tensor
     * @param alpha
     * TODO make a better name
     */
    template<typename TMatrixType1, typename TMatrixType2>
    static void SpecialProduct1FourthOrderTensor(double alpha, const TMatrixType1& A, const TMatrixType2& B, Fourth_Order_Tensor& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                    {
                        Result[i][j](k, l) += alpha * A(i, k) * B(j, l);
                    }
                }
            }
        }
    }

    /**
     * In the indices notation: Result(i,j,k,l) += alpha * A(i, l) * B(j, k)
     * @param C the given Tensor
     * @param alpha
     * TODO make a better name
     */
    template<typename TMatrixType1, typename TMatrixType2>
    static void SpecialProduct2FourthOrderTensor(double alpha, const TMatrixType1& A, const TMatrixType2& B, Fourth_Order_Tensor& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                    {
                        Result[i][j](k, l) += alpha * A(i, l) * B(j, k);
                    }
                }
            }
        }
    }

    // C += alpha A
    static inline void AddFourthOrderTensor(double alpha, const Fourth_Order_Tensor& A, Fourth_Order_Tensor& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                noalias(Result[i][j]) += alpha * A[i][j];
            }
        }
    }

    // C += alpha A * B
    static inline void ProductFourthOrderTensor(double alpha, const Fourth_Order_Tensor& A, const Fourth_Order_Tensor& B, Fourth_Order_Tensor& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                    {
                        for(unsigned int m = 0; m < 3; ++m)
                        {
                            for(unsigned int n = 0; n < 3; ++n)
                            {
                                Result[i][j](k, l) += alpha * A[i][j](m, n) * B[m][n](k, l);
                            }
                        }
                    }
                }
            }
        }
    }

    // // inverse a fourth order tensor
    // REMARK: this function is likely incorrect
    // static inline void InverseFourthOrderTensor(const Fourth_Order_Tensor& A, Fourth_Order_Tensor& InvA)
    // {
    //     Matrix T(6, 6), InvT(6, 6);

    //     TensorToMatrix(A, T);
    //     InvertMatrix(T, InvT);
    //     MatrixToTensor(InvT, InvA);
    // }

    /**
     * Computes the derivatives of the inverse of matrix
     * d (A^-1) / dA = -A^(-1)_ki A^(-1)_lj
     * Reference:
     * + Le Khanh Chau's lecture note
     * + https://www.quora.com/What-is-the-derivative-of-inverse-matrix
     */
    template<typename TMatrixType>
    static void InverseDerivatives(const TMatrixType& InvA, Fourth_Order_Tensor& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                    {
                        Result[i][j](k, l) = -InvA(i, k) * InvA(l, j);
                    }
                }
            }
        }
    }

    /**
     * Computes the derivatives of the inverse of matrix and multiply with a factor
     * d (A^-1) / dA = -A^(-1)_ki A^(-1)_lj
     * Reference:
     * + Le Khanh Chau's lecture note
     * + https://www.quora.com/What-is-the-derivative-of-inverse-matrix
     */
    template<typename TMatrixType>
    static void AddInverseDerivatives(double alpha, const TMatrixType& InvA, Fourth_Order_Tensor& Result)
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                    {
                        Result[i][j](k, l) += -alpha*InvA(i, k) * InvA(l, j);
                    }
                }
            }
        }
    }

    /**
    * Generates the fourth order deviatoric unity tensor
    * @param Unity the deviatoric unity (will be overwritten)
    */
    static void DeviatoricUnity(std::vector<std::vector<MatrixType> >& Unity)
    {
        Unity.resize(3);

        MatrixType kronecker(3,3);
        noalias(kronecker)=ZeroMatrix(3,3);
        for(unsigned int i=0; i<3; i++)
        {
            kronecker(i,i)=1;
        }

        for(unsigned int i=0; i<3; i++)
        {
            Unity[i].resize(3);
            for(unsigned int j=0; j<3; j++)
            {
                Unity[i][j].resize(3,3,false);
                noalias(Unity[i][j])= ZeroMatrix(3,3);

                for(unsigned int k=0; k<3; k++)
                {
                    for(unsigned int l=0; l<3; l++)
                    {
                        Unity[i][j](k,l)=kronecker(i,k)*kronecker(j,l)
                                         -1.0/3.0*kronecker(i,j)*kronecker(k,l);
                    }
                }
            }
        }
    }

    /**
     * Generates the fourth order deviatoric unity tensor
     * @param Unity the deviatoric unity (will be overwritten)
     */
    static void DeviatoricUnity(array_1d<double,81>& Unity)
    {
        MatrixType kronecker(3,3);
        noalias(kronecker)=ZeroMatrix(3,3);
        for(unsigned int i=0; i<3; i++)
        {
            kronecker(i,i)=1;
        }

        for(unsigned int i=0; i<3; i++)
            for(unsigned int j=0; j<3; j++)
                for(unsigned int k=0; k<3; k++)
                    for(unsigned int l=0; l<3; l++)
                        Unity[27*i+9*j+3*k+l]=kronecker(i,k)*kronecker(j,l)
                                              -1.0/3.0*kronecker(i,j)*kronecker(k,l);
    }

    /**
    * Compute the elasticity tensor
    * @param C Elastic tensor
    * @param E Young modulus
    * @param NU Poisson ratio
    */
    static inline void CalculateElasticTensor( Fourth_Order_Tensor& C, double E, double NU )
    {
        MatrixType kronecker(3, 3);
        noalias(kronecker) = ZeroMatrix(3, 3);
        for(unsigned int i = 0; i < 3; ++i)
            kronecker(i, i) = 1.0;

        double lambda = NU * E / ((1 + NU) * (1 - 2 * NU));
        double mu     = E / (2 * (1 + NU));

        C.resize(3);
        for(unsigned int i = 0; i < 3; ++i)
        {
            C[i].resize(3);
            for(unsigned int j = 0; j < 3; ++j)
            {
                C[i][j].resize(3, 3);
                noalias(C[i][j]) = ZeroMatrix(3, 3);
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                        C[i][j](k,l) = lambda * kronecker(i, j) * kronecker(k, l)
                                     + mu * (kronecker(i, k) * kronecker(j, l)
                                           + kronecker(i, l) * kronecker(j, k));
                  }
             }
        }
    }

    // return A:B
    template<typename TMatrixType1, typename TMatrixType2>
    static double mat_inner_prod(const TMatrixType1& A, const TMatrixType2& B)
    {
        double res = 0.0;
        unsigned int m = A.size1();
        unsigned int n = A.size2();
        for(unsigned int i = 0; i < m; ++i)
            for(unsigned int j = 0; j < n; ++j)
                res += A(i, j) * B(i, j);
        return res;
    }

    // perform the Fortran vector product, i.e. element-by-element product (C[i] = A[i]*B[i])
    template<typename TVectorType1, typename TVectorType2>
    static TVectorType1 vec_fortran_prod(const TVectorType1& A, const TVectorType2& B)
    {
        TVectorType1 C = A;
        for (std::size_t i = 0; i < A.size(); ++i)
            C[i] *= B[i];
        return C;
    }

    // perform the Fortran matrix product, i.e. element-by-element product (C[i] = A[i]*B[i])
    template<typename TMatrixType1, typename TMatrixType2>
    static TMatrixType1 mat_fortran_prod(const TMatrixType1& A, const TMatrixType2& B)
    {
        TMatrixType1 C = A;
        for (std::size_t i = 0; i < A.size1(); ++i)
            for (std::size_t j = 0; j < A.size2(); ++j)
                C(i, j) *= B(i, j);
        return C;
    }

    template<typename TMatrixType>
    static inline double Trace(const TMatrixType& A)
    {
        double tr = 0.0;
        for(unsigned int i = 0; i < A.size1(); ++i)
            tr += A(i, i);
        return tr;
    }

    /*
     * Compute the isotropic function of the type
     *       Y(X) = sum{ y(x_i) E_i }
     * WHERE Y AND X ARE SYMMETRIC TENSORS, x_i AND E_i ARE, RESPECTIVELY
     * THE EIGENVALUES AND EIGENPROJECTIONS OF X, AND y(.) IS A SCALAR
     * FUNCTION.
     * X must be 3 x 3 matrix
     * Y will be resized accordingly
     * e must be sorted
     * Rererence: Section A.5.2, Computational Plasticity, de Souza Neto.
     */
    template<typename TMatrixType>
    static void ComputeIsotropicTensorFunction(
        unitary_func_t func,
        TMatrixType& Y,
        const std::vector<double>& e,
        const std::vector<Matrix>& eigprj
    )
    {
        if ((Y.size1() != 3) || (Y.size2() != 3))
            Y.resize(3, 3, false);

        Y.clear();
        for (int d = 0; d < 3; ++d)
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
     * e must be sorted
     * Rererence: Section A.5.2, Computational Plasticity, de Souza Neto.
     */
    template<typename TMatrixType>
    static void ComputeIsotropicTensorFunction(
        unitary_func_t func,
        TMatrixType& Y,
        const std::vector<double>& e,
        const std::vector<Matrix>& eigprj,
        double TOL // tolerance to compare eigenvalues
    )
    {
        if ((Y.size1() != 3) || (Y.size2() != 3))
            Y.resize(3, 3, false);

        Y.clear();
        if (fabs(e[0] - e[1]) > TOL && fabs(e[1] - e[2]) > TOL && fabs(e[0] - e[2]) > TOL)
        {
            for (int d = 0; d < 3; ++d)
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

    /*
     * Compute the derivative of the isotropic function of the type
     *       Y(X) = sum{ y(x_i) E_i }
     * WHERE Y AND X ARE SYMMETRIC TENSORS, x_i AND E_i ARE, RESPECTIVELY
      * THE EIGENVALUES AND EIGENPROJECTIONS OF X, AND y(.) IS A SCALAR
     * FUNCTION.
     * X must be 3 x 3 matrix
     * e and eigprj are the principle values and eigenprojections of X
     * dYdX will be resized accordingly
     * Rererence: Section A.5.2, Computational Plasticity, de Souza Neto.
     */
    template<typename TMatrixType>
    static void ComputeDerivativeIsotropicTensorFunction(
        unitary_func_t func,
        unitary_func_t dfunc,
        Fourth_Order_Tensor& dYdX,
        const TMatrixType& X,
        const std::vector<double>& e,
        const std::vector<TMatrixType>& eigprj,
        double TOL // tolerance to compare eigenvalues
    )
    {
        Fourth_Order_Tensor dX2dX;
        CalculateFourthOrderZeroTensor(dX2dX);
        TMatrixType I = IdentityMatrix(3);
        for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
                for(int k = 0; k < 3; ++k)
                    for(int l = 0; l < 3; ++l)
                        dX2dX[i][j](k, l) = 0.5 * (I(i, k) * X(l, j)
                                                 + I(i, l) * X(k, j)
                                                 + X(i, k) * I(j, l)
                                                 + X(i, l) * I(k, j));

        Fourth_Order_Tensor Is;
        CalculateFourthOrderSymmetricTensor(Is);

        CalculateFourthOrderZeroTensor(dYdX);

        if (fabs(e[0] - e[1]) > TOL && fabs(e[1] - e[2]) > TOL && fabs(e[0] - e[2]) > TOL)
        {
            int a, b, c;
            for (a = 0; a < 3; ++a)
            {
                if (a == 0) {b = 1; c = 2;}
                else if (a == 1) {b = 2; c = 0;}
                else if (a == 2) {b = 0; c = 1;}

                double aux1 = func(e[a]) / ((e[a] - e[b]) * (e[a] - e[c]));
                double aux2 = -aux1 * (e[b] + e[c]);
                double aux3 = -aux1 * (e[a] - e[b] + e[a] - e[c]);
                double aux4 = -aux1 * (e[b] - e[c]);

                AddFourthOrderTensor(aux1, dX2dX, dYdX);
                AddFourthOrderTensor(aux2, Is, dYdX);
                OuterProductFourthOrderTensor(aux3, eigprj[a], eigprj[a], dYdX);
                OuterProductFourthOrderTensor(aux4, eigprj[b], eigprj[b], dYdX);
                OuterProductFourthOrderTensor(-aux4, eigprj[c], eigprj[c], dYdX);
                OuterProductFourthOrderTensor(dfunc(e[a]), eigprj[a], eigprj[a], dYdX);
            }
        }
        else
        {
            if (fabs(e[0] - e[1]) < TOL && fabs(e[1] - e[2]) < TOL)
            {
                AddFourthOrderTensor(dfunc(e[0]), Is, dYdX);
            }
            else
            {
                const Matrix eye = IdentityMatrix(3);
                int a, c;

                if (fabs(e[0] - e[1]) < TOL)
                {
                    a = 2;
                    c = 0;
                }
                else // (fabs(e[1] - e[2]) < TOL)
                {
                    a = 0;
                    c = 2;
                }

                const double s1 = (func(e[a]) - func(e[c])) / pow(e[a] - e[c], 2) - dfunc(e[c]) / (e[a] - e[c]);
                const double s2 = 2.0 * e[c] * (func(e[a]) - func(e[c])) / pow(e[a] - e[c], 2) - dfunc(e[c]) * (e[a] + e[c]) / (e[a] - e[c]);
                const double s3 = 2.0 * (func(e[a]) - func(e[c])) / pow(e[a] - e[c], 3) - (dfunc(e[a]) + dfunc(e[c])) / pow(e[a] - e[c], 2);
                const double s4 = e[c] * s3;
                const double s5 = s4;
                const double s6 = e[c]*e[c] * s3;

                AddFourthOrderTensor(s1, dX2dX, dYdX);
                AddFourthOrderTensor(-s2, Is, dYdX);
                OuterProductFourthOrderTensor(-s3, X, X, dYdX);
                OuterProductFourthOrderTensor(s4, X, eye, dYdX);
                OuterProductFourthOrderTensor(s5, eye, X, dYdX);
                OuterProductFourthOrderTensor(-s6, eye, eye, dYdX);
            }
        }
    }

    /// Computation of the derivative of a general isotropic tensor function in three dimensions.
    /// Ref: Box A.6, Computational Plasticity, de Souza Neto.
    /// {ii, mm, jj} are the indices of the sorted pstress
    static inline void ComputeDerivativeGeneralIsotropicTensorFunction(
        int ii, int mm, int jj,
        const std::vector<double>& pstress,
        const std::vector<double>& pstrain,
        const MatrixType& dpstrs,
        const std::vector<MatrixType>& E,
        Fourth_Order_Tensor& Dep,
        double TOL = 1.0e-10
    )
    {
        MatrixType X = pstrain[0] * E[0] + pstrain[1] * E[1] + pstrain[2] * E[2];
        Fourth_Order_Tensor dX2dX;
        CalculateFourthOrderZeroTensor(dX2dX);
        MatrixType I = IdentityMatrix(3);
        for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
                for(int k = 0; k < 3; ++k)
                    for(int l = 0; l < 3; ++l)
                        dX2dX[i][j](k, l) = 0.5 * (I(i, k) * X(l, j)
                                                 + I(i, l) * X(k, j)
                                                 + X(i, k) * I(j, l)
                                                 + X(i, l) * I(k, j));

        #ifdef DEBUG_MOHR_COULOMB
        KRATOS_WATCH(pstrain[ii])
        KRATOS_WATCH(pstrain[mm])
        KRATOS_WATCH(pstrain[jj])
        #endif

        // box A.6
        Fourth_Order_Tensor Is;
        CalculateFourthOrderSymmetricTensor(Is);
        CalculateFourthOrderZeroTensor(Dep);
        if( ( fabs(pstrain[ii] - pstrain[mm]) > TOL )
         && ( fabs(pstrain[mm] - pstrain[jj]) > TOL )
         && ( fabs(pstrain[ii] - pstrain[jj]) > TOL ) )
        {
            #ifdef DEBUG_MOHR_COULOMB_WARNING_BOX_A6
            std::cout << "Box A.6 first case" << std::endl;
            #endif
            int cyc[] = {ii, mm, jj};
            for(int i0 = 0; i0 < 3; ++i0)
            {
                int i1 = (i0 == 0 ? 1 : (i0 == 1 ? 2 : 0)); // i1 = i0 + 1
                int i2 = (i0 == 0 ? 2 : (i0 == 1 ? 0 : 1)); // i2 = i0 + 2

                int a = cyc[i0];
                int b = cyc[i1];
                int c = cyc[i2];

                double aux1 = pstress[a] / (pstrain[a] - pstrain[b]) / (pstrain[a] - pstrain[c]);
                double aux2 = pstrain[b] + pstrain[c];
                double aux3 = pstrain[a] - pstrain[b] + pstrain[a] - pstrain[c];
                double aux4 = pstrain[b] - pstrain[c];
                AddFourthOrderTensor(aux1, dX2dX, Dep);
                AddFourthOrderTensor(-aux1 * aux2, Is, Dep);
                OuterProductFourthOrderTensor(-aux1 * aux3, E[a], E[a], Dep);
                OuterProductFourthOrderTensor(-aux1 * aux4, E[b], E[b], Dep);
                OuterProductFourthOrderTensor(aux1 * aux4, E[c], E[c], Dep);
            }

            for(int i = 0; i < 3; ++i)
            {
                for(int j = 0; j < 3; ++j)
                {
                    OuterProductFourthOrderTensor(dpstrs(i, j), E[i], E[j], Dep);
                }
            }
        }
        else
        {
            int a, b, c;
            bool isSame = false;
            if( ( fabs(pstrain[jj] - pstrain[mm]) < TOL ) && ( fabs(pstrain[ii] - pstrain[jj]) > TOL ))
            { a = ii; b = mm; c = jj; }
            else if( ( fabs(pstrain[jj] - pstrain[ii]) < TOL ) && ( fabs(pstrain[ii] - pstrain[mm]) > TOL ))
            { a = mm; b = ii; c = jj; }
            else if( ( fabs(pstrain[ii] - pstrain[mm]) < TOL ) && ( fabs(pstrain[jj] - pstrain[mm]) > TOL ))
            { a = jj; b = ii; c = mm; }
            else isSame = true;

            if(!isSame)
            {
                #ifdef DEBUG_MOHR_COULOMB_WARNING_BOX_A6
                std::cout << "Box A.6 second case" << std::endl;
                #endif
                double s1 = (pstress[a] - pstress[c]) / pow(pstrain[a] - pstrain[c], 2) + 1.0 / (pstrain[a] - pstrain[c]) * (dpstrs(c, b) - dpstrs(c, c));
                double s2 = 2 * pstrain[c] * (pstress[a] - pstress[c]) / pow(pstrain[a] - pstrain[c], 2) + (pstrain[a] + pstrain[c]) / (pstrain[a] - pstrain[c]) * (dpstrs(c, b) - dpstrs(c, c));
                double s3 = 2 * (pstress[a] - pstress[c]) / pow(pstrain[a] - pstrain[c], 3) + 1.0 / pow(pstrain[a] - pstrain[c], 2) * (dpstrs(a, c) + dpstrs(c, a) - dpstrs(a, a) - dpstrs(c, c));
                double s4 = 2 * pstrain[c] * (pstress[a] - pstress[c]) / pow(pstrain[a] - pstrain[c], 3) + 1.0 / (pstrain[a] - pstrain[c]) * (dpstrs(a, c) - dpstrs(c, b)) + pstrain[c] / pow(pstrain[a] - pstrain[c], 2) * (dpstrs(a, c) + dpstrs(c, a) - dpstrs(a, a) - dpstrs(c, c));
                double s5 = 2 * pstrain[c] * (pstress[a] - pstress[c]) / pow(pstrain[a] - pstrain[c], 3) + 1.0 / (pstrain[a] - pstrain[c]) * (dpstrs(c, a) - dpstrs(c, b)) + pstrain[c] / pow(pstrain[a] - pstrain[c], 2) * (dpstrs(a, c) + dpstrs(c, a) - dpstrs(a, a) - dpstrs(c, c));
                double s6 = 2 * pow(pstrain[c], 2) * (pstress[a] - pstress[c]) / pow(pstrain[a] - pstrain[c], 3) + pstrain[a] * pstrain[c] / pow(pstrain[a] - pstrain[c], 2) * (dpstrs(a, c) + dpstrs(c, a)) - pow(pstrain[c], 2) / pow(pstrain[a] - pstrain[c], 2) * (dpstrs(a, a) + dpstrs(c, c)) - (pstrain[a] + pstrain[c]) / (pstrain[a] - pstrain[c]) * dpstrs(c, b);
                AddFourthOrderTensor(s1, dX2dX, Dep);
                AddFourthOrderTensor(-s2, Is, Dep);
                OuterProductFourthOrderTensor(-s3, X, X, Dep);
                OuterProductFourthOrderTensor(s4, X, I, Dep);
                OuterProductFourthOrderTensor(s5, I, X, Dep);
                OuterProductFourthOrderTensor(-s6, I, I, Dep);
            }
            else
            {
                #ifdef DEBUG_MOHR_COULOMB_WARNING_BOX_A6
                std::cout << "Box A.6 third case" << std::endl;
                #endif
                AddFourthOrderTensor(dpstrs(ii, ii) - dpstrs(ii, mm), Is, Dep);
                OuterProductFourthOrderTensor(dpstrs(ii, mm), I, I, Dep);
            }
        }
    }

    /**
    * Performs clipping on the two polygons clipping_points and subjected_points (the technique used i
    * Sutherland-Hodgman clipping) and returns the overlapping polygon result_points. The method works
    * in 3D. Both polygons have to be convex, but they can be slightly perturbated in 3D space, this
    * allows for performing clipping on two interpolated interfaces
    * @param clipping_points vertices of clipping polygon
    * @param subjected_points vertices of subjected polygon
    * @param result_points vertices of overlapping polygon
    * @return false= no overlapping polygon, true= overlapping polygon found
    */
    static bool Clipping(std::vector<PointType* >& clipping_points,std::vector<PointType* >& subjected_points, std::vector<PointType* >& result_points, double tolerance)
    {
        result_points= subjected_points;
        VectorType actual_edge(3);
        VectorType actual_normal(3);
        std::vector<PointType* > temp_results;
        bool is_visible= false;
        for(unsigned int clipp_edge=0; clipp_edge<clipping_points.size(); clipp_edge++)
        {
            temp_results.clear();
            unsigned int    index_clipp_2=0;
            if(clipp_edge< (clipping_points.size()-1))
                index_clipp_2= clipp_edge+1;
            //define clipping edge vector
            noalias(actual_edge)= *(clipping_points[clipp_edge])-*(clipping_points[index_clipp_2]);
            noalias(actual_edge)= actual_edge/sqrt(inner_prod(actual_edge,actual_edge));

            //define normal on clipping-edge vector towards visible side
            if(clipp_edge< (clipping_points.size()-2))
                actual_normal=*(clipping_points[clipp_edge+2])-(*(clipping_points[clipp_edge])+actual_edge*inner_prod((*(clipping_points[clipp_edge+2])-*(clipping_points[clipp_edge])),actual_edge));
            else if(clipp_edge< (clipping_points.size()-1))
                actual_normal=*(clipping_points[0])-(*(clipping_points[clipp_edge])+actual_edge*inner_prod((*(clipping_points[0])-*(clipping_points[clipp_edge])),actual_edge));
            else
                actual_normal=*(clipping_points[1])-(*(clipping_points[clipp_edge])+actual_edge*inner_prod((*(clipping_points[1])-*(clipping_points[clipp_edge])),actual_edge));

            noalias(actual_normal)=actual_normal/(sqrt(inner_prod(actual_normal,actual_normal)));

            //test if the first point is visible or unvisible
            if( inner_prod((*(result_points[0])-*(clipping_points[clipp_edge])), actual_normal)> tolerance)
            {
                is_visible= true;
            }
            else
            {
                is_visible= false;
            }

            for(unsigned int subj_edge=0; subj_edge< result_points.size(); subj_edge++)
            {
                unsigned int    index_subj_2=0;

                if(subj_edge< (result_points.size()-1))
                    index_subj_2= subj_edge+1;

                //Test whether the points of the actual subj_edge lay on clipp_edge
                if(fabs(inner_prod((*(result_points[subj_edge])-(*(clipping_points[clipp_edge]))), actual_normal))<= tolerance)
                {
                    temp_results.push_back(result_points[subj_edge]);
                    if( inner_prod((*(result_points[index_subj_2])-*(clipping_points[clipp_edge])), actual_normal)> tolerance)
                        is_visible= true;
                    else
                        is_visible= false;

                    continue;
                }
                //Calculate minimal distance between the two points
                VectorType b(2);
                b(0)= -inner_prod((*(result_points[index_subj_2])-*(result_points[subj_edge])),(*(result_points[subj_edge])-*(clipping_points[clipp_edge])));
                b(1)= inner_prod((*(clipping_points[index_clipp_2])-*(clipping_points[clipp_edge])),(*(result_points[subj_edge])-*(clipping_points[clipp_edge])));
                MatrixType A(2,2);
                A(0,0)=inner_prod((*(result_points[index_subj_2])-*(result_points[subj_edge])),(*(result_points[index_subj_2])-*(result_points[subj_edge])));
                A(0,1)=-inner_prod((*(result_points[index_subj_2])-*(result_points[subj_edge])),(*(clipping_points[index_clipp_2])-*(clipping_points[clipp_edge])));
                A(1,0)= A(0,1);
                A(1,1)=inner_prod(*(clipping_points[index_clipp_2])-*(clipping_points[clipp_edge]),*(clipping_points[index_clipp_2])-*(clipping_points[clipp_edge]));
                VectorType coeff(2);
                coeff(0)=1.0/A(0,0)*(b(0)-A(0,1)/(A(1,1)-A(0,1)*A(1,0)/A(0,0))*(b(1)-b(0)*A(1,0)/A(0,0)));
                coeff(1)=1.0/(A(1,1)-A(0,1)*A(1,0)/A(0,0))*(b(1)-b(0)*A(1,0)/A(0,0));


                //TEST on distance to endpoints of the line
                VectorType dist_vec(3);
                noalias(dist_vec)= *(result_points[subj_edge])+coeff(0)*(*(result_points[index_subj_2])-*(result_points[subj_edge]))-(*(clipping_points[clipp_edge])+coeff(1)*(*(clipping_points[index_clipp_2])-*(clipping_points[clipp_edge])));

                if( coeff(0) > tolerance && coeff(0) < (1-tolerance)&& (sqrt(inner_prod(dist_vec,dist_vec))< tolerance))
                {
                    if(is_visible)
                        temp_results.push_back(result_points[subj_edge]);

                    temp_results.push_back(new PointType((*(clipping_points[clipp_edge])+coeff(1)*(*(clipping_points[index_clipp_2])-*(clipping_points[clipp_edge])))));

                    is_visible= !is_visible;

                    continue;
                }
                if(is_visible)
                    temp_results.push_back(result_points[subj_edge]);
            }
            result_points=temp_results;
        }
        if(result_points.size()==0)
            return false;
        else
            return true;
    }

    /**
     * Solve a1 * x + b1 * y = c1
     *       a2 * x + b2 * y = c2
     */
    static void Solve2Unknowns(
        const TDataType& a1, const TDataType& b1, const TDataType& c1,
        const TDataType& a2, const TDataType& b2, const TDataType& c2,
        TDataType& x, TDataType& y
    )
    {
        x = (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
        y = (c2 * a1 - c1 * a2) / (a1 * b2 - a2 * b1);
    }

    /**
     * Solve a1 * x + b1 * y + c1 * z = d1
     *       a2 * x + b2 * y + c2 * z = d2
     *       a3 * x + b3 * y + c3 * z = d3
     */
    static void Solve3Unknowns(
        const TDataType& a1, const TDataType& b1, const TDataType& c1, const TDataType& d1,
        const TDataType& a2, const TDataType& b2, const TDataType& c2, const TDataType& d2,
        const TDataType& a3, const TDataType& b3, const TDataType& c3, const TDataType& d3,
        TDataType& x, TDataType& y, TDataType& z
    )
    {
        TDataType DetA = Det3(a1, b1, c1, a2, b2, c2, a3, b3, c3);
        x = Det3(d1, b1, c1, d2, b2, c2, d3, b3, c3) / DetA;
        y = Det3(a1, d1, c1, a2, d2, c2, a3, d3, c3) / DetA;
        z = Det3(a1, b1, d1, a2, b2, d2, a3, b3, d3) / DetA;
    }

    /*
     * Solve the quadratic equation a*x^2 + b*x + c = 0
     * In return:   0: no solution
     *              1: one solution
     *              2: two solutions
     */
    static int SolveQuadratic(
        const TDataType& a, const TDataType& b, const TDataType& c,
        TDataType x[2]
    )
    {
        TDataType b2 = b / 2;
        TDataType delta = b2 * b2 - a * c;

        if(delta < 0.0)
        {
            return 0; // no solution
        }
        else if(delta == 0.00)
        {
            if(a == 0.0)
            {
                return 0; // no solution
            }
            else
            {
                x[0] = -b2 / a;
                x[1] = x[0];
                return 1; // one solution
            }
        }
        else
        {
            if(a == 0.0)
            {
                x[0] = - c / b;
                x[1] = x[0];
                return 1; // one solution
            }
            else
            {
                x[0] = (-b2 - sqrt(delta)) / a;
                x[1] = (-b2 + sqrt(delta)) / a;
                return 2; // two solution
            }
        }

        return 0;
    }

    /**
     * Calculate the second derivatives providing the Jacobian and the local gradients
     * REMARKS: the arrangement is
     * in 2D: [d^2/dx^2 d^2/dy^2 d^2/dxdy]
     * in 3D: [d^2/dx^2 d^2/dy^2 d^2/dz^2 d^2/dxdy d^2/dydz d^2/dxdz]
     */
    template<class TGeometryType>
    static void CalculateSecondDerivatives(const TGeometryType& rGeometry, std::vector<Vector>& rD2N_DX2,
        const Matrix& J, const Matrix& DN_DX, const typename TGeometryType::CoordinatesArrayType& rPoint)
    {
        const unsigned int dim = rGeometry.WorkingSpaceDimension();
        const unsigned int number_of_nodes = rGeometry.size();

        // compute the shape function local second derivatives
        typename TGeometryType::ShapeFunctionsSecondDerivativesType D2N_De2;
        D2N_De2 = rGeometry.ShapeFunctionsSecondDerivatives(D2N_De2, rPoint);

        // compute the second derivatives of physical coordinates
        Matrix D2X_De2(dim, dim);
        Matrix D2Y_De2(dim, dim);
        Matrix D2Z_De2(dim, dim);
        for(unsigned int i = 0; i < dim; ++i)
        {
            for(unsigned int j = 0; j < dim; ++j)
            {
                D2X_De2(i, j) = 0.0;
                for(unsigned int k = 0; k < number_of_nodes; ++k)
                    D2X_De2(i, j) += rGeometry[k].X0() * D2N_De2[k](i, j);

                D2Y_De2(i, j) = 0.0;
                for(unsigned int k = 0; k < number_of_nodes; ++k)
                    D2Y_De2(i, j) += rGeometry[k].Y0() * D2N_De2[k](i, j);

                if(dim > 2)
                {
                    D2Z_De2(i, j) = 0.0;
                    for(unsigned int k = 0; k < number_of_nodes; ++k)
                        D2Z_De2(i, j) += rGeometry[k].Z0() * D2N_De2[k](i, j);
                }
            }
        }

        // compute the second derivatives in physical space
        unsigned int mat_size = dim * (dim + 1) / 2;
        Matrix D2(mat_size, mat_size);
        if(dim == 2)
        {
            D2(0, 0) = pow(J(0, 0), 2);
            D2(0, 1) = pow(J(0, 1), 2);
            D2(0, 2) = J(0, 0) * J(0, 1);

            D2(1, 0) = pow(J(1, 0), 2);
            D2(1, 1) = pow(J(1, 1), 2);
            D2(1, 2) = J(1, 0) * J(1, 1);

            D2(2, 0) = 2.0 * J(0, 0) * J(1, 0);
            D2(2, 1) = 2.0 * J(0, 1) * J(1, 1);
            D2(2, 2) = J(0, 0) * J(1, 1) + J(1, 0) * J(0, 1);
        }
        else if(dim == 3)
        {
            D2(0, 0) = pow(J(0, 0), 2);
            D2(0, 1) = pow(J(0, 1), 2);
            D2(0, 2) = pow(J(0, 2), 2);
            D2(0, 3) = J(0, 0) * J(0, 1);
            D2(0, 4) = J(0, 1) * J(0, 2);
            D2(0, 5) = J(0, 0) * J(0, 2);

            D2(1, 0) = pow(J(1, 0), 2);
            D2(1, 1) = pow(J(1, 1), 2);
            D2(1, 2) = pow(J(1, 2), 2);
            D2(1, 3) = J(1, 0) * J(1, 1);
            D2(1, 4) = J(1, 1) * J(1, 2);
            D2(1, 5) = J(1, 0) * J(1, 2);

            D2(2, 0) = pow(J(2, 0), 2);
            D2(2, 1) = pow(J(2, 1), 2);
            D2(2, 2) = pow(J(2, 2), 2);
            D2(2, 3) = J(2, 0) * J(2, 1);
            D2(2, 4) = J(2, 1) * J(2, 2);
            D2(2, 5) = J(2, 0) * J(2, 2);

            D2(3, 0) = 2.0 * J(0, 0) * J(1, 0);
            D2(3, 1) = 2.0 * J(0, 1) * J(1, 1);
            D2(3, 2) = 2.0 * J(0, 2) * J(1, 2);
            D2(3, 3) = J(0, 0) * J(1, 1) + J(0, 1) * J(1, 0);
            D2(3, 4) = J(0, 1) * J(1, 2) + J(0, 2) * J(1, 1);
            D2(3, 5) = J(0, 0) * J(1, 2) + J(0, 2) * J(1, 0);

            D2(4, 0) = 2.0 * J(1, 0) * J(2, 0);
            D2(4, 1) = 2.0 * J(1, 1) * J(2, 1);
            D2(4, 2) = 2.0 * J(1, 2) * J(2, 2);
            D2(4, 3) = J(1, 0) * J(2, 1) + J(1, 1) * J(2, 0);
            D2(4, 4) = J(1, 1) * J(2, 2) + J(1, 2) * J(2, 1);
            D2(4, 5) = J(1, 0) * J(2, 2) + J(1, 2) * J(2, 0);

            D2(5, 0) = 2.0 * J(0, 0) * J(2, 0);
            D2(5, 1) = 2.0 * J(0, 1) * J(2, 1);
            D2(5, 2) = 2.0 * J(0, 2) * J(2, 2);
            D2(5, 3) = J(0, 0) * J(2, 1) + J(0, 1) * J(2, 1);
            D2(5, 4) = J(0, 1) * J(2, 2) + J(0, 2) * J(2, 1);
            D2(5, 5) = J(0, 0) * J(2, 2) + J(0, 2) * J(2, 0);
        }

        Matrix InvD2(mat_size, mat_size);
        InvertMatrix(D2, InvD2);

        if(rD2N_DX2.size() != number_of_nodes)
            rD2N_DX2.resize(number_of_nodes);

        Vector R2(mat_size);
        for(unsigned int i = 0; i < number_of_nodes; ++i)
        {
            if(rD2N_DX2[i].size() != mat_size)
                rD2N_DX2[i].resize(mat_size, false);

            if(dim == 2)
            {
                R2(0) = D2N_De2[i](0, 0) - DN_DX(i, 0) * D2X_De2(0, 0) - DN_DX(i, 1) * D2Y_De2(0, 0);
                R2(1) = D2N_De2[i](1, 1) - DN_DX(i, 0) * D2X_De2(1, 1) - DN_DX(i, 1) * D2Y_De2(1, 1);
                R2(2) = D2N_De2[i](0, 1) - DN_DX(i, 0) * D2X_De2(0, 1) - DN_DX(i, 1) * D2Y_De2(0, 1);
            }
            else if(dim == 3)
            {
                R2(0) = D2N_De2[i](0, 0) - DN_DX(i, 0) * D2X_De2(0, 0) - DN_DX(i, 1) * D2Y_De2(0, 0) - DN_DX(i, 2) * D2Z_De2(0, 0);
                R2(1) = D2N_De2[i](1, 1) - DN_DX(i, 0) * D2X_De2(1, 1) - DN_DX(i, 1) * D2Y_De2(1, 1) - DN_DX(i, 2) * D2Z_De2(1, 1);
                R2(2) = D2N_De2[i](2, 2) - DN_DX(i, 0) * D2X_De2(2, 2) - DN_DX(i, 1) * D2Y_De2(2, 2) - DN_DX(i, 2) * D2Z_De2(2, 2);
                R2(3) = D2N_De2[i](0, 1) - DN_DX(i, 0) * D2X_De2(0, 1) - DN_DX(i, 1) * D2Y_De2(0, 1) - DN_DX(i, 2) * D2Z_De2(0, 1);
                R2(4) = D2N_De2[i](1, 2) - DN_DX(i, 0) * D2X_De2(1, 2) - DN_DX(i, 1) * D2Y_De2(1, 2) - DN_DX(i, 2) * D2Z_De2(1, 2);
                R2(5) = D2N_De2[i](0, 2) - DN_DX(i, 0) * D2X_De2(0, 2) - DN_DX(i, 1) * D2Y_De2(0, 2) - DN_DX(i, 2) * D2Z_De2(0, 2);
            }

            noalias(rD2N_DX2[i]) = prod(R2, InvD2);
        }
    }

    /**
     * Calculate the second derivatives of the geometry at one point
     * REMARKS: the arrangement is
     * in 2D: [d^2/dx^2 d^2/dy^2 d^2/dxdy]
     * in 3D: [d^2/dx^2 d^2/dy^2 d^2/dz^2 d^2/dxdy d^2/dydz d^2/dxdz]
     */
    template<class TGeometryType, std::size_t Configuration>
    static void CalculateSecondDerivatives(
        const TGeometryType& rGeometry,
        std::vector<Vector>& rD2N_DX2,
        const typename TGeometryType::CoordinatesArrayType& rPoint)
    {
        unsigned int dim = rGeometry.WorkingSpaceDimension();
        unsigned int number_of_nodes = rGeometry.size();

        // compute the Jacobian
        Matrix J;
        if(Configuration == 0) // compute the Jacobian in undeformed configuration
        {
            Matrix DeltaPosition(rGeometry.size(), 3);
            for ( unsigned int node = 0; node < rGeometry.size(); ++node )
                noalias( row( DeltaPosition, node ) ) = rGeometry[node].Coordinates() - rGeometry[node].GetInitialPosition();
            J = rGeometry.Jacobian( J, rPoint, DeltaPosition );
        }
        else if(Configuration == 1) // compute the Jacobian in deformed configuration
        {
            J = rGeometry.Jacobian(J, rPoint);
        }

        // compute inverse of Jacobian
        Matrix InvJ;
        double DetJ;
        MathUtils<double>::InvertMatrix(J, InvJ, DetJ);

        // compute the shape function local gradients
        Matrix DN_De;
        DN_De = rGeometry.ShapeFunctionsLocalGradients(DN_De, rPoint);

        // compute the shape function gradients w.r.t physical coordinates
        Matrix DN_DX(number_of_nodes, dim);
        noalias(DN_DX) = prod(DN_De, InvJ);

        // compute the second derivatives
        CalculateSecondDerivatives<TGeometryType>(rGeometry, rD2N_DX2, J, DN_DX, rPoint);
    }

    /**
     * Assemble a sub-vector to a big vector
     * num*dim is the size of sub-vector
     * stride is the shift within dim of the big vector
     * begin is where to assemble the values in the big vector. This can be used to shift the chunk.
     * Example:
     *  Let say the sub-vector is [1 1.1 1.2 2 2.1 2.2 ...]
     *  dim in this case is 3
     *  with stride = 2 and begin = 0, the big vector will look like
     *    [1 1.1 1.2 0 0 2 2.1 2.2 0 0 ....]
     */
    static void AssembleVectorFromSubVector( VectorType& rRightHandSideVector,
            const unsigned int& dim,
            const unsigned int& num,
            const unsigned int& stride,
            const unsigned int& begin,
            const Vector& R_P )
    {
        for ( unsigned int prim = 0; prim < num; ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                rRightHandSideVector( begin + prim*(dim+stride) + i ) += R_P( prim*dim + i );
            }
        }
    }

    /**
     * Assemble a sub-matrix to a big matrix
     * The index notation in each dimension is similar to AssembleVectorFromSubVector
     */
    static void AssembleMatrixFromSubMatrix( MatrixType& rLeftHandSideMatrix,
            const unsigned int& dim1,
            const unsigned int& num1,
            const unsigned int& stride1,
            const unsigned int& begin1,
            const unsigned int& dim2,
            const unsigned int& num2,
            const unsigned int& stride2,
            const unsigned int& begin2,
            const Matrix& K_PQ)
    {
        for ( unsigned int prim = 0; prim < num1; ++prim )
        {
            for ( unsigned int i = 0; i < dim1; ++i )
            {
                for ( unsigned int sec = 0; sec < num2; ++sec )
                {
                    for ( unsigned int j = 0; j < dim2; ++j )
                    {
                        rLeftHandSideMatrix( begin1 + prim*(dim1+stride1) + i, begin2 + sec*(dim2+stride2) + j )
                                += K_PQ( prim*dim1 + i, sec*dim2 + j );
                    }
                }
            }
        }
    }

};// class SD_MathUtils

}

#endif /* SD_MATH_UTILS defined */
