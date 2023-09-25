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
/* *********************************************************
*
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2009-01-14 17:14:12 $
*   Revision:            $Revision: 1.13 $
*
* ***********************************************************/

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws/isotropic_3d.h"
#include "structural_application_variables.h"

namespace Kratos
{


/**
 * TO BE TESTED!!!
 */
Isotropic3D::Isotropic3D()
    : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
Isotropic3D::~Isotropic3D()
{
}

bool Isotropic3D::Has( const Variable<int>& rThisVariable )
{
    return false;
}

bool Isotropic3D::Has( const Variable<double>& rThisVariable )
{
    if ( rThisVariable == PRESTRESS_FACTOR || rThisVariable == YOUNG_MODULUS || rThisVariable == POISSON_RATIO )
        return true;

    return false;
}

bool Isotropic3D::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == INSITU_STRESS )
        return true;

    if ( rThisVariable == PRESTRESS )
        return true;

    if ( rThisVariable == STRESSES )
        return true;

    return false;
}

bool Isotropic3D::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}


int& Isotropic3D::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    if (rThisVariable == IS_SHAPE_FUNCTION_REQUIRED)
        rValue = 0;

    return rValue;
}

double& Isotropic3D::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
    {
        rValue = mPrestressFactor;
        return rValue;
    }

    if(rThisVariable == YOUNG_MODULUS )
    {
       rValue = mE;
       return rValue;
    }

    if ( rThisVariable == POISSON_RATIO )
    {
        rValue = mNU;
        return rValue;
    }

    if(rThisVariable == DAMAGE)
    {
        rValue = 0.00;
        return rValue;
    }

    if (rThisVariable == DELTA_TIME)
    {
        rValue = sqrt(mE/mDE);
        return rValue;
    }

    if (rThisVariable == PRESSURE_P)
    {
        rValue = -(mCurrentStress[0] + mCurrentStress[1] + mCurrentStress[2]) / 3;
        return rValue;
    }

    if (rThisVariable == PRESSURE_Q || rThisVariable == VON_MISES_STRESS)
    {
        double p = (mCurrentStress[0] + mCurrentStress[1] + mCurrentStress[2]) / 3;
        double sxx = mCurrentStress[0] - p;
        double syy = mCurrentStress[1] - p;
        double szz = mCurrentStress[2] - p;
        double sxy = mCurrentStress[3];
        double syz = mCurrentStress[4];
        double sxz = mCurrentStress[5];

        rValue = sqrt( 3.0 * ( 0.5*(sxx*sxx + syy*syy + szz*szz) + sxy*sxy + syz*syz + sxz*sxz ) );
        return rValue;
    }

    rValue = 0.00;
    return rValue;
}

Vector& Isotropic3D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        const unsigned int size = mPrestress.size();
        if(rValue.size() != size)
            rValue.resize(size, false );
        noalias(rValue) = mPrestressFactor * mPrestress;
        return rValue;
    }

    if ( rThisVariable == STRESSES || rThisVariable == THREED_STRESSES )
    {
        const unsigned int size = mCurrentStress.size();
        rValue.resize(size, false );
        noalias(rValue)  = mCurrentStress;
        return rValue;
    }

    if ( rThisVariable == THREED_STRAIN )
    {
        // REF: http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_stress.cfm
        const unsigned int size = mCurrentStress.size();
        rValue.resize(size, false );
        Vector Stress = mCurrentStress + mPrestressFactor * mPrestress;
        rValue(0) = (Stress(0) - mNU*Stress(1) - mNU*Stress(2)) / mE;
        rValue(1) = (Stress(1) - mNU*Stress(0) - mNU*Stress(2)) / mE;
        rValue(2) = (Stress(2) - mNU*Stress(0) - mNU*Stress(1)) / mE;
        rValue(3) = 2.0*(1.0+mNU)*Stress(3) / mE;
        rValue(4) = 2.0*(1.0+mNU)*Stress(4) / mE;
        rValue(5) = 2.0*(1.0+mNU)*Stress(5) / mE;
        return rValue;
    }

    if ( rThisVariable == PLASTIC_STRAIN_VECTOR )
    {
        rValue = ZeroVector( 6 );
        return( rValue );
    }

    if ( rThisVariable == INTERNAL_VARIABLES )
    {
        rValue = ZeroVector( 1 );
        return( rValue );
    }

    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}

Matrix& Isotropic3D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
        KRATOS_THROW_ERROR( std::logic_error, "Matrix Variable case not considered", "" );
}

void Isotropic3D::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void Isotropic3D::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
    if ( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;
}

void Isotropic3D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void Isotropic3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        noalias(mPrestress) = rValue;
    }
    else if ( rThisVariable == STRESSES || rThisVariable == INITIAL_STRESS )
    {
        noalias(mCurrentStress) = rValue;
    }
}

void Isotropic3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}


void Isotropic3D::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mCurrentStress = ZeroVector( 6 );
    mPrestress = ZeroVector( 6 );
    mPrestressFactor = 1.0;
    mE = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mDE = props[DENSITY];
}

void Isotropic3D::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    noalias(mCurrentStress) = mPrestressFactor*mPrestress;
}

void Isotropic3D::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void Isotropic3D::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void Isotropic3D::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void Isotropic3D::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void Isotropic3D::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        if(AlgorithmicTangent.size1() != 6 || AlgorithmicTangent.size2() != 6)
            AlgorithmicTangent.resize(6, 6, false);
        CalculateElasticMatrix( AlgorithmicTangent, mE, mNU );
    }

    if (rValues.IsSetStressVector())
    {
        const Vector& StrainVector = rValues.GetStrainVector();
        Vector& StressVector = rValues.GetStressVector();
        if(StressVector.size() != 6)
            StressVector.resize(6, false);

        if (rValues.IsSetConstitutiveMatrix())
        {
            Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
            CalculateStress( StrainVector, AlgorithmicTangent, StressVector );
        }
        else
        {
            Matrix AlgorithmicTangent(6, 6);
            CalculateElasticMatrix( AlgorithmicTangent, mE, mNU );
            CalculateStress( StrainVector, AlgorithmicTangent, StressVector );
        }
    }
}

void Isotropic3D::CalculateMaterialResponsePK2 (Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void Isotropic3D::CalculateMaterialResponse( const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo,
        const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        bool CalculateStresses,
        int CalculateTangent,
        bool SaveInternalVariables )
{
    ConstitutiveLaw::Parameters const_params;
    Vector ThisStrainVector = StrainVector;
    const_params.SetStrainVector(ThisStrainVector);
    if (CalculateStresses)
        const_params.SetStressVector(StressVector);
    if (CalculateTangent)
        const_params.SetConstitutiveMatrix(AlgorithmicTangent);

    this->CalculateMaterialResponseCauchy(const_params);
}

void Isotropic3D::CalculateElasticMatrix( Matrix& C, const double& E, const double& NU )
{
    //setting up material matrix
    double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
    double c2 = c1 * ( 1 - NU );
    double c3 = c1 * NU;
    double c4 = c1 * 0.5 * ( 1 - 2 * NU );
    //filling material matrix
    noalias(C) = ZeroMatrix(6, 6);
    C( 0, 0 ) = c2;
    C( 0, 1 ) = c3;
    C( 0, 2 ) = c3;
    C( 1, 0 ) = c3;
    C( 1, 1 ) = c2;
    C( 1, 2 ) = c3;
    C( 2, 0 ) = c3;
    C( 2, 1 ) = c3;
    C( 2, 2 ) = c2;
    C( 3, 3 ) = c4;
    C( 4, 4 ) = c4;
    C( 5, 5 ) = c4;
}

/**
 * TO BE TESTED!!!
 */
void Isotropic3D::CalculateStress( const Vector& StrainVector, Matrix& AlgorithmicTangent, Vector& StressVector )
{
    if ( StressVector.size() != 6 )
    {
        StressVector.resize( 6 );
    }

    noalias( StressVector ) = prod( AlgorithmicTangent, StrainVector ) - mPrestressFactor * mPrestress;

    noalias(mCurrentStress) = StressVector;
}

/**
 * TO BE TESTED!!!
 */
void Isotropic3D::CalculateStress( const double& E, const double& NU, const Vector& StrainVector, Vector& StressVector ) const
{
    if ( StressVector.size() != 6 )
    {
        StressVector.resize( 6 );
    }

    Matrix Ce( 6, 6 );

    this->CalculateElasticMatrix( Ce, E, NU );

    noalias( StressVector ) = prod( Ce, StrainVector );
}


//**********************************************************************
void Isotropic3D::CalculateCauchyStresses(
    Vector& rCauchy_StressVector,
    const Matrix& rF,
    const Vector& rPK2_StressVector,
    const Vector& rGreenLagrangeStrainVector )
{
    Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

    double J = MathUtils<double>::Det3( rF );
    boost::numeric::ublas::bounded_matrix<double, 3, 3> mstemp;
    boost::numeric::ublas::bounded_matrix<double, 3, 3> msaux;

    noalias( mstemp ) = prod( rF, S );
    noalias( msaux ) = prod( mstemp, trans( rF ) );
    msaux *= J;

    if ( rCauchy_StressVector.size() != 6 )
        rCauchy_StressVector.resize( 6 );

    rCauchy_StressVector[0] = msaux( 0, 0 );

    rCauchy_StressVector[1] = msaux( 1, 1 );

    rCauchy_StressVector[2] = msaux( 2, 2 );

    rCauchy_StressVector[3] = msaux( 0, 1 );

    rCauchy_StressVector[4] = msaux( 0, 2 );

    rCauchy_StressVector[5] = msaux( 1, 2 );
}

//**********************************************************************
int Isotropic3D::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY

    if ( !props.Has( YOUNG_MODULUS ) || !props.Has(POISSON_RATIO) )
    {
        KRATOS_THROW_ERROR( std::logic_error, "this constitutive law requires YOUNG_MODULUS and POISSON_RATIO given as KRATOS variables", "" );
    }

    double nu = props[POISSON_RATIO];

    if ( nu > 0.499 && nu < 0.501 )
    {
        KRATOS_THROW_ERROR( std::logic_error, "invalid poisson ratio in input, close to incompressibility", "" );
        return -1;
    }

    return 0;

    KRATOS_CATCH( "" );
}
} // Namespace Kratos
