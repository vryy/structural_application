/*
LICENSE: see soil_mechanics_application/LICENSE.txt
*/
/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 8 Aug 2018 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes

#include "constitutive_laws/neo_hookean_2d.h"
#include "utilities/math_utils.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "structural_application.h"

namespace Kratos
{

//**********************************************************************
NeoHookean2D::NeoHookean2D()
    : ConstitutiveLaw()
{
}

//**********************************************************************
NeoHookean2D::~NeoHookean2D()
{
}

//**********************************************************************
bool NeoHookean2D::Has( const Variable<int>& rThisVariable )
{
    return false;
}

//**********************************************************************
bool NeoHookean2D::Has( const Variable<double>& rThisVariable )
{
    if ( rThisVariable == PRESTRESS_FACTOR || rThisVariable == YOUNG_MODULUS || rThisVariable == POISSON_RATIO )
        return true;

    return false;
}

//**********************************************************************
bool NeoHookean2D::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == INSITU_STRESS )
        return true;

    if ( rThisVariable == PRESTRESS )
        return true;

    if ( rThisVariable == STRESSES )
        return true;

    return false;
}

//**********************************************************************
bool NeoHookean2D::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

//**********************************************************************
int& NeoHookean2D::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return rValue;
}

//**********************************************************************
double& NeoHookean2D::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if ( rThisVariable == PRESTRESS_FACTOR ){
        rValue = mPrestressFactor;
    return rValue;
    }
    if(rThisVariable == YOUNG_MODULUS ){
       rValue = mE;
       return rValue;
     }


    if ( rThisVariable == POISSON_RATIO ){
        rValue = mNU;
        return rValue;
    }

    if(rThisVariable==DAMAGE){
        rValue = 0.00;
        return rValue;
    }

    if (rThisVariable==DELTA_TIME){
        rValue = sqrt(mE/mDE);
        return rValue;
    }

    rValue = 0.00;
    return rValue;
}

//**********************************************************************
Vector& NeoHookean2D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        const unsigned int size = mPrestress.size();
        if(rValue.size() != size)
            rValue.resize(size, false);
        noalias(rValue) = mPrestressFactor * mPrestress;
        return rValue;
    }

    if ( rThisVariable == STRESSES )
    {
        const unsigned int size = mCurrentStress.size();
        if(rValue.size() != size)
            rValue.resize(size, false);
        noalias(rValue)  = mCurrentStress;
        return rValue;
    }

    if ( rThisVariable == PLASTIC_STRAIN_VECTOR )
    {
        noalias(rValue) = ZeroVector( rValue.size() );
        return( rValue );
    }

    if ( rThisVariable == INTERNAL_VARIABLES )
    {
        rValue = ZeroVector( 1 );
        return( rValue );
    }

    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}

//**********************************************************************
Matrix& NeoHookean2D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
        KRATOS_THROW_ERROR( std::logic_error, "Matrix Variable case not considered", "" );
}

//**********************************************************************
void NeoHookean2D::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

//**********************************************************************
void NeoHookean2D::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
    if ( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;
}

//**********************************************************************
void NeoHookean2D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

//**********************************************************************
void NeoHookean2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        noalias(mPrestress) = rValue;
    }
}

//**********************************************************************
void NeoHookean2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

//**********************************************************************
void NeoHookean2D::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mCurrentStress = ZeroVector( 3 );
    mPrestress = ZeroVector( 3 );
    mPrestressFactor = 1.0;
    mE = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mDE = props[DENSITY];
}

//**********************************************************************
void NeoHookean2D::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    mPrestress = ZeroVector( 3 );
    mPrestressFactor = 1.0;
}

//**********************************************************************
void NeoHookean2D::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

//**********************************************************************
void NeoHookean2D::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

//**********************************************************************
void NeoHookean2D::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

//**********************************************************************
void NeoHookean2D::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

//**********************************************************************
void  NeoHookean2D::CalculateMaterialResponse( const Vector& StrainVector,
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
    CalculateStress( StressVector, StrainVector );
    CalculateTangentMatrix( AlgorithmicTangent, StrainVector );
}

//**********************************************************************
void NeoHookean2D::CalculateTangentMatrix( Matrix& C, const Vector& StrainVector )
{
    if ( C.size1() != 3 || C.size2() != 3 )
    {
        C.resize( 3, 3 );
    }

    double e_11 = StrainVector(0);
    double e_22 = StrainVector(1);
    double e_12 = StrainVector(2);

    double mu = mE / (2 * (1.0 + mNU));
    double lambda = mE * mNU / ((1.0 + mNU) * (1.0 - 2*mNU));

    double aux = 4*e_11*e_22+ 2*e_11 - pow(e_12, 2) + 2*e_22 + 1;

    C( 0, 0 ) = pow(2*e_22 + 1, 2)*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 0, 1 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 1)*(2*e_22 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 0, 2 ) = e_12*(2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);

    C( 1, 0 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 1)*(2*e_22 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 1, 1 ) = pow(2*e_11 + 1, 2)*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 1, 2 ) = e_12*(2*e_11 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);

    C( 2, 0 ) = e_12*(2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 2, 1 ) = e_12*(2*e_11 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 2, 2 ) = (2*pow(e_12, 2)*(-lambda*log(aux) + lambda + 2*mu) - lambda*(aux)*log(aux) + 2*mu*(aux))/(2*pow(aux, 2));
}

//**********************************************************************
void NeoHookean2D::CalculateStress( Vector& StressVector, const Vector& StrainVector )
{
    if ( StressVector.size() != 3 )
    {
        StressVector.resize( 3 );
    }

    double mu = mE / (2 * (1.0 + mNU));
    double lambda = mE * mNU / ((1.0 + mNU) * (1.0 - 2*mNU));

    double e_11 = StrainVector(0);
    double e_22 = StrainVector(1);
    double e_12 = StrainVector(2);

    double aux = 4*e_11*e_22 + 2*e_11 - pow(e_12, 2) + 2*e_22 + 1;

    StressVector(0) = (lambda*(2*e_22 + 1)*log(aux)/2 - mu*(2*e_22 + 1) + mu*aux) / aux;

    StressVector(1) = (lambda*(2*e_11 + 1)*log(aux)/2 - mu*(2*e_11 + 1) + mu*aux) / aux;

    StressVector(2) = e_12*(-lambda*log(aux) + 2*mu) / (2*aux);

    noalias( StressVector ) -= mPrestressFactor * mPrestress;

    noalias(mCurrentStress) = StressVector;
}

//**********************************************************************
void NeoHookean2D::CalculateCauchyStresses(
    Vector& rCauchy_StressVector,
    const Matrix& rF,
    const Vector& rPK2_StressVector,
    const Vector& rGreenLagrangeStrainVector )
{
    Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

    double J = MathUtils<double>::Det2( rF );
    boost::numeric::ublas::bounded_matrix<double, 2, 2> mstemp;
    boost::numeric::ublas::bounded_matrix<double, 2, 2> msaux;

    noalias( mstemp ) = prod( rF, S );
    noalias( msaux ) = prod( mstemp, trans( rF ) );
    msaux *= J;

    if ( rCauchy_StressVector.size() != 6 )
        rCauchy_StressVector.resize( 6 );

    rCauchy_StressVector[0] = msaux( 0, 0 );

    rCauchy_StressVector[1] = msaux( 1, 1 );

    rCauchy_StressVector[2] = msaux( 0, 1 );
}

//**********************************************************************
int NeoHookean2D::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
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
