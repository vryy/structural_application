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
#include "utilities/math_utils.h"
#include "constitutive_laws/neo_hookean_2d.h"
#include "structural_application_variables.h"

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
    if (rThisVariable == IS_SHAPE_FUNCTION_REQUIRED)
        rValue = 0;
    if (rThisVariable == PARENT_ELEMENT_ID)
        rValue = mElemId;
    if (rThisVariable == INTEGRATION_POINT_INDEX)
        rValue = mGaussId;

    return rValue;
}

//**********************************************************************
double& NeoHookean2D::GetValue( const Variable<double>& rThisVariable, double& rValue )
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

    if(rThisVariable==DAMAGE)
    {
        rValue = 0.00;
        return rValue;
    }

    if (rThisVariable==DELTA_TIME)
    {
        rValue = sqrt(mE/mDE);
        return rValue;
    }

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

    if ( rThisVariable == THREED_STRESSES )
    {
        if(rValue.size() != 6)
            rValue.resize(6, false);
        rValue(0) = mCurrentStress(0);
        rValue(1) = mCurrentStress(1);
        rValue(2) = mNU * (mCurrentStress(0) + mCurrentStress(1));
        rValue(3) = mCurrentStress(2);
        rValue(4) = 0.0;
        rValue(5) = 0.0;
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

    return rValue;
}

//**********************************************************************
Matrix& NeoHookean2D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return rValue;
}

//**********************************************************************
void NeoHookean2D::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == PARENT_ELEMENT_ID)
        mElemId = rValue;
    if (rThisVariable == INTEGRATION_POINT_INDEX)
        mGaussId = rValue;
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
void NeoHookean2D::CalculateMaterialResponsePK2(Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();
    Vector& StressVector = rValues.GetStressVector();
    Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();

    const int cmap[] = {0, 1, 3};

    Vector StrainVector3D(6);
    noalias(StrainVector3D) = ZeroVector(6);

    for (int i = 0; i < 3; ++i)
        StrainVector3D(cmap[i]) = StrainVector(i);

    if (rValues.IsSetStressVector())
    {
        Vector StressVector3D(6);

        CalculateStress( StressVector3D, StrainVector3D );

        for (int i = 0; i < 3; ++i)
            StressVector(i) = StressVector3D(cmap[i]);
    }

    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix AlgorithmicTangent3D(6, 6);
        CalculateTangentMatrix( AlgorithmicTangent3D, StrainVector3D );

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                AlgorithmicTangent(i, j) = AlgorithmicTangent3D(cmap[i], cmap[j]);
            }
        }
    }
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
    Vector StrainVector3D(6), StressVector3D(6);
    Matrix AlgorithmicTangent3D(6, 6);

    noalias(StrainVector3D) = ZeroVector(6);
    StrainVector3D(0) = StrainVector(0);
    StrainVector3D(1) = StrainVector(1);
    StrainVector3D(3) = StrainVector(2);

    CalculateStress( StressVector3D, StrainVector3D );
    CalculateTangentMatrix( AlgorithmicTangent3D, StrainVector3D );

    StressVector(0) = StressVector3D(0);
    StressVector(1) = StressVector3D(1);
    StressVector(2) = StressVector3D(3);

    std::vector<int> cmap = {0, 1, 3};
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            AlgorithmicTangent(i, j) = AlgorithmicTangent3D(cmap[i], cmap[j]);
        }
    }
}

//**********************************************************************
void NeoHookean2D::CalculateTangentMatrix( Matrix& C, const Vector& StrainVector )
{
    if ( C.size1() != 6 || C.size2() != 6 )
    {
        C.resize( 6, 6 );
    }

    double e_11 = StrainVector(0);
    double e_22 = StrainVector(1);
    double e_33 = StrainVector(2);
    double e_12 = StrainVector(3);
    double e_23 = StrainVector(4);
    double e_13 = StrainVector(5);

    double mu = mE / (2 * (1.0 + mNU));
    double lambda = mE * mNU / ((1.0 + mNU) * (1.0 - 2*mNU));

    double aux = 4*e_11*e_22 + 4*e_11*e_33 + 2*e_11 - pow(e_12, 2) - pow(e_13, 2) + 4*e_22*e_33 + 2*e_22 - pow(e_23, 2) + 2*e_33 + 1;

    C( 0, 0 ) = pow(2*e_22 + 2*e_33 + 1, 2)*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 0, 1 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_33 + 1)*(2*e_22 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 0, 2 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_22 + 1)*(2*e_22 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 0, 3 ) = e_12*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 0, 4 ) = e_23*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 0, 5 ) = e_13*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);

    C( 1, 0 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_33 + 1)*(2*e_22 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 1, 1 ) = pow(2*e_11 + 2*e_33 + 1, 2)*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 1, 2 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_22 + 1)*(2*e_11 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 1, 3 ) = e_12*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 1, 4 ) = e_23*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 1, 5 ) = e_13*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);

    C( 2, 0 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_22 + 1)*(2*e_22 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 2, 1 ) = ((lambda*log(aux) - 2*mu)*(aux) + (2*e_11 + 2*e_22 + 1)*(2*e_11 + 2*e_33 + 1)*(-lambda*log(aux) + lambda + 2*mu))/pow(aux, 2);
    C( 2, 2 ) = pow(2*e_11 + 2*e_22 + 1, 2)*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 2, 3 ) = e_12*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 2, 4 ) = e_23*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 2, 5 ) = e_13*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);

    C( 3, 0 ) = e_12*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 3, 1 ) = e_12*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 3, 2 ) = e_12*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 3, 3 ) = (2*pow(e_12, 2)*(-lambda*log(aux) + lambda + 2*mu) - lambda*(aux)*log(aux) + 2*mu*(aux))/(2*pow(aux, 2));
    C( 3, 4 ) = e_12*e_23*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 3, 5 ) = e_12*e_13*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);

    C( 4, 0 ) = e_23*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 4, 1 ) = e_23*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 4, 2 ) = e_23*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 4, 3 ) = e_12*e_23*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 4, 4 ) = (2*pow(e_23, 2)*(-lambda*log(aux) + lambda + 2*mu) - lambda*(aux)*log(aux) + 2*mu*(aux))/(2*pow(aux, 2));
    C( 4, 5 ) = e_13*e_23*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);

    C( 5, 0 ) = e_13*(2*e_22 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 5, 1 ) = e_13*(2*e_11 + 2*e_33 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 5, 2 ) = e_13*(2*e_11 + 2*e_22 + 1)*(lambda*log(aux) - lambda - 2*mu)/pow(aux, 2);
    C( 5, 3 ) = e_12*e_13*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 5, 4 ) = e_13*e_23*(-lambda*log(aux) + lambda + 2*mu)/pow(aux, 2);
    C( 5, 5 ) = (2*pow(e_13, 2)*(-lambda*log(aux) + lambda + 2*mu) - lambda*(aux)*log(aux) + 2*mu*(aux))/(2*pow(aux, 2));
}

//**********************************************************************
void NeoHookean2D::CalculateStress( Vector& StressVector, const Vector& StrainVector )
{
    if ( StressVector.size() != 6 )
    {
        StressVector.resize( 6 );
    }

    double e_11 = StrainVector(0);
    double e_22 = StrainVector(1);
    double e_33 = StrainVector(2);
    double e_12 = StrainVector(3);
    double e_23 = StrainVector(4);
    double e_13 = StrainVector(5);

    double mu = mE / (2 * (1.0 + mNU));
    double lambda = mE * mNU / ((1.0 + mNU) * (1.0 - 2*mNU));


    double aux = 4*e_11*e_22 + 4*e_11*e_33 + 2*e_11 - pow(e_12, 2) - pow(e_13, 2) + 4*e_22*e_33 + 2*e_22 - pow(e_23, 2) + 2*e_33 + 1;

    if (aux < 0)
    {
        KRATOS_WATCH(mE)
        KRATOS_WATCH(mNU)
        KRATOS_WATCH(StrainVector)
        KRATOS_WATCH(aux)
        KRATOS_THROW_ERROR(std::logic_error, "Error encounting NaN values", "")
    }

    StressVector(0) = (lambda*(2*e_22 + 2*e_33 + 1)*log(aux)/2 - mu*(2*e_22 + 2*e_33 + 1) + mu*aux) / aux;

    StressVector(1) = (lambda*(2*e_11 + 2*e_33 + 1)*log(aux)/2 - mu*(2*e_11 + 2*e_33 + 1) + mu*aux) / aux;

    StressVector(2) = (lambda*(2*e_11 + 2*e_22 + 1)*log(aux)/2 - mu*(2*e_11 + 2*e_22 + 1) + mu*aux) / aux;

    StressVector(3) = e_12*(-lambda*log(aux) + 2*mu) / (2*aux);

    StressVector(4) = e_23*(-lambda*log(aux) + 2*mu) / (2*aux);

    StressVector(5) = e_13*(-lambda*log(aux) + 2*mu) / (2*aux);

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
int NeoHookean2D::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
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
