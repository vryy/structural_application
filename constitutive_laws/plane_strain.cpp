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
*   Last Modified by:    $Author: kazem $
*   Date:                $Date: 2009-01-16 10:50:04 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws/plane_strain.h"
#include "structural_application_variables.h"

namespace Kratos
{
/**
 * TO BE TESTED!!!
 */
PlaneStrain::PlaneStrain()
    : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
PlaneStrain::~PlaneStrain()
{
}

bool PlaneStrain::Has(const Variable<int>& rThisVariable)
{
    return false;
}

bool PlaneStrain::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool PlaneStrain::Has( const Variable<Vector>& rThisVariable )
{
    if(rThisVariable == STRESSES)
        return true;

    return false;
}

bool PlaneStrain::Has( const Variable<Matrix>& rThisVariable )
{
    if(rThisVariable == ALGORITHMIC_TANGENT)
        return true;
    return false;
}

int& PlaneStrain::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    if (rThisVariable == IS_SHAPE_FUNCTION_REQUIRED)
        rValue = 0;

    return rValue;
}

double& PlaneStrain::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if(rThisVariable==DELTA_TIME)
        rValue = sqrt(mE/mDE);

    if(rThisVariable == PRESTRESS_FACTOR)
        rValue = mPrestressFactor;

    if(rThisVariable == YOUNG_MODULUS)
        rValue = mE;

    if(rThisVariable == POISSON_RATIO)
        rValue = mNU;

    if (rThisVariable == PRESSURE_P)
    {
        double o_zz = mNU * (mCurrentStress[0] + mCurrentStress[1]);
        rValue = -(mCurrentStress[0] + mCurrentStress[1] + o_zz) / 3;
        return rValue;
    }

    if (rThisVariable == PRESSURE_Q || rThisVariable == VON_MISES_STRESS)
    {
        double o_zz = mNU * (mCurrentStress[0] + mCurrentStress[1]);
        double p = (mCurrentStress[0] + mCurrentStress[1] + o_zz) / 3;
        double sxx = mCurrentStress[0] - p;
        double syy = mCurrentStress[1] - p;
        double szz = o_zz - p;
        double sxy = mCurrentStress[2];
        double syz = 0.0;
        double sxz = 0.0;

        rValue = sqrt( 3.0 * ( 0.5*(sxx*sxx + syy*syy + szz*szz) + sxy*sxy + syz*syz + sxz*sxz ) );
        return rValue;
    }

    return rValue;
}

Vector& PlaneStrain::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == STRESSES || rThisVariable == STRESSES_OLD )
    {
        rValue = mCurrentStress;
    }

    if ( rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS )
    {
        rValue = mPreStress;
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
    }

    if ( rThisVariable == THREED_PRESTRESS )
    {
        if(rValue.size() != 6)
            rValue.resize(6, false);
        rValue(0) = mPreStress(0);
        rValue(1) = mPreStress(1);
        rValue(2) = mNU * (mPreStress(0) + mPreStress(1));
        rValue(3) = mPreStress(2);
        rValue(4) = 0.0;
        rValue(5) = 0.0;
    }

    // REF: http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_strain.cfm
    if ( rThisVariable == THREED_STRAIN )
    {
        if(rValue.size() != 6)
            rValue.resize(6, false);
        double aux = (1.0+mNU)/mE;
        rValue(0) = aux * ((1.0-mNU)*mCurrentStress(0) - mNU*mCurrentStress(1));
        rValue(1) = aux * ((1.0-mNU)*mCurrentStress(1) - mNU*mCurrentStress(0));
        rValue(2) = 0.0;
        rValue(3) = 2.0*aux*mCurrentStress(2);
        rValue(4) = 0.0;
        rValue(5) = 0.0;
    }

    return rValue;
}

Matrix& PlaneStrain::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    if (rThisVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
    {
        for(unsigned int i = 0; i< rValue.size2(); ++i)
            rValue(0, i) = 0.00;
    }

    if(rThisVariable == ALGORITHMIC_TANGENT || rThisVariable == ELASTIC_TANGENT)
    {
        if(rValue.size1() != 3 || rValue.size2() != 3)
            rValue.resize(3, 3, false);
        CalculateElasticMatrix( rValue, mE, mNU );
    }

    return( rValue );
}

void PlaneStrain::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void PlaneStrain::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
    if ( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;
}

void PlaneStrain::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS )
    {
        noalias(mPreStress) = rValue;
    }
    if ( rThisVariable == STRESSES || rThisVariable == INITIAL_STRESS )
    {
        if(mCurrentStress.size() != rValue.size())
            mCurrentStress.resize(rValue.size(), false);
        noalias(mCurrentStress) = rValue;
    }
    if ( rThisVariable == THREED_STRESSES )
    {
        mCurrentStress(0) = rValue(0);
        mCurrentStress(1) = rValue(1);
        mCurrentStress(2) = rValue(3);
    }
}

void PlaneStrain::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void PlaneStrain::Calculate( const Variable<Matrix>& rVariable, Matrix& rResult,
                             const ProcessInfo& rCurrentProcessInfo )
{
}

void PlaneStrain::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    if (rValues.IsSetStressVector())
    {
        const Vector& StrainVector = rValues.GetStrainVector();
        Vector& StressVector = rValues.GetStressVector();

        if(StressVector.size() != 3)
            StressVector.resize(3, false);
        CalculateStress( StrainVector, StressVector );
    }

    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        if(AlgorithmicTangent.size1() != 3 || AlgorithmicTangent.size2() != 3)
            AlgorithmicTangent.resize(3, 3, false);
        CalculateConstitutiveMatrix( AlgorithmicTangent );
    }
}

void PlaneStrain::CalculateMaterialResponsePK2 (Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void PlaneStrain::CalculateMaterialResponse( const Vector& StrainVector,
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
    const_params.SetStressVector(StressVector);
    const_params.SetConstitutiveMatrix(AlgorithmicTangent);

    this->CalculateMaterialResponseCauchy(const_params);
}

/**
 * TO BE TESTED!!!
 */
void PlaneStrain::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mCurrentStress = ZeroVector( 3 );
    mPreStress = ZeroVector( 3 );
    mPrestressFactor = 1.0;
    mE  = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mDE = props[DENSITY];
}

void PlaneStrain::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    noalias(mCurrentStress) = ZeroVector(3);
}

void PlaneStrain::InitializeNonLinearIteration( const Properties& rMaterialProperties,
                                                const GeometryType& rElementGeometry,
                                                const Vector& rShapeFunctionsValues,
                                                const ProcessInfo& rCurrentProcessInfo )
{
}

void PlaneStrain::FinalizeNonLinearIteration( const Properties& rMaterialProperties,
                                              const GeometryType& rElementGeometry,
                                              const Vector& rShapeFunctionsValues,
                                              const ProcessInfo& rCurrentProcessInfo )
{
}

/**
 * TO BE TESTED!!!
 */
void PlaneStrain::CalculateElasticMatrix( Matrix& C, const double& E, const double& NU ) const
{
    double c1 = E * ( 1.00 - NU ) / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
    double c2 = E * NU / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
    double c3 = 0.5 * E / ( 1 + NU );

    C( 0, 0 ) = c1;
    C( 0, 1 ) = c2;
    C( 0, 2 ) = 0.0;
    C( 1, 0 ) = c2;
    C( 1, 1 ) = c1;
    C( 1, 2 ) = 0.0;
    C( 2, 0 ) = 0.0;
    C( 2, 1 ) = 0.0;
    C( 2, 2 ) = c3;
}

/**
 * TO BE TESTED!!!
 */
void PlaneStrain::CalculateStress( const Vector& StrainVector, Vector& StressVector )
{
    double c1 = mE * ( 1.00 - mNU ) / (( 1.00 + mNU ) * ( 1.00 - 2 * mNU ) );
    double c2 = mE * mNU / (( 1.00 + mNU ) * ( 1.00 - 2 * mNU ) );
    double c3 = 0.5 * mE / ( 1 + mNU );

    // compute the stress based on strain
    StressVector[0] = c1 * StrainVector[0] + c2 * StrainVector[1];
    StressVector[1] = c1 * StrainVector[1] + c2 * StrainVector[0];
    StressVector[2] = c3 * StrainVector[2];

    noalias(StressVector) -= mPrestressFactor*mPreStress;

    noalias( mCurrentStress ) = StressVector;
}

/**
 * TO BE TESTED!!!
 */
void PlaneStrain::CalculateStress( const double& E, const double& NU, const Vector& StrainVector, Vector& StressVector ) const
{
    double c1 = E * ( 1.00 - NU ) / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
    double c2 = E * NU / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
    double c3 = 0.5 * E / ( 1 + NU );

    // compute the stress based on strain
    StressVector[0] = c1 * StrainVector[0] + c2 * StrainVector[1];
    StressVector[1] = c1 * StrainVector[1] + c2 * StrainVector[0];
    StressVector[2] = c3 * StrainVector[2];
}


/**
 * TO BE REVIEWED!!!
 */
void PlaneStrain::CalculateConstitutiveMatrix( Matrix& rResult ) const
{
    CalculateElasticMatrix( rResult, mE, mNU );
}


std::size_t PlaneStrain::GetStrainSize() const
{
    return 3;
}

//**********************************************************************
//**********************************************************************

void PlaneStrain::CalculateCauchyStresses(
    Vector& rCauchy_StressVector,
    const Matrix& rF,
    const Vector& rPK2_StressVector,
    const Vector& rGreenLagrangeStrainVector )
{
    Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

    double J = MathUtils<double>::Det2( rF );

    boost::numeric::ublas::bounded_matrix<double, 2, 2> temp;
    boost::numeric::ublas::bounded_matrix<double, 2, 2> aux;

    noalias( temp ) = prod( rF, S );
    noalias( aux ) = prod( temp, trans( rF ) );
    aux *= J;

    if ( rCauchy_StressVector.size() != 3 )
        rCauchy_StressVector.resize( 3 );

    rCauchy_StressVector[0] = aux( 0, 0 );

    rCauchy_StressVector[1] = aux( 1, 1 );

    rCauchy_StressVector[2] = aux( 0,1 );
}

//**********************************************************************

int PlaneStrain::Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo) const
{
    if(YOUNG_MODULUS.Key() == 0 || props[YOUNG_MODULUS]<= 0.00)
        KRATOS_THROW_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");

    const double& nu = props[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );
    if(POISSON_RATIO.Key() == 0 || check==true) // props[POISSON_RATIO] == 1.00 || props[POISSON_RATIO] == -1.00)
        KRATOS_THROW_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");

    if(DENSITY.Key() == 0 || props[DENSITY]<0.00)
        KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");

    return 0;

}

} // Namespace Kratos
