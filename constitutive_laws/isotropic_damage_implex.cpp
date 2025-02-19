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
*   Last Modified by:    $Author: seyedali $
*   Date:                $Date: 2014-05-15 00:00:00 $
*   Revision:            $Revision: 0.54 $
*
* ***********************************************************/

// Remark: with reference to J. Oliver et al. (2008) and the
// Java version of IMPL-EX algorithm by hbui

// System includes
#include<cmath>
#include <iostream>

// External includes

// Project includes

#include "utilities/math_utils.h"
#include "constitutive_laws/isotropic_damage_implex.h"
#include "structural_application_variables.h"

namespace Kratos
{

/**
 * TO BE TESTED!!!
 */
IsotropicDamageIMPLEX::IsotropicDamageIMPLEX()
    : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
IsotropicDamageIMPLEX::~IsotropicDamageIMPLEX()
{
}

bool IsotropicDamageIMPLEX::Has( const Variable<int>& rThisVariable )
{
    return false;
}

bool IsotropicDamageIMPLEX::Has( const Variable<double>& rThisVariable )
{
    if( rThisVariable == YOUNG_MODULUS || rThisVariable == POISSON_RATIO || rThisVariable == DAMAGE || rThisVariable == ALPHA)
        return true;
    return false;
}

bool IsotropicDamageIMPLEX::Has( const Variable<Vector>& rThisVariable )
{
    if( rThisVariable == STRESSES )
        return true;
    return false;
}

bool IsotropicDamageIMPLEX::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

double& IsotropicDamageIMPLEX::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if( rThisVariable == YOUNG_MODULUS )
        rValue = mE;
    if( rThisVariable == POISSON_RATIO )
        rValue = mNU;
    if( rThisVariable == DAMAGE )
        rValue = mD;
    if( rThisVariable == ALPHA )
        rValue = mAlpha;
    if( rThisVariable == EQUIVALENT_STRAIN )
    {
        // Compute elastic constitutive tensor
        const unsigned int strain_size = mCurrentStrain.size();
        Matrix C = ZeroMatrix(strain_size, strain_size);
        CalculateElasticMatrix(C, mE, mNU);

        // Compute effective stresses
        Vector EffectiveStressVector = prod(C, mCurrentStrain);

        // Equivalent strain based on the energy norm
        rValue = sqrt((1.0/mE)*inner_prod(trans(mCurrentStrain), EffectiveStressVector));
    }

    return rValue;
}

Vector& IsotropicDamageIMPLEX::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if( rThisVariable == STRESSES )
        noalias(rValue) = mCurrentStress;

    return rValue;
}

Matrix& IsotropicDamageIMPLEX::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return rValue;
}

void IsotropicDamageIMPLEX::SetValue( const Variable<bool>& rThisVariable, const bool& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == FIX_DAMAGE )
    {
        if (rValue)
            mDamageFlag = -1;
        else
            mDamageFlag = 0;
    }
}

void IsotropicDamageIMPLEX::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                                      const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == PARENT_ELEMENT_ID )
        mElemId = rValue;
    if( rThisVariable == INTEGRATION_POINT_INDEX )
        mGaussId = rValue;
}

void IsotropicDamageIMPLEX::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                      const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if( rThisVariable == POISSON_RATIO )
        mNU = rValue;
    if( rThisVariable == DAMAGE )
    {
        mD = rValue;
        mAlpha = this->ComputeKappa(mD);
        mAlpha_old = mAlpha;
        mAlpha_old_old = mAlpha;
        mq = mAlpha;
    }
    if( rThisVariable == INITIAL_DAMAGE )
    {
        mInitialDamage = rValue;
        mD = mInitialDamage;
        mAlpha = this->ComputeKappa(mD);
        mAlpha_old = mAlpha;
        mAlpha_old_old = mAlpha;
    }
    if( rThisVariable == EQUIVALENT_STRAIN )
    {
        mInitialEps = rValue;
    }
}

void IsotropicDamageIMPLEX::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                                      const array_1d<double, 3>& rValue,
                                      const ProcessInfo& rCurrentProcessInfo )
{
}

void IsotropicDamageIMPLEX::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                      const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == STRESSES )
        noalias(mCurrentStress) = rValue;
}

void IsotropicDamageIMPLEX::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                      const ProcessInfo& rCurrentProcessInfo )
{
}

void IsotropicDamageIMPLEX::InitializeMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = dim * (dim+1) / 2;

    mCurrentStress = ZeroVector(strain_size);
    mCurrentStrain = ZeroVector(strain_size);
    mE  = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mE_0 = props[DAMAGE_E0];
    mE_f = props[DAMAGE_EF];
    mD   = 0.0;
    mAlpha = mE_0;
    mq     = mAlpha;
    mAlpha_old_old = mAlpha_old;
    mAlpha_old     = mAlpha;

    mInitialDamage = 0.0;
    mInitialEps = 0.0;
    mDamageFlag = 0;

    mDeltaTime_old = 1.0;
    mDeltaTime = 1.0;
}

void IsotropicDamageIMPLEX::ResetMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    mD = mInitialDamage;
    mAlpha = mE_0;
    mq     = mAlpha;
    mAlpha_old_old = mAlpha_old;
    mAlpha_old     = mAlpha;
}

void IsotropicDamageIMPLEX::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        const ProcessInfo& CurrentProcessInfo )
{
    // Initialize algorithmic IMPL-EX variables
    mAlpha_old_old = mAlpha_old;
    mAlpha_old = mAlpha;
    mq_old = mq;

//     KRATOS_WATCH(mAlpha_old_old);
//     KRATOS_WATCH(mAlpha_old);
//     KRATOS_WATCH(mAlpha);

    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = dim * (dim+1) / 2;

    //==========================================================================================
    //-----------------------------------{ EXPLICIT ALGORITHM }---------------------------------
    //==========================================================================================

    // Compute elastic constitutive tensor
    Matrix C = ZeroMatrix(strain_size, strain_size);
    CalculateElasticMatrix(C, mE, mNU);

    // Explicit extrapolation of the internal variable alpha
    mdAlpha = mAlpha_old - mAlpha_old_old;
    mAlpha_alg = mAlpha_old + mdAlpha;

    // Calculate softening function
    double H = SofteningLaw(mAlpha_old);

    // Compute stress-like internal variable
    mq_alg = mq_old + H*mdAlpha;

    // Compute extrapolated damage variable
    mDamage_alg = mInitialDamage + (1.0 - mInitialDamage) * ( 1.0 - ( mq_alg / mAlpha_alg ) );

    // Compute algorithmic tangent operator
    mC_alg = (1.0 - mDamage_alg)*C;
}

void IsotropicDamageIMPLEX::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        const ProcessInfo& CurrentProcessInfo )
{
    mDeltaTime_old = CurrentProcessInfo[DELTA_TIME];
    if (std::abs(mDeltaTime_old) < 1e-13) mDeltaTime_old = 1.0;
}

void IsotropicDamageIMPLEX::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void IsotropicDamageIMPLEX::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        const ProcessInfo& CurrentProcessInfo )
{
    // KRATOS_WATCH("NONLINEAR ITERATION ACTIVATED------------------------------------------------------------------------");
}

void IsotropicDamageIMPLEX::CalculateMaterialResponseCauchy (Parameters& parameters)
{
    this->CalculateMaterialResponse( parameters.GetStrainVector()
        , parameters.GetDeformationGradientF()
        , parameters.GetStressVector()
        , parameters.GetConstitutiveMatrix()
        , parameters.GetProcessInfo()
        , parameters.GetMaterialProperties()
        , parameters.GetElementGeometry()
        , parameters.GetShapeFunctionsValues()
        , parameters.IsSetStressVector()
        , parameters.IsSetConstitutiveMatrix()
        , true
    );
}

void IsotropicDamageIMPLEX::CalculateMaterialResponse( const Vector& StrainVector,
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
    if (CurrentProcessInfo[SET_CALCULATE_REACTION])
    {
        if (CalculateStresses)
        {
            noalias(StressVector) = mCurrentStress;
            return;
        }
    }

    const unsigned int strain_size = StrainVector.size();

    double TOL = 1.0e-8;
    if (props.Has(LOCAL_ERROR_TOLERANCE))
        TOL = props[LOCAL_ERROR_TOLERANCE];

    //==========================================================================================
    //----------------------------------{ IMPLICIT ALGORITHM }----------------------------------
    //==========================================================================================

    // store the current strain
    noalias(mCurrentStrain) = StrainVector;

    // Compute elastic constitutive tensor
    Matrix C = ZeroMatrix(strain_size, strain_size);
    CalculateElasticMatrix(C, mE, mNU);

    // Compute effective stresses
    Vector EffectiveStressVector = prod(C, StrainVector);

    // Equivalent strain based on the energy norm
    double eps_eq = sqrt((1.0/mE)*inner_prod(trans(StrainVector), EffectiveStressVector)) - mInitialEps;

    // Compute trial state
    double alpha_trial = mAlpha_old;

    // Determine loading function
    double f = eps_eq - alpha_trial;

    // Check loading function and compute damage multiplier
    double dlambda;
    if (mDamageFlag != -1)
    {
        if( f > TOL )
        {
            // Damage/loading state
            dlambda = f;
            mDamageFlag = 1;
            // std::cout << "Damage is updated at element " << mElemId << ", point " << mGaussId << std::endl;
        }
        else
        {
            // Elastic/unloading state
            dlambda = 0.0;
            mDamageFlag = 0;
        }
    }
    else
    {
        dlambda = 0.0; // forcing steady damage
    }

//     KRATOS_WATCH(dlambda);

    // Update internal variable
    mAlpha = mAlpha_old + dlambda;

    // Calculate softening function
    double H = SofteningLaw(mAlpha_old);

    // Update stress-like internal variable
    mq = mq_old + H*dlambda;

    // Update damage variable
    mD = mInitialDamage + (1.0 - mInitialDamage) * (1.0 - (mq/mAlpha));

//     KRATOS_WATCH(mD);
//     KRATOS_WATCH(mq);
//     KRATOS_WATCH(mAlpha);

//     if (mD > 1.0)
//     {
//        KRATOS_WATCH("DAMAGE VARIABLE IS NOT VALID.");
//        KRATOS_WATCH(mq/mAlpha);
//     }
//
//     if (mD < 0)
//        KRATOS_WATCH("DAMAGE VARIABLE IS NEGATIVE ...");

    // Update corrected stresses
    Matrix C_t = (1.0 - mD)*C;
    mCurrentStress = prod( C_t, StrainVector );

    //==========================================================================================
    //----------------------------------{ EXPLICIT ALGORITHM }----------------------------------
    //==========================================================================================

    mDeltaTime = CurrentProcessInfo[DELTA_TIME];

    // Explicit extrapolation of the internal variable alpha
    mdAlpha = mAlpha_old - mAlpha_old_old;
    mAlpha_alg = mAlpha_old + mDeltaTime / mDeltaTime_old * mdAlpha;

    // Compute stress-like internal variable
    mq_alg = mq_old + H*mdAlpha;

    // Compute extrapolated damage variable
    mDamage_alg = mInitialDamage + (1.0 - mInitialDamage) * ( 1.0 - ( mq_alg / mAlpha_alg ) );

    if (std::abs(mDamage_alg - 1.0) < 1e-3)
    {
        mDamage_alg = 0.999;
        std::cout << "Damage is 0.999 at element " << mElemId << ", point " << mGaussId << std::endl;
    }

    // Compute algorithmic tangent operator
    mC_alg = (1.0 - mDamage_alg)*C;

    if (CalculateTangent)
    {
        // Set algorithmic tangent operator
        noalias(AlgorithmicTangent) = mC_alg;
    }

    // Compute extrapolated stresses
    if (CalculateStresses)
        noalias(StressVector) = prod( mC_alg, StrainVector );
}

void IsotropicDamageIMPLEX::CalculateElasticMatrix( Matrix& C, const double E, const double NU ) const
{
    // Compute elastic constitutive tensor
    if (C.size1() == 6) // 3D
    {
        double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
        double c2 = c1 * ( 1 - NU );
        double c3 = c1 * NU;
        double c4 = c1 * 0.5 * ( 1 - 2 * NU );
        C(0,0) = c2;  C(0,1) = c3;  C(0,2) = c3;  C(0,3) = 0.0; C(0,4) = 0.0; C(0,5) = 0.0;
        C(1,0) = c3;  C(1,1) = c2;  C(1,2) = c3;  C(1,3) = 0.0; C(1,4) = 0.0; C(1,5) = 0.0;
        C(2,0) = c3;  C(2,1) = c3;  C(2,2) = c2;  C(2,3) = 0.0; C(2,4) = 0.0; C(2,5) = 0.0;
        C(3,0) = 0.0; C(3,1) = 0.0; C(3,2) = 0.0; C(3,3) = c4;  C(3,4) = 0.0; C(3,5) = 0.0;
        C(4,0) = 0.0; C(4,1) = 0.0; C(4,2) = 0.0; C(4,3) = 0.0; C(4,4) = c4;  C(4,5) = 0.0;
        C(5,0) = 0.0; C(5,1) = 0.0; C(5,2) = 0.0; C(5,3) = 0.0; C(5,4) = 0.0; C(5,5) = c4;
    }
    else if (C.size1() == 3)    // plane strain
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
}

double IsotropicDamageIMPLEX::SofteningLaw( const double alpha) const
{
     double H = -(mE_0/mE_f)*exp(-(alpha-mE_0)/mE_f);

     return H;
}

double IsotropicDamageIMPLEX::DamageFunction(const double kappa) const
{
    return 1.0 - mE_0/kappa * std::exp(-(kappa-mE_0)/mE_f);
}

double IsotropicDamageIMPLEX::DamageFunctionDerivative(const double kappa) const
{
    return mE_0/kappa * (1.0/kappa + 1.0/mE_f) * std::exp(-(kappa-mE_0)/mE_f);
}

double IsotropicDamageIMPLEX::ComputeKappa(const double d, const double TOL, const int max_iters) const
{
    double kappa = mAlpha_old; // starting value

    int it = 0;
    do
    {
        double f = d - this->DamageFunction(kappa);

        if (std::abs(f) < TOL)
            return kappa;

        double df = this->DamageFunctionDerivative(kappa);
        kappa += f/df;

        ++it;

    } while (it < max_iters);

    KRATOS_ERROR << "The Newton-Raphson iteration does not converge within " << max_iters << " iterations";
}

int IsotropicDamageIMPLEX::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH( "" );
}

} // Namespace Kratos
