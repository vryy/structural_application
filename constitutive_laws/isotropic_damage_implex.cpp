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
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 11 Feb 2025 $
*   Revision:            $Revision: 1.13 $
*
* ***********************************************************/

// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#ifndef SD_APP_FORWARD_COMPATIBILITY
#include "includes/c2c_variables.h" // contain KAPPA
#endif
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "constitutive_laws/isotropic_damage_implex.h"
#include "structural_application_variables.h"

// #define DEBUG_DAMAGE

namespace Kratos
{

IsotropicDamageIMPLEX::IsotropicDamageIMPLEX()
    : ConstitutiveLaw()
{
}

IsotropicDamageIMPLEX::~IsotropicDamageIMPLEX()
{
}

bool IsotropicDamageIMPLEX::Has( const Variable<double>& rThisVariable )
{
    if( rThisVariable == YOUNG_MODULUS
     || rThisVariable == POISSON_RATIO
     || rThisVariable == DAMAGE_E0
     || rThisVariable == DAMAGE_EF
     || rThisVariable == DAMAGE )
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
    if( rThisVariable == DAMAGE_E0 )
        rValue = me0;
    if( rThisVariable == DAMAGE_EF )
        rValue = mef;
    if( rThisVariable == DAMAGE )
        rValue = mCurrentDamage;
    if( rThisVariable == KAPPA )
        rValue = mKappa;
    if( rThisVariable == EQUIVALENT_STRAIN )
    {
        SD_MathUtils<double>::Fourth_Order_Tensor De;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(De);
        SD_MathUtils<double>::CalculateElasticTensor(De, mE, mNU);

        Matrix effective_stress(3, 3);
        effective_stress.clear();
        SD_MathUtils<double>::ContractFourthOrderTensor(1.0, De, m_strain_n1, effective_stress);

        rValue = std::sqrt(SD_MathUtils<double>::mat_inner_prod(effective_stress, m_strain_n1) / mE) - mInitialEps;
    }

    return rValue;
}

Vector& IsotropicDamageIMPLEX::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if( rThisVariable == STRESSES )
        SD_MathUtils<double>::StressTensorToVector(m_stress_n1, rValue);

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
    if( rThisVariable == DAMAGE_E0 )
        me0 = rValue;
    if( rThisVariable == DAMAGE_EF )
        mef = rValue;
    if( rThisVariable == DAMAGE )
    {
        mCurrentDamage = rValue;
        mKappa = this->ComputeKappa(mCurrentDamage);
        // if (mKappa < me0) mKappa = me0; // this is not necessary, since the kappa will be > e0 when damage is in range [0, 1]
        mKappa_old = mKappa;
        // KRATOS_WATCH(rValue)
        // KRATOS_WATCH(mKappa)
        // if (mKappa == 0.0)
        //     KRATOS_ERROR << "Kappa should not be zero";
        if (mCurrentDamage > 1.0 || mCurrentDamage < 0.0)
        {
            KRATOS_ERROR << "The damage " << mCurrentDamage << " does not fall in range [0 1]";
        }
        // std::cout << "Damage is set to " << mCurrentDamage << std::endl;
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
        SD_MathUtils<double>::StressVectorToTensor(rValue, m_stress_n1);
}

void IsotropicDamageIMPLEX::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
}

void IsotropicDamageIMPLEX::InitializeMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    m_stress_n  = ZeroMatrix( 3, 3 );
    m_stress_n1 = ZeroMatrix( 3, 3 );
    m_stress_alg = ZeroMatrix( 3, 3 );
    m_strain_n1 = ZeroMatrix( 3, 3 );
    mE   = props[YOUNG_MODULUS];
    mNU  = props[POISSON_RATIO];
    me0 = props[DAMAGE_E0];
    mef = props[DAMAGE_EF];
    mKappa     = me0;
    mKappa_old = mKappa;
    mKappa_old_old = mKappa;
    mq         = mKappa;
    mq_old     = mKappa;
    mCurrentDamage = 0.0;
    mCurrentDamageAlg = 0.0;
    mInitialEps    = 0.0;
    mDamageFlag    = 0;
}

void IsotropicDamageIMPLEX::ResetMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    this->ResetState();
}

void IsotropicDamageIMPLEX::ResetState()
{
    mKappa = me0;
    mKappa_old = mKappa;
    mq_old = mq;
    m_stress_n1.clear();
    m_stress_alg.clear();
    mCurrentDamage = 0.0;
    mCurrentDamageAlg = 0.0;
    mDamageFlag = 0;
}

void IsotropicDamageIMPLEX::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void IsotropicDamageIMPLEX::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void IsotropicDamageIMPLEX::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void IsotropicDamageIMPLEX::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mKappa_old_old = mKappa_old;
    mKappa_old = mKappa;
    mq_old = mq;
    noalias(m_stress_n) = m_stress_n1;
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

void  IsotropicDamageIMPLEX::CalculateMaterialResponse( const Vector& StrainVector,
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
            SD_MathUtils<double>::StressTensorToVector(m_stress_n1, StressVector);
            return;
        }
    }

    double TOL = 1.0e-8;
    if (props.Has(LOCAL_ERROR_TOLERANCE))
        TOL = props[LOCAL_ERROR_TOLERANCE];

    if (CalculateStresses)
    {
        this->StressIntegration(StrainVector, TOL, CurrentProcessInfo, props);

        SD_MathUtils<double>::StressTensorToVector(m_stress_alg, StressVector);

        #ifdef DEBUG_DAMAGE
        // if (mElemId == 1 && mGaussId == 0)
        {
            KRATOS_WATCH(StressVector)
        }
        #endif
    }

    if (CalculateTangent)
    {
        this->ComputeTangent(AlgorithmicTangent, CurrentProcessInfo, props);

        #ifdef DEBUG_DAMAGE
        // if (mElemId == 1 && mGaussId == 0)
        {
            KRATOS_WATCH(AlgorithmicTangent)
        }
        #endif
    }

    #ifdef DEBUG_DAMAGE
    // if (mElemId == 1 && mGaussId == 0)
    {
        KRATOS_WATCH("------------")
    }
    #endif

    // KRATOS_WATCH(StressVector)
    // KRATOS_WATCH(AlgorithmicTangent)
}

void IsotropicDamageIMPLEX::StressIntegration(const Vector& StrainVector, const double TOL,
            const ProcessInfo& CurrentProcessInfo, const Properties& rProperties)
{
    /***************** Implicit integration *****************/

    // equivalent strain
    SD_MathUtils<double>::StrainVectorToTensor(StrainVector, m_strain_n1);

    SD_MathUtils<double>::Fourth_Order_Tensor De;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(De);
    SD_MathUtils<double>::CalculateElasticTensor(De, mE, mNU);

    Matrix effective_stress(3, 3);
    effective_stress.clear();
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0, De, m_strain_n1, effective_stress);

    double eps = std::sqrt(SD_MathUtils<double>::mat_inner_prod(effective_stress, m_strain_n1) / mE) - mInitialEps;

    // loading check
    const double H = this->SofteningLaw(mKappa_old);
    mKappa = mKappa_old;
    mq = mq_old;
    if (mDamageFlag != -1) // -1 is the special value. It will stop the damage update
    {
        // if( eps > mKappa + TOL )
        if( eps > mKappa )
        {
            mKappa = eps;

            mq = mq_old + H*(mKappa - mKappa_old);

            mDamageFlag = 1;
        }
        else
        {
            mDamageFlag = 0;
        }

        // damage parameter
        mCurrentDamage = (1.0 - mq/mKappa);
    }

    // stress update
    noalias(m_stress_n1) = (1.0 - mCurrentDamage) * effective_stress;

    /***************** Explicit integration *****************/

    // compute extrapolated damage variable
    double dkappa = mKappa_old - mKappa_old_old;
    // KRATOS_WATCH(dkappa)
    double kappa_alg = mKappa_old + dkappa;

    double q_alg = mq_old + H*dkappa;

    mCurrentDamageAlg = 1.0 - ( q_alg / kappa_alg );

    // algorithmic stress update
    noalias(m_stress_alg) = (1.0 - mCurrentDamageAlg) * effective_stress;
}

void IsotropicDamageIMPLEX::ComputeTangent(Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo, const Properties& rProperties) const
{
    SD_MathUtils<double>::Fourth_Order_Tensor D;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(D);
    SD_MathUtils<double>::CalculateElasticTensor(D, mE, mNU);

    SD_MathUtils<double>::ScaleFourthOrderTensor(D, 1.0 - mCurrentDamageAlg);

    SD_MathUtils<double>::TensorToMatrix(D, AlgorithmicTangent);
}

double IsotropicDamageIMPLEX::DamageFunction(const double kappa) const
{
    return 1.0 - me0/kappa * std::exp(-(kappa-me0)/mef); // TODO check this damage function. It is too sensitive.
}

double IsotropicDamageIMPLEX::DamageFunctionDerivative(const double kappa) const
{
    return me0/kappa * (1.0/kappa + 1.0/mef) * std::exp(-(kappa-me0)/mef);
}

double IsotropicDamageIMPLEX::SofteningLaw( const double kappa) const
{
    return -(me0/mef) * std::exp(-(kappa-me0)/mef);
}

double IsotropicDamageIMPLEX::ComputeKappa(const double d, const double TOL, const int max_iters) const
{
    double kappa = mKappa_old; // starting value

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

//**********************************************************************
int IsotropicDamageIMPLEX::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH( "" );
}

} // Namespace Kratos

#undef DEBUG_DAMAGE
