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
#include "constitutive_laws/isotropic_damage_model.h"
#include "structural_application_variables.h"

namespace Kratos
{

IsotropicDamageModel::IsotropicDamageModel()
    : ConstitutiveLaw()
{
}

IsotropicDamageModel::~IsotropicDamageModel()
{
}

bool IsotropicDamageModel::Has( const Variable<double>& rThisVariable )
{
    if( rThisVariable == YOUNG_MODULUS
     || rThisVariable == POISSON_RATIO
     || rThisVariable == DAMAGE_E0
     || rThisVariable == DAMAGE_EF
     || rThisVariable == DAMAGE )
        return true;
    return false;
}

bool IsotropicDamageModel::Has( const Variable<Vector>& rThisVariable )
{
    if( rThisVariable == STRESSES )
        return true;
    return false;
}

bool IsotropicDamageModel::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

double& IsotropicDamageModel::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if( rThisVariable == YOUNG_MODULUS )
        rValue = mE;
    if( rThisVariable == POISSON_RATIO )
        rValue = mNU;
    if( rThisVariable == DAMAGE_E0 )
        rValue = mE_0;
    if( rThisVariable == DAMAGE_EF )
        rValue = mE_f;
    if( rThisVariable == DAMAGE )
        rValue = mCurrentDamage;
    if( rThisVariable == INITIAL_DAMAGE )
        rValue = mInitialDamage;
    if( rThisVariable == KAPPA )
        rValue = mKappa;
    if( rThisVariable == EQUIVALENT_STRAIN )
    {
        SD_MathUtils<double>::Fourth_Order_Tensor De;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(De);
        SD_MathUtils<double>::CalculateElasticTensor(De, mE, mNU);

        Matrix stress_no_damage(3, 3);
        stress_no_damage.clear();
        SD_MathUtils<double>::ContractFourthOrderTensor(1.0, De, m_strain_n1, stress_no_damage);

        rValue = std::sqrt(SD_MathUtils<double>::mat_inner_prod(stress_no_damage, m_strain_n1) / mE) - mInitialEps;
    }

    return rValue;
}

Vector& IsotropicDamageModel::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if( rThisVariable == STRESSES )
        SD_MathUtils<double>::StressTensorToVector(m_stress_n1, rValue);

    return rValue;
}

Matrix& IsotropicDamageModel::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return rValue;
}

void IsotropicDamageModel::SetValue( const Variable<bool>& rThisVariable, const bool& rValue,
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

void IsotropicDamageModel::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == PARENT_ELEMENT_ID )
        mElemId = rValue;
    if( rThisVariable == INTEGRATION_POINT_INDEX )
        mGaussId = rValue;
}

void IsotropicDamageModel::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if( rThisVariable == POISSON_RATIO )
        mNU = rValue;
    if( rThisVariable == DAMAGE_E0 )
        mE_0 = rValue;
    if( rThisVariable == DAMAGE_EF )
        mE_f = rValue;
    if( rThisVariable == INITIAL_DAMAGE )
    {
        mInitialDamage = rValue;
        mCurrentDamage = mInitialDamage;
        // std::cout << "Initial damage is set to " << mInitialDamage << std::endl;
    }
    if( rThisVariable == DAMAGE )
    {
        mCurrentDamage = rValue;
        mKappa = this->ComputeKappa(mCurrentDamage);
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

void IsotropicDamageModel::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                                     const array_1d<double, 3>& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
}

void IsotropicDamageModel::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == STRESSES )
        SD_MathUtils<double>::StressVectorToTensor(rValue, m_stress_n1);
}

void IsotropicDamageModel::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
}

void IsotropicDamageModel::InitializeMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    m_stress_n  = ZeroMatrix( 3, 3 );
    m_stress_n1 = ZeroMatrix( 3, 3 );
    m_strain_n1 = ZeroMatrix( 3, 3 );
    mE   = props[YOUNG_MODULUS];
    mNU  = props[POISSON_RATIO];
    mE_0 = props[DAMAGE_E0];
    mE_f = props[DAMAGE_EF];
    mKappa     = mE_0;
    mKappa_old = mKappa;
    mCurrentDamage = 0.0;
    mInitialDamage = 0.0;
    mInitialEps    = 0.0;
    mDamageFlag    = 0;
}

void IsotropicDamageModel::ResetMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    this->ResetState();
}

void IsotropicDamageModel::ResetState()
{
    mKappa = mE_0;
    mKappa_old = mKappa;
    m_stress_n1.clear();
    mCurrentDamage = mInitialDamage;
    mDamageFlag = 0;
}

void IsotropicDamageModel::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void IsotropicDamageModel::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    // restore mKappa
    mKappa = mKappa_old;
}

void IsotropicDamageModel::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
}

void IsotropicDamageModel::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mKappa_old = mKappa;
    noalias(m_stress_n) = m_stress_n1;
}

void IsotropicDamageModel::CalculateMaterialResponseCauchy (Parameters& parameters)
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

void  IsotropicDamageModel::CalculateMaterialResponse( const Vector& StrainVector,
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

        SD_MathUtils<double>::StressTensorToVector(m_stress_n1, StressVector);

        {
        }
    }

    if (CalculateTangent)
    {
        this->ComputeTangent(AlgorithmicTangent, CurrentProcessInfo, props);
    }


        // KRATOS_WATCH(mE_0)
        // KRATOS_WATCH(mE_f)
        // KRATOS_WATCH(mKappa)
        // KRATOS_WATCH(d)
        // KRATOS_WATCH(this->DamageFunctionDerivative(mKappa))
}

void IsotropicDamageModel::StressIntegration(const Vector& StrainVector, const double TOL,
            const ProcessInfo& CurrentProcessInfo, const Properties& rProperties)
{
    // equivalent strain
    SD_MathUtils<double>::StrainVectorToTensor(StrainVector, m_strain_n1);

    SD_MathUtils<double>::Fourth_Order_Tensor De;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(De);
    SD_MathUtils<double>::CalculateElasticTensor(De, mE, mNU);

    Matrix stress_no_damage(3, 3);
    stress_no_damage.clear();
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0, De, m_strain_n1, stress_no_damage);

    double eps = std::sqrt(SD_MathUtils<double>::mat_inner_prod(stress_no_damage, m_strain_n1) / mE) - mInitialEps;

    // loading check
    if (mDamageFlag != -1) // -1 is the special value. It will stop the damage update
    {
        if( eps > mKappa + TOL )
        {
            mKappa = eps;
            mDamageFlag = 1;
        }
        else
        {
            mDamageFlag = 0;
        }

        // damage parameter
        // here we do not allow the damage to exceed 1 by multiplying it with the rest of the initial damage
        mCurrentDamage = mInitialDamage + (1.0 - mInitialDamage) * this->DamageFunction(mKappa);
    }

    // KRATOS_WATCH(StressVector)
    // KRATOS_WATCH(AlgorithmicTangent)
    // stress update
    noalias(m_stress_n1) = (1.0 - mCurrentDamage) * stress_no_damage;

}

void IsotropicDamageModel::ComputeTangent(Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo, const Properties& rProperties) const
{
    SD_MathUtils<double>::Fourth_Order_Tensor D;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(D);
    SD_MathUtils<double>::CalculateElasticTensor(D, mE, mNU);

    Matrix stress_no_damage(3, 3);
    stress_no_damage.clear();
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0, D, m_strain_n1, stress_no_damage);

    SD_MathUtils<double>::ScaleFourthOrderTensor(D, 1.0 - mCurrentDamage);
    if (mDamageFlag > 0)
    {
        double eps = std::sqrt(SD_MathUtils<double>::mat_inner_prod(stress_no_damage, m_strain_n1) / mE);

        Matrix dd_de(3, 3); // derivatives of mCurrentDamage w.r.t strain tensor
        noalias(dd_de) = (1.0 - mInitialDamage) * this->DamageFunctionDerivative(mKappa) * stress_no_damage / (mE*eps);
        SD_MathUtils<double>::OuterProductFourthOrderTensor(-1.0, stress_no_damage, dd_de, D);
    }

    SD_MathUtils<double>::TensorToMatrix(D, AlgorithmicTangent);
}

double IsotropicDamageModel::DamageFunction(const double kappa) const
{
    return 1.0 - mE_0/kappa * std::exp(-(kappa-mE_0)/mE_f);
}

double IsotropicDamageModel::DamageFunctionDerivative(const double kappa) const
{
    return mE_0/kappa * (1.0/kappa + 1.0/mE_f) * std::exp(-(kappa-mE_0)/mE_f);
}

double IsotropicDamageModel::ComputeKappa(const double d, const double TOL, const int max_iters) const
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
int IsotropicDamageModel::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH( "" );
}

} // Namespace Kratos
