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
*   Date:                $Date: 21 Sep 2023 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/


// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws/isotropic_3d_dc.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"

namespace Kratos
{
/**
 * TO BE TESTED!!!
 */
Isotropic3DDC::Isotropic3DDC()
    : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
Isotropic3DDC::~Isotropic3DDC()
{
}

bool Isotropic3DDC::Has(const Variable<int>& rThisVariable)
{
    return false;
}

bool Isotropic3DDC::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool Isotropic3DDC::Has( const Variable<Vector>& rThisVariable )
{
    if(rThisVariable == STRESSES)
        return true;
    if(rThisVariable == CURRENT_STRAIN_VECTOR)
        return true;

    return false;
}

bool Isotropic3DDC::Has( const Variable<Matrix>& rThisVariable )
{
    if(rThisVariable == ALGORITHMIC_TANGENT || rThisVariable == THREED_ALGORITHMIC_TANGENT)
        return true;
    return false;
}

int& Isotropic3DDC::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    if (rThisVariable == IS_SHAPE_FUNCTION_REQUIRED)
        rValue = 0;

    return rValue;
}

double& Isotropic3DDC::GetValue( const Variable<double>& rThisVariable, double& rValue )
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
        rValue = -(m_stress_n(0, 0) + m_stress_n(1, 1) + m_stress_n(2, 2)) / 3;
        return rValue;
    }

    if (rThisVariable == PRESSURE_Q || rThisVariable == VON_MISES_STRESS)
    {
        double p = (m_stress_n(0, 0) + m_stress_n(1, 1) + m_stress_n(2, 2)) / 3;
        double sxx = m_stress_n(0, 0) - p;
        double syy = m_stress_n(1, 1) - p;
        double szz = m_stress_n(2, 2) - p;
        double sxy = m_stress_n(0, 1);
        double syz = m_stress_n(1, 2);
        double sxz = m_stress_n(0, 2);

        rValue = sqrt( 3.0 * ( 0.5*(sxx*sxx + syy*syy + szz*szz) + sxy*sxy + syz*syz + sxz*sxz ) );
        return rValue;
    }

    return rValue;
}

Vector& Isotropic3DDC::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == STRESSES || rThisVariable == STRESSES_OLD )
    {
        SD_MathUtils<double>::StressTensorToVector(m_stress_n, rValue);
    }
    else if ( rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS )
    {
        rValue = mPreStress;
    }
    else if ( rThisVariable == THREED_STRESSES )
    {
        if(rValue.size() != 6)
            rValue.resize(6, false);
        rValue(0) = m_stress_n(0, 0);
        rValue(1) = m_stress_n(1, 1);
        rValue(2) = m_stress_n(2, 2);
        rValue(3) = m_stress_n(0, 1);
        rValue(4) = m_stress_n(1, 2);
        rValue(5) = m_stress_n(0, 2);
    }
    else if ( rThisVariable == THREED_STRAIN )
    {
        // REF: http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_strain.cfm
        if(rValue.size() != 6)
            rValue.resize(6, false);
        rValue(0) = m_strain_n(0, 0);
        rValue(1) = m_strain_n(1, 1);
        rValue(2) = m_strain_n(2, 2);
        rValue(3) = 2*m_strain_n(0, 1);
        rValue(4) = 2*m_strain_n(1, 2);
        rValue(5) = 2*m_strain_n(0, 2);
    }
    else if ( rThisVariable == ELASTIC_STRAIN_VECTOR || rThisVariable == STRAIN )
    {
        SD_MathUtils<double>::StrainTensorToVector(m_strain_n, rValue);
    }

    return rValue;
}

Matrix& Isotropic3DDC::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    if (rThisVariable == CAUCHY_STRESS_TENSOR)
    {
        noalias(rValue) = m_stress_n1;
    }
    else if(rThisVariable == ALGORITHMIC_TANGENT || rThisVariable == ELASTIC_TANGENT)
    {
        SD_MathUtils<double>::Fourth_Order_Tensor De;
        SD_MathUtils<double>::CalculateElasticTensor(De, mE, mNU);
        SD_MathUtils<double>::TensorToMatrix(De, rValue);
    }
    else if(rThisVariable == THREED_ALGORITHMIC_TANGENT)
    {
        if(rValue.size1() != 6 || rValue.size2() != 6)
            rValue.resize(6, 6, false);
        SD_MathUtils<double>::Fourth_Order_Tensor De;
        SD_MathUtils<double>::CalculateElasticTensor(De, mE, mNU);
        SD_MathUtils<double>::TensorToMatrix(De, rValue);
    }
    else if(rThisVariable == ELASTIC_STRAIN_TENSOR)
    {
        noalias(rValue) = m_strain_n;
    }

    return( rValue );
}

void Isotropic3DDC::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                              const ProcessInfo& rCurrentProcessInfo )
{
}

void Isotropic3DDC::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                              const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
    if ( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;
}

void Isotropic3DDC::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                              const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS )
    {
        noalias(mPreStress) = rValue;
    }
    else if ( rThisVariable == STRESSES || rThisVariable == INITIAL_STRESS )
    {
        SD_MathUtils<double>::StressVectorToTensor(rValue, m_stress_n);
        SD_MathUtils<double>::StressVectorToTensor(rValue, m_stress_n1);
    }
    else if ( rThisVariable == THREED_STRESSES )
    {
        SD_MathUtils<double>::StressVectorToTensor(rValue, m_stress_n);
        SD_MathUtils<double>::StressVectorToTensor(rValue, m_stress_n1);
    }
    else if ( rThisVariable == ELASTIC_STRAIN_VECTOR || rThisVariable == CURRENT_STRAIN_VECTOR )
    {
        SD_MathUtils<double>::StrainVectorToTensor(rValue, m_strain_n1);
    }
}

void Isotropic3DDC::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                              const ProcessInfo& rCurrentProcessInfo )
{
}

void Isotropic3DDC::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    if (rValues.IsSetStressVector())
    {
        Vector& StressVector = rValues.GetStressVector();
        SD_MathUtils<double>::StressTensorToVector(m_stress_n1, StressVector);
    }

    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        SD_MathUtils<double>::Fourth_Order_Tensor De;
        SD_MathUtils<double>::CalculateElasticTensor(De, mE, mNU);
        SD_MathUtils<double>::TensorToMatrix(De, AlgorithmicTangent);
    }
}

void Isotropic3DDC::CalculateMaterialResponsePK2 (Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void Isotropic3DDC::CalculateMaterialResponse( const Vector& StrainVector,
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
void Isotropic3DDC::InitializeMaterial( const Properties& props,
                                        const GeometryType& geom,
                                        const Vector& ShapeFunctionsValues )
{
    m_strain_n = ZeroMatrix(3, 3);
    m_strain_n1 = ZeroMatrix(3, 3);
    m_stress_n = ZeroMatrix(3, 3);
    m_stress_n1 = ZeroMatrix(3, 3);
    mPreStress = ZeroVector( 3 );
    mPrestressFactor = 1.0;
    mE  = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mDE = props[DENSITY];
}

void Isotropic3DDC::ResetMaterial( const Properties& props,
                                   const GeometryType& geom,
                                   const Vector& ShapeFunctionsValues )
{
    noalias(m_strain_n) = ZeroMatrix(3, 3);
    noalias(m_strain_n1) = ZeroMatrix(3, 3);

    Matrix prestress(3, 3);
    SD_MathUtils<double>::StressVectorToTensor(-mPreStress, prestress);

    noalias(m_stress_n) = mPrestressFactor*prestress;
    noalias(m_stress_n1) = m_stress_n;
}

void Isotropic3DDC::InitializeNonLinearIteration( const Properties& rMaterialProperties,
                                                  const GeometryType& rElementGeometry,
                                                  const Vector& rShapeFunctionsValues,
                                                  const ProcessInfo& rCurrentProcessInfo )
{
}

void Isotropic3DDC::FinalizeNonLinearIteration( const Properties& rMaterialProperties,
                                                const GeometryType& rElementGeometry,
                                                const Vector& rShapeFunctionsValues,
                                                const ProcessInfo& rCurrentProcessInfo )
{
    CalculateStress(m_strain_n1, m_stress_n1);
}

void Isotropic3DDC::FinalizeSolutionStep ( const Properties& props,
                                           const GeometryType& geom, //this is just to give the array of nodes
                                           const Vector& ShapeFunctionsValues ,
                                           const ProcessInfo& CurrentProcessInfo )
{
    noalias( m_stress_n ) = m_stress_n1;
    noalias( m_strain_n ) = m_strain_n1;
}

void Isotropic3DDC::CalculateStress( const Matrix& strain_tensor, Matrix& stress_tensor ) const
{
    noalias(stress_tensor) = m_stress_n;
    Matrix incremental_strain = strain_tensor - m_strain_n;

    SD_MathUtils<double>::Fourth_Order_Tensor De;
    SD_MathUtils<double>::CalculateElasticTensor(De, mE, mNU);
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0, De, incremental_strain, stress_tensor);

    // KRATOS_WATCH(m_strain_n)
    // KRATOS_WATCH(m_stress_n)
    // KRATOS_WATCH(strain_tensor)
    // KRATOS_WATCH(stress_tensor)
}

//**********************************************************************

int Isotropic3DDC::Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo) const
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
