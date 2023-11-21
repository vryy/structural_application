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
*   Date:                $Date: 17/11/2023 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_STRUCTURAL_APPLICATION_ST_VENANT_KIRCHHOFF_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APPLICATION_ST_VENANT_KIRCHHOFF_H_INCLUDED

// System includes

// External includes

// Project includes
#include "constitutive_laws/plane_strain.h"
#include "constitutive_laws/isotropic_3d.h"


namespace Kratos
{

template<int TDim>
class StVenantKirchhoff;

template<>
class StVenantKirchhoff<2> : public PlaneStrain
{
    typedef PlaneStrain BaseType;

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone(new StVenantKirchhoff<2>());
        return p_clone;
    }

    ConstitutiveLaw::StrainMeasure GetStrainMeasure() override
    {
        return StrainMeasure_GreenLagrange;
    }

    ConstitutiveLaw::StressMeasure GetStressMeasure() override
    {
        return StressMeasure_PK2;
    }

    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue ) override
    {
        if( rThisVariable == PK2_STRESS_TENSOR )
        {
            if (rValue.size1() != 3 || rValue.size2() != 3)
                rValue.resize(3, 3, false);
            SD_MathUtils<double>::StressVectorToTensor(mCurrentStress, rValue);
        }

        return BaseType::GetValue(rThisVariable, rValue);
    }

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(Parameters& rValues) final
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }

    std::string Info() const override
    {
        return "StVenantKirchhoff_PlaneStrain";
    }
};

template<>
class StVenantKirchhoff<3> : public Isotropic3D
{
    typedef Isotropic3D BaseType;

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone(new StVenantKirchhoff<3>());
        return p_clone;
    }

    ConstitutiveLaw::StrainMeasure GetStrainMeasure() override
    {
        return StrainMeasure_GreenLagrange;
    }

    ConstitutiveLaw::StressMeasure GetStressMeasure() override
    {
        return StressMeasure_PK2;
    }

    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue ) override
    {
        if( rThisVariable == PK2_STRESS_TENSOR )
        {
            if (rValue.size1() != 3 || rValue.size2() != 3)
                rValue.resize(3, 3, false);
            SD_MathUtils<double>::StressVectorToTensor(mCurrentStress, rValue);
        }

        return BaseType::GetValue(rThisVariable, rValue);
    }

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(Parameters& rValues) final
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }

    std::string Info() const override
    {
        return "StVenantKirchhoff_3D";
    }
};

}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_APPLICATION_ST_VENANT_KIRCHHOFF_H_INCLUDED  defined
