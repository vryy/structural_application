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
*   Date:                $Date: 5 Mar 2020 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws/values_container_constitutive_law.h"
#include "structural_application_variables.h"

namespace Kratos
{

ValuesContainerConstitutiveLaw::ValuesContainerConstitutiveLaw()
: ConstitutiveLaw(), mpConstitutiveLaw(new ConstitutiveLaw())
{
}

ValuesContainerConstitutiveLaw::ValuesContainerConstitutiveLaw(ConstitutiveLaw::Pointer pOther)
: mpConstitutiveLaw(pOther)
{
}

/**
 * TO BE TESTED!!!
 */
ValuesContainerConstitutiveLaw::~ValuesContainerConstitutiveLaw()
{
}

bool ValuesContainerConstitutiveLaw::Has( const Variable<int>& rThisVariable )
{
    if (mValuesPosition.find(rThisVariable.Key()) != mValuesPosition.end())
        return true;

    return mpConstitutiveLaw->Has( rThisVariable );
}

bool ValuesContainerConstitutiveLaw::Has( const Variable<double>& rThisVariable )
{
    if (mValuesPosition.find(rThisVariable.Key()) != mValuesPosition.end())
        return true;

    return mpConstitutiveLaw->Has( rThisVariable );
}

bool ValuesContainerConstitutiveLaw::Has( const Variable<array_1d<double, 3> >& rThisVariable )
{
    if (mValuesPosition.find(rThisVariable.Key()) != mValuesPosition.end())
        return true;

    return mpConstitutiveLaw->Has( rThisVariable );
}

bool ValuesContainerConstitutiveLaw::Has( const Variable<Vector>& rThisVariable )
{
    if (mValuesPosition.find(rThisVariable.Key()) != mValuesPosition.end())
        return true;

    return mpConstitutiveLaw->Has( rThisVariable );
}

bool ValuesContainerConstitutiveLaw::Has( const Variable<Matrix>& rThisVariable )
{
    if (mValuesPosition.find(rThisVariable.Key()) != mValuesPosition.end())
        return true;

    return mpConstitutiveLaw->Has( rThisVariable );
}

int& ValuesContainerConstitutiveLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    ValuesKeyMapType::iterator it = mValuesPosition.find(rThisVariable.Key());
    if (it != mValuesPosition.end())
    {
        rValue = mIntValues[it->second];
        return rValue;
    }

    return mpConstitutiveLaw->GetValue( rThisVariable, rValue );
}

double& ValuesContainerConstitutiveLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    ValuesKeyMapType::iterator it = mValuesPosition.find(rThisVariable.Key());
    if (it != mValuesPosition.end())
    {
        rValue = mDoubleValues[it->second];
        return rValue;
    }

    return mpConstitutiveLaw->GetValue( rThisVariable, rValue );
}

array_1d<double, 3>& ValuesContainerConstitutiveLaw::GetValue( const Variable<array_1d<double, 3> >& rThisVariable, array_1d<double, 3>& rValue )
{
    return mpConstitutiveLaw->GetValue( rThisVariable, rValue );
}

Vector& ValuesContainerConstitutiveLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return mpConstitutiveLaw->GetValue( rThisVariable, rValue );
}

Matrix& ValuesContainerConstitutiveLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return mpConstitutiveLaw->GetValue( rThisVariable, rValue );
}

void ValuesContainerConstitutiveLaw::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue( rThisVariable, rValue, rCurrentProcessInfo );
}

void ValuesContainerConstitutiveLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    ValuesKeyMapType::iterator it = mValuesPosition.find(rThisVariable.Key());
    if (it == mValuesPosition.end())
    {
        mDoubleValues.push_back(rValue);
        mValuesPosition[rThisVariable.Key()] = mDoubleValues.size()-1;
    }
    else
    {
        mDoubleValues[it->second] = rValue;
    }

    mpConstitutiveLaw->SetValue( rThisVariable, rValue, rCurrentProcessInfo );
}

void ValuesContainerConstitutiveLaw::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue( rThisVariable, rValue, rCurrentProcessInfo );
}

void ValuesContainerConstitutiveLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue( rThisVariable, rValue, rCurrentProcessInfo );
}

void ValuesContainerConstitutiveLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue( rThisVariable, rValue, rCurrentProcessInfo );
}

void ValuesContainerConstitutiveLaw::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->InitializeMaterial( props, geom, ShapeFunctionsValues );
}

void ValuesContainerConstitutiveLaw::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->ResetMaterial( props, geom, ShapeFunctionsValues );
}

void ValuesContainerConstitutiveLaw::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeSolutionStep( props, geom, ShapeFunctionsValues, CurrentProcessInfo );
}

void ValuesContainerConstitutiveLaw::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeNonLinearIteration( props, geom, ShapeFunctionsValues, CurrentProcessInfo );
}

void ValuesContainerConstitutiveLaw::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeNonLinearIteration( props, geom, ShapeFunctionsValues, CurrentProcessInfo );
}

void ValuesContainerConstitutiveLaw::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeSolutionStep( props, geom, ShapeFunctionsValues, CurrentProcessInfo );
}

void  ValuesContainerConstitutiveLaw::CalculateMaterialResponse( const Vector& StrainVector,
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
    mpConstitutiveLaw->CalculateMaterialResponse( StrainVector,
        DeformationGradient,
        StressVector,
        AlgorithmicTangent,
        CurrentProcessInfo,
        props,
        geom,
        ShapeFunctionsValues,
        CalculateStresses,
        CalculateTangent,
        SaveInternalVariables );
}

//**********************************************************************
int ValuesContainerConstitutiveLaw::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    return mpConstitutiveLaw->Check(props, geom, CurrentProcessInfo);
}

} // Namespace Kratos
