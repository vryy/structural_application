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
#include "constitutive_laws/dummy_constitutive_law.h"
#include "structural_application_variables.h"

namespace Kratos
{


/**
 * TO BE TESTED!!!
 */
DummyConstitutiveLaw::DummyConstitutiveLaw()
    : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
DummyConstitutiveLaw::~DummyConstitutiveLaw()
{
}

bool DummyConstitutiveLaw::Has( const Variable<int>& rThisVariable )
{
    return false;
}

bool DummyConstitutiveLaw::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool DummyConstitutiveLaw::Has( const Variable<Vector>& rThisVariable )
{
    return false;
}

bool DummyConstitutiveLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

int& DummyConstitutiveLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
        rValue = mElemId;
    if(rThisVariable == INTEGRATION_POINT_INDEX)
        rValue = mGaussId;
    return rValue;
}

double& DummyConstitutiveLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    rValue = 0.0;
    return rValue;
}

Vector& DummyConstitutiveLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    rValue.clear();
    return rValue;
}

Matrix& DummyConstitutiveLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    rValue.clear();
    return rValue;
}

void DummyConstitutiveLaw::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
        mElemId = rValue;
    if(rThisVariable == INTEGRATION_POINT_INDEX)
        mGaussId = rValue;
}

void DummyConstitutiveLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
}

void DummyConstitutiveLaw::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                                     const array_1d<double, 3>& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
}

void DummyConstitutiveLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
}

void DummyConstitutiveLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
}


void DummyConstitutiveLaw::InitializeMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
}

void DummyConstitutiveLaw::ResetMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
}

void DummyConstitutiveLaw::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    // std::stringstream ss;
    // ss << "Error calling DummyConstitutiveLaw::" << __FUNCTION__
    //    << ", Properties: " << props.Id()
    //    << ", Element: " << mElemId
    //    << ", Point: " << mGaussId;
    // KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
}

void DummyConstitutiveLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    std::stringstream ss;
    ss << "Error calling DummyConstitutiveLaw::" << __FUNCTION__
       << ", Properties: " << rMaterialProperties.Id()
       << ", Element: " << mElemId
       << ", Point: " << mGaussId;
    KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
}

void DummyConstitutiveLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    std::stringstream ss;
    ss << "Error calling DummyConstitutiveLaw::" << __FUNCTION__
       << ", Properties: " << rMaterialProperties.Id()
       << ", Element: " << mElemId
       << ", Point: " << mGaussId;
    KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
}

void DummyConstitutiveLaw::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    // std::stringstream ss;
    // ss << "Error calling DummyConstitutiveLaw::" << __FUNCTION__
    //    << ", Properties: " << props.Id()
    //    << ", Element: " << mElemId
    //    << ", Point: " << mGaussId;
    // KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
}

void  DummyConstitutiveLaw::CalculateMaterialResponse( const Vector& StrainVector,
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
    std::stringstream ss;
    ss << "Error calling DummyConstitutiveLaw::" << __FUNCTION__
       << ", Properties: " << props.Id()
       << ", Element: " << mElemId
       << ", Point: " << mGaussId;
    KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
}

//**********************************************************************
int DummyConstitutiveLaw::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH( "" );
}
} // Namespace Kratos
