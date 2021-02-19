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
*   Last Modified by:    $Author: hurga $
*   Date:                $Date: 2008-02-15 10:37:19 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

// System includes


// External includes

// Project includes
#include "custom_elements/truss_element.h"
#include "structural_application_variables.h"

namespace Kratos
{

TrussElement::TrussElement(IndexType NewId,
                           GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
    //THIS IS THE DEFAULT CONSTRUCTOR
}

//************************************************************************************
//************************************************************************************
//THIS IS A TRUSS ELEMENT FOR FINITE STRAINS (ALSO CALLED CRISFIELD-ELEMENT)
// it consists of a 2Node Element with neo-Hook Material behaviour
// all operations are done inside the elements
// it does not contain works only with line2d and does not contain a material
//************************************************************************************
//************************************************************************************
TrussElement::TrussElement(IndexType NewId,
                           GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    //DEFINE HERE THE INTEGRATION METHOD (the number of integration points to be used)
    //see /kratos/source/integration_rules.cpp for further details
    if(GetGeometry().size() != 2)
    {
        std::cout<<"this element works only with a 2 node line"<<std::endl;
    }
}

Element::Pointer TrussElement::Create(IndexType NewId,
                                      NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TrussElement(NewId, GetGeometry().Create(ThisNodes),
                            pProperties));
}

//************************************************************************************
//************************************************************************************
//THIS IS THE DESTRUCTOR, as we use UBLAS classes, we don't have to care about garbage collecting
//************************************************************************************
//************************************************************************************
TrussElement::~TrussElement()
{
}
//************************************************************************************
//THIS IS THE INITIALIZATION OF THE ELEMENT (CALLED AT THE BEGIN OF EACH CALCULATION)
//************************************************************************************
void TrussElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

//************************************************************************************
//THIS is the main method here the integration in space (loop over the integration points) is done
//************************************************************************************
void TrussElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo,
                                bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    //resizing and calculate the LHS contribution if needed
    if (CalculateStiffnessMatrixFlag || CalculateResidualVectorFlag)
    {
        if(rLeftHandSideMatrix.size1() != 6 || rLeftHandSideMatrix.size2() != 6)
            rLeftHandSideMatrix.resize(6, 6, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(6, 6);

        double lx = GetGeometry()[1].X0() - GetGeometry()[0].X0();
        double ly = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
        double lz = GetGeometry()[1].Z0() - GetGeometry()[0].Z0();
        double l = sqrt(lx*lx + ly*ly + lz*lz);
        double E = GetProperties()[YOUNG_MODULUS];
        double A = GetProperties()[AREA];

        Matrix T(2, 6);
        noalias(T) = ZeroMatrix(2, 6);
        T(0, 0) = lx/l;
        T(0, 1) = ly/l;
        T(0, 2) = lz/l;
        T(1, 3) = lx/l;
        T(1, 4) = ly/l;
        T(1, 5) = lz/l;

        Matrix ke(2, 2);
        double k = E*A/l;
        ke(0, 0) = k;
        ke(0, 1) = -k;
        ke(1, 0) = -k;
        ke(1, 1) = k;

        noalias(rLeftHandSideMatrix) += prod(trans(T), Matrix(prod(ke, T)));
    }

    //resizing and calculate the RHS contribution if needed
    if (CalculateResidualVectorFlag)
    {
        if(rRightHandSideVector.size() != 6)
            rRightHandSideVector.resize(6, false);
        noalias(rRightHandSideVector) = ZeroVector(6);

        // get the current displacement
        Vector CurrentDisplacement(6);
        CurrentDisplacement(0) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X);
        CurrentDisplacement(1) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y);
        CurrentDisplacement(2) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z);
        CurrentDisplacement(3) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X);
        CurrentDisplacement(4) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y);
        CurrentDisplacement(5) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z);

        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, CurrentDisplacement);
    }

    // KRATOS_WATCH(rLeftHandSideMatrix)
    // KRATOS_WATCH(rRightHandSideVector)

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
//This method is called from outside the element
//************************************************************************************
//************************************************************************************
void TrussElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************
//This method is called from outside the element
//************************************************************************************
//************************************************************************************

void TrussElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************
//This method is called from outside the element
//************************************************************************************
//************************************************************************************

void TrussElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if(rMassMatrix.size1() != 6 || rMassMatrix.size2() != 6)
        rMassMatrix.resize(6, 6, false);
    noalias(rMassMatrix) = ZeroMatrix(6, 6);

    double lx = GetGeometry()[1].X0() - GetGeometry()[0].X0();
    double ly = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
    double lz = GetGeometry()[1].Z0() - GetGeometry()[0].Z0();
    double l = sqrt(lx*lx + ly*ly + lz*lz);
    double A = GetProperties()[AREA];
    double rho = GetProperties()[DENSITY];
    double m = rho*l*A/6;

    subrange(rMassMatrix, 0, 3, 0, 3) = 2.0*m*IdentityMatrix(3);
    subrange(rMassMatrix, 0, 3, 3, 6) = m*IdentityMatrix(3);
    subrange(rMassMatrix, 3, 6, 0, 3) = m*IdentityMatrix(3);
    subrange(rMassMatrix, 3, 6, 3, 6) = 2.0*m*IdentityMatrix(3);
}

//************************************************************************************
//************************************************************************************
//This method is called from outside the element
//************************************************************************************
//************************************************************************************

void TrussElement::CalculateDampingMatrix(MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    double alpha = GetProperties()[RAYLEIGH_DAMPING_ALPHA];
    double beta = GetProperties()[RAYLEIGH_DAMPING_BETA];

    this->CalculateMassMatrix(rDampMatrix, rCurrentProcessInfo);

    rDampMatrix *= alpha;

    Matrix K;
    Vector dummy;
    this->CalculateAll(K, dummy, rCurrentProcessInfo, true, false);

    noalias(rDampMatrix) += beta*K;
}

//************************************************************************************
//************************************************************************************
// returns the used integration method
//************************************************************************************
//************************************************************************************
TrussElement::IntegrationMethod TrussElement::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_1;
}

//************************************************************************************
//Informations to assemble the global vectors and matrices
//************************************************************************************

void TrussElement::EquationIdVector(EquationIdVectorType& rResult,
                                    const ProcessInfo& CurrentProcessInfo) const
{
    if(rResult.size() != 6)
        rResult.resize(6);

    rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[2] = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
    rResult[3] = GetGeometry()[1].GetDof(DISPLACEMENT_X).EquationId();
    rResult[4] = GetGeometry()[1].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[5] = GetGeometry()[1].GetDof(DISPLACEMENT_Z).EquationId();
}

//************************************************************************************
//************************************************************************************

void TrussElement::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo&
                              CurrentProcessInfo) const
{
    ElementalDofList.resize(0);

    ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_X));
    ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Y));
    ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Z));
    ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_X));
    ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Y));
    ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Z));
}

//************************************************************************************
//************************************************************************************

void TrussElement::GetValuesVector(Vector& values, int Step) const
{
    if(values.size() != 6)
        values.resize(6, false);

    values(0) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X, Step);
    values(1) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y, Step);
    values(2) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z, Step);
    values(3) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X, Step);
    values(4) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y, Step);
    values(5) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z, Step);
}

//************************************************************************************
//************************************************************************************

void TrussElement::GetFirstDerivativesVector(Vector& values, int Step) const
{
    if(values.size() != 6)
        values.resize(6, false);

    values(0) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_DT_X, Step);
    values(1) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_DT_Y, Step);
    values(2) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_DT_Z, Step);
    values(3) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_DT_X, Step);
    values(4) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_DT_Y, Step);
    values(5) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_DT_Z, Step);
}

//************************************************************************************
//************************************************************************************

void TrussElement::GetSecondDerivativesVector(Vector& values, int Step) const
{
    if(values.size() != 6)
        values.resize(6, false);

    values(0) = GetGeometry()[0].GetSolutionStepValue(ACCELERATION_X, Step);
    values(1) = GetGeometry()[0].GetSolutionStepValue(ACCELERATION_Y, Step);
    values(2) = GetGeometry()[0].GetSolutionStepValue(ACCELERATION_Z, Step);
    values(3) = GetGeometry()[1].GetSolutionStepValue(ACCELERATION_X, Step);
    values(4) = GetGeometry()[1].GetSolutionStepValue(ACCELERATION_Y, Step);
    values(5) = GetGeometry()[1].GetSolutionStepValue(ACCELERATION_Z, Step);
}

} // Namespace Kratos


