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
*   Date:                $Date: 15 Dec 2017 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "geometries/point_3d.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_conditions/inclined_constraint_3d.h"
#include "structural_application_variables.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
InclinedConstraint3D::InclinedConstraint3D( IndexType NewId,
        GeometryType::Pointer pGeometry) :
    Condition( NewId, pGeometry )
{
}

InclinedConstraint3D::InclinedConstraint3D( IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) :
    Condition( NewId, pGeometry, pProperties )
{
}

InclinedConstraint3D::InclinedConstraint3D( IndexType NewId, Node<3>::Pointer const& node,
        const double& a, const double& b, const double& c, const double& d,
        PropertiesType::Pointer pProperties )
    : Condition( NewId, GeometryType::Pointer( new Point3D<Node<3> >( node ) ), pProperties )
    , mA(a), mB(b), mC(c), mD(d)
{
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************

Condition::Pointer InclinedConstraint3D::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new InclinedConstraint3D(NewId, GetGeometry().Create(ThisNodes),
                                pProperties));
}

Condition::Pointer InclinedConstraint3D::Create( IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new InclinedConstraint3D(NewId, pGeom,
                                pProperties));
}

/**
 * Destructor. Never to be called manually
 */
InclinedConstraint3D::~InclinedConstraint3D()
{
}


//************************************************************************************
//************************************************************************************

/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void InclinedConstraint3D::CalculateRightHandSide( VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{

    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag  = true;

    MatrixType matrix = Matrix();
    CalculateAll(matrix, rRightHandSideVector,
                 rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************

/**
 * calculates this contact element's local contributions
 */
void InclinedConstraint3D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag  = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
                 rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************
/**
 * calculates the contact related contributions to the system
 * Does nothing as assembling is to be switched to linking objects
 */
void InclinedConstraint3D::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag)
{
    unsigned int MatSize = 4;

    //resize as needed the LHS
    if ( CalculateResidualVectorFlag == true || CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        //resize the RHS=force vector if its size is not correct
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    double dX = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X);
    double dY = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y);
    double dZ = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z);

    double lambda = GetGeometry()[0].GetSolutionStepValue(LAGRANGE_MULTIPLIER_CONSTRAINT);

    rRightHandSideVector(0) = lambda * mA;
    rRightHandSideVector(1) = lambda * mB;
    rRightHandSideVector(2) = lambda * mC;
    rRightHandSideVector(3) = mA*dX + mB*dY + mC*dZ + mD;

    rLeftHandSideMatrix(0, 0) = 0.0;
    rLeftHandSideMatrix(0, 1) = 0.0;
    rLeftHandSideMatrix(0, 2) = 0.0;
    rLeftHandSideMatrix(0, 3) = -mA;

    rLeftHandSideMatrix(1, 0) = 0.0;
    rLeftHandSideMatrix(1, 1) = 0.0;
    rLeftHandSideMatrix(1, 2) = 0.0;
    rLeftHandSideMatrix(1, 3) = -mB;

    rLeftHandSideMatrix(2, 0) = 0.0;
    rLeftHandSideMatrix(2, 1) = 0.0;
    rLeftHandSideMatrix(2, 2) = 0.0;
    rLeftHandSideMatrix(2, 3) = -mC;

    rLeftHandSideMatrix(3, 0) = -mA;
    rLeftHandSideMatrix(3, 1) = -mB;
    rLeftHandSideMatrix(3, 2) = -mC;
    rLeftHandSideMatrix(3, 3) = 0.0;

//    KRATOS_WATCH( rLeftHandSideMatrix )
//    KRATOS_WATCH( rRightHandSideVector )
}

//************************************************************************************
//************************************************************************************

void InclinedConstraint3D::EquationIdVector( EquationIdVectorType& rResult,
        const ProcessInfo& CurrentProcessInfo) const
{
    //determining size of DOF list
    //dimension of space
    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    rResult.resize(GetGeometry().size()*dim+1, false);
    for( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        rResult[dim*i  ] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[dim*i+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[dim*i+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }
    rResult[GetGeometry().size()*dim] = GetGeometry()[0].GetDof(LAGRANGE_MULTIPLIER_CONSTRAINT).EquationId();
}

//************************************************************************************
//************************************************************************************

void InclinedConstraint3D::GetDofList( DofsVectorType& ConditionalDofList,
                                        const ProcessInfo& CurrentProcessInfo) const
{
    //determining size of DOF list
    //dimension of space
    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    ConditionalDofList.resize(GetGeometry().size()*dim+1);
    for( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        ConditionalDofList[dim*i  ] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        ConditionalDofList[dim*i+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        ConditionalDofList[dim*i+2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
    }
    ConditionalDofList[GetGeometry().size()*dim] = GetGeometry()[0].pGetDof(LAGRANGE_MULTIPLIER_CONSTRAINT);
}

} // Namespace Kratos
