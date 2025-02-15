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
*   Date:                $Date: 17 Oct 2017 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "geometries/line_2d_3.h"
#include "structural_application_variables.h"
#include "custom_conditions/collinear_constraint_2d.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
CollinearConstraint2D::CollinearConstraint2D( IndexType NewId,
        GeometryType::Pointer pGeometry) :
    Condition( NewId, pGeometry )
{
}

CollinearConstraint2D::CollinearConstraint2D( IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) :
    Condition( NewId, pGeometry, pProperties )
{
}

CollinearConstraint2D::CollinearConstraint2D( IndexType NewId, Node<3>::Pointer const& nodeA, Node<3>::Pointer const& nodeC, Node<3>::Pointer const& nodeB,
        PropertiesType::Pointer pProperties )
    : Condition( NewId, GeometryType::Pointer( new Line2D3<Node<3> >( nodeA, nodeB, nodeC ) ), pProperties )
{
}

CollinearConstraint2D::CollinearConstraint2D( IndexType NewId, NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties )
    : Condition( NewId, GeometryType::Pointer( new Line2D3<Node<3> >( ThisNodes ) ), pProperties )
{
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************

Condition::Pointer CollinearConstraint2D::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new CollinearConstraint2D(NewId, GetGeometry().Create(ThisNodes),
                                pProperties));
}

Condition::Pointer CollinearConstraint2D::Create( IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new CollinearConstraint2D(NewId, pGeom,
                                pProperties));
}

/**
 * Destructor. Never to be called manually
 */
CollinearConstraint2D::~CollinearConstraint2D()
{
}

//************************************************************************************
//************************************************************************************
void CollinearConstraint2D::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
void CollinearConstraint2D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
void CollinearConstraint2D::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag)
{
    unsigned int MatSize = 7;

    //resize as needed the LHS
    if ( CalculateResidualVectorFlag || CalculateStiffnessMatrixFlag ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) //calculation of the matrix is required
    {
        //resize the RHS=force vector if its size is not correct
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    //compute constraint values
    double xA = GetGeometry()[0].X0() + GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X);
    double yA = GetGeometry()[0].Y0() + GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y);

    double xB = GetGeometry()[1].X0() + GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X);
    double yB = GetGeometry()[1].Y0() + GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y);

    double xC = GetGeometry()[2].X0() + GetGeometry()[2].GetSolutionStepValue(DISPLACEMENT_X);
    double yC = GetGeometry()[2].Y0() + GetGeometry()[2].GetSolutionStepValue(DISPLACEMENT_Y);

    double lambda = GetGeometry()[2].GetSolutionStepValue(LAGRANGE_MULTIPLIER_CONSTRAINT);

    rRightHandSideVector(0) = lambda * (yB - yC);
    rRightHandSideVector(1) = lambda * (xC - xB);
    rRightHandSideVector(2) = lambda * (yC - yA);
    rRightHandSideVector(3) = lambda * (xA - xC);
    rRightHandSideVector(4) = lambda * (yA - yB);
    rRightHandSideVector(5) = lambda * (xB - xA);
    rRightHandSideVector(6) = (xA*yB + xB*yC + xC*yA) - (xB*yA + xC*yB + xA*yC);

    rLeftHandSideMatrix(0, 3) = -lambda;
    rLeftHandSideMatrix(0, 5) = lambda;
    rLeftHandSideMatrix(0, 6) = -(yB - yC);

    rLeftHandSideMatrix(1, 2) = lambda;
    rLeftHandSideMatrix(1, 4) = -lambda;
    rLeftHandSideMatrix(1, 6) = -(xC - xB);

    rLeftHandSideMatrix(2, 1) = lambda;
    rLeftHandSideMatrix(2, 5) = -lambda;
    rLeftHandSideMatrix(2, 6) = -(yC - yA);

    rLeftHandSideMatrix(3, 0) = -lambda;
    rLeftHandSideMatrix(3, 4) = lambda;
    rLeftHandSideMatrix(3, 6) = -(xA - xC);

    rLeftHandSideMatrix(4, 1) = -lambda;
    rLeftHandSideMatrix(4, 3) = lambda;
    rLeftHandSideMatrix(4, 6) = -(yA - yB);

    rLeftHandSideMatrix(5, 0) = lambda;
    rLeftHandSideMatrix(5, 2) = -lambda;
    rLeftHandSideMatrix(5, 6) = -(xB - xA);

    rLeftHandSideMatrix(6, 0) = -(yB - yC);
    rLeftHandSideMatrix(6, 1) = -(xC - xB);
    rLeftHandSideMatrix(6, 2) = -(yC - yA);
    rLeftHandSideMatrix(6, 3) = -(xA - xC);
    rLeftHandSideMatrix(6, 4) = -(yA - yB);
    rLeftHandSideMatrix(6, 5) = -(xB - xA);
    rLeftHandSideMatrix(6, 6) = 0.0;

//    KRATOS_WATCH( rLeftHandSideMatrix )
//    KRATOS_WATCH( rRightHandSideVector )
}

//************************************************************************************
//************************************************************************************

void CollinearConstraint2D::EquationIdVector( EquationIdVectorType& rResult,
        const ProcessInfo& CurrentProcessInfo ) const
{
    //determining size of DOF list
    rResult.resize(GetGeometry().size()*2+1, false);
    for( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        rResult[2*i] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[2*i+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
    }
    rResult[GetGeometry().size()*2] = GetGeometry()[2].GetDof(LAGRANGE_MULTIPLIER_CONSTRAINT).EquationId();
}

//************************************************************************************
//************************************************************************************

void CollinearConstraint2D::GetDofList( DofsVectorType& ConditionalDofList,
                                        const ProcessInfo& CurrentProcessInfo) const
{
    //determining size of DOF list
    ConditionalDofList.resize(GetGeometry().size()*2+1);
    for( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        ConditionalDofList[2*i] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        ConditionalDofList[2*i+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
    }
    ConditionalDofList[GetGeometry().size()*2] = GetGeometry()[2].pGetDof(LAGRANGE_MULTIPLIER_CONSTRAINT);
}

} // Namespace Kratos
