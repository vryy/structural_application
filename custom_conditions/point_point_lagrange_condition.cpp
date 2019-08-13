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
*   Last Modified by:    $Author: Jelena $
*   Date:                $Date: 2014 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/point_point_lagrange_condition.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "geometries/line_3d_2.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
PointPointLagrangeCondition::PointPointLagrangeCondition( IndexType NewId,
        GeometryType::Pointer pGeometry) :
    Condition( NewId, pGeometry )
{
}

PointPointLagrangeCondition::PointPointLagrangeCondition( IndexType NewId,
        GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
    Condition( NewId, pGeometry, pProperties )
{
}

PointPointLagrangeCondition::PointPointLagrangeCondition(
    IndexType NewId,
    Node<3>::Pointer const& node1,
    Node<3>::Pointer const& node2,
    PropertiesType::Pointer pProperties
) : Condition( NewId, GeometryType::Pointer( new Line3D2<Node<3> >( node1, node2 ) ), pProperties )
{
}

PointPointLagrangeCondition::PointPointLagrangeCondition(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
) : Condition( NewId, GeometryType::Pointer( new Line3D2<Node<3> >( ThisNodes ) ), pProperties )
{
}


//
Condition::Pointer PointPointLagrangeCondition::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new PointPointLagrangeCondition(NewId, GetGeometry().Create(ThisNodes),
                                pProperties));
}




/**
 * Destructor. Never to be called manually
 */
PointPointLagrangeCondition::~PointPointLagrangeCondition()
{
}



//************************************************************************************
//************************************************************************************

/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void PointPointLagrangeCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
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
void PointPointLagrangeCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
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
void PointPointLagrangeCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag)
{
    unsigned int MatSize = 9;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        //resize the LHS=StiffnessMatrix if its size is not correct
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

    Vector displacements = ZeroVector(6);
    for( unsigned int i=0; i<1; i++ )
    {
        displacements[3*i] = GetGeometry()[i].X();
        displacements[3*i+1] = GetGeometry()[i].Y();
        displacements[3*i+2] = GetGeometry()[i].Z();
    }

    for( unsigned int i=1; i<2; i++ )
    {
        displacements[3*i] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X);
        displacements[3*i+1] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y);
        displacements[3*i+2] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z);
    }

	Vector rel_Disp = ZeroVector(3);
    for( unsigned int i=0; i<2; i++ )
    {
        rel_Disp[0] = displacements[0]-displacements[3];
        rel_Disp[1] = displacements[1]-displacements[4];
        rel_Disp[2] = displacements[2]-displacements[5];
    }
 //       KRATOS_WATCH (rel_Disp)

	Vector relDisp = ZeroVector(3);
	noalias(relDisp) += GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT) - GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT);
 //       KRATOS_WATCH (relDisp)

    double mStiffness = GetProperties()[LINING_JOINT_STIFFNESS];
	rRightHandSideVector[0] -= relDisp[0];
	rRightHandSideVector[1] -= relDisp[1];
	rRightHandSideVector[2] -= relDisp[2];
	rRightHandSideVector[3] += relDisp[0];
	rRightHandSideVector[3+1] += relDisp[1];
	rRightHandSideVector[3+2] += relDisp[2];
	rRightHandSideVector[2*3] += relDisp[0]/mStiffness;
	rRightHandSideVector[2*3+1] += relDisp[1]/mStiffness;
	rRightHandSideVector[2*3+2] += relDisp[2]/mStiffness;

    for ( unsigned int i = 0; i < 3; i++ )
    {
        rLeftHandSideMatrix(  i,  2*3 + i ) -= 1.0;
        rLeftHandSideMatrix( 2*3 + i,  i ) -= 1.0;
        rLeftHandSideMatrix( 3 + i,  2*3  + i ) += 1.0;
        rLeftHandSideMatrix( 2*3 + i, 3 + i ) += 1.0;
    }
}

void PointPointLagrangeCondition::EquationIdVector( EquationIdVectorType& rResult,
        ProcessInfo& CurrentProcessInfo)
{

    //determining size of DOF list
    //dimension of space
    unsigned int dim = 3;
	unsigned int index;
    rResult.resize(3*dim,false);
    for( unsigned int i=0; i<2; i++ )
    {
        rResult[dim*i] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[dim*i+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[dim*i+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

    }

    index = dim*2;
	rResult[index] = GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_X).EquationId();
	rResult[index+1] = GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Y).EquationId();
	rResult[index+2] = GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Z).EquationId();
}

//************************************************************************************
//************************************************************************************
void PointPointLagrangeCondition::GetDofList( DofsVectorType& ConditionalDofList,
                                        ProcessInfo& CurrentProcessInfo)
{
//determining size of DOF list
    //dimension of space
    unsigned int dim = 3;
	unsigned int index;
    ConditionalDofList.resize(3*dim);
    for( unsigned int i=0; i<2; i++ )
    {
        ConditionalDofList[dim*i] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        ConditionalDofList[dim*i+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        ConditionalDofList[dim*i+2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);

    }

    index = dim*2;
    ConditionalDofList[index] = GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_X);
    ConditionalDofList[index+1] = GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Y);
    ConditionalDofList[index+2] = GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Z);
}

} // Namespace Kratos
