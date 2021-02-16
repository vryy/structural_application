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
*   Date:                $Date: 8 Jun 2017 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/point_bedding_condition.h"
#include "structural_application_variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "geometries/point_3d.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
PointBeddingCondition::PointBeddingCondition( IndexType NewId,
        GeometryType::Pointer pGeometry) :
    Condition( NewId, pGeometry )
{
}

PointBeddingCondition::PointBeddingCondition( IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) :
    Condition( NewId, pGeometry, pProperties )
{
}

PointBeddingCondition::PointBeddingCondition( IndexType NewId, Node<3>::Pointer const& pNode,
        PropertiesType::Pointer pProperties )
    : Condition( NewId, GeometryType::Pointer( new Point3D<Node<3> >( pNode ) ), pProperties )
{
}

PointBeddingCondition::PointBeddingCondition( IndexType NewId, NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties )
    : Condition( NewId, GeometryType::Pointer( new Point3D<Node<3> >( ThisNodes ) ), pProperties )
{
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************

Condition::Pointer PointBeddingCondition::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new PointBeddingCondition(NewId, GetGeometry().Create(ThisNodes),
                                pProperties));
}

Condition::Pointer PointBeddingCondition::Create( IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new PointBeddingCondition(NewId, pGeom,
                                pProperties));
}

/**
 * Destructor. Never to be called manually
 */
PointBeddingCondition::~PointBeddingCondition()
{
}

//************************************************************************************
//************************************************************************************

/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void PointBeddingCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
void PointBeddingCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
void PointBeddingCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag)
{
    unsigned int MatSize = 3;

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

    const Vector& stiffVect = this->GetValue(ELASTIC_BEDDING_STIFFNESS);
    rLeftHandSideMatrix(0, 0) = stiffVect(0);
    rLeftHandSideMatrix(1, 1) = stiffVect(1);
    rLeftHandSideMatrix(2, 2) = stiffVect(2);

    //compute internal forces
    Vector displacements = ZeroVector(6);
    for( unsigned int i = 0; i < 2; i++ )
    {
        displacements[3*i  ] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X);
        displacements[3*i+1] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y);
        displacements[3*i+2] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z);
    }
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, displacements);

//    KRATOS_WATCH( rLeftHandSideMatrix );
//    KRATOS_WATCH( rRightHandSideVector );
}

//************************************************************************************
//************************************************************************************

//     void PointBeddingCondition::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo )
//     {
//         KRATOS_WATCH( rValue );
//         if( rThisVariable == JOINT_STIFFNESS )
//         {
//             if( rValue.size1() != 6 || rValue.size2() != 6 )
//                 KRATOS_THROW_ERROR( std::logic_error, "Stiffness of the joint must be 6x6", "" );
//             noalias(mStiffnessMatrix) = rValue;
//         }
//     }

//************************************************************************************
//************************************************************************************

void PointBeddingCondition::EquationIdVector( EquationIdVectorType& rResult,
        const ProcessInfo& CurrentProcessInfo ) const
{

    //determining size of DOF list
    //dimension of space
    unsigned int dim = 3;
    rResult.resize(2*dim,false);
    for( unsigned int i=0; i<2; i++ )
    {
        rResult[dim*i] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[dim*i+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[dim*i+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
void PointBeddingCondition::GetDofList( DofsVectorType& ConditionalDofList,
                                        const ProcessInfo& CurrentProcessInfo) const
{
//determining size of DOF list
    //dimension of space
    unsigned int dim = 3;
    ConditionalDofList.resize(2*dim);
    for( unsigned int i=0; i<2; i++ )
    {
        ConditionalDofList[dim*i] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        ConditionalDofList[dim*i+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        ConditionalDofList[dim*i+2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
    }
}

} // Namespace Kratos
