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
//
//   Project Name:        KratosStructuralApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 4 Apr 2026 $
//
//


// System includes


// External includes


// Project includes
#include "geometries/point_3d.h"
#include "custom_conditions/rotational_springs.h"
#include "structural_application_variables.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
RotationalSprings::RotationalSprings(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

RotationalSprings::RotationalSprings(IndexType NewId, typename GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

RotationalSprings::RotationalSprings(IndexType NewId, typename NodeType::Pointer const& pNode, typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, typename GeometryType::Pointer( new Point3D<NodeType>( pNode ) ), pProperties)
{
}

RotationalSprings::~RotationalSprings()
{
}

//************************************************************************************
//************************************************************************************
typename RotationalSprings::ConditionType::Pointer RotationalSprings::Create(IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const
{
    return typename BaseType::Pointer(new RotationalSprings(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

typename RotationalSprings::ConditionType::Pointer RotationalSprings::Create(IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const
{
    return typename BaseType::Pointer(new RotationalSprings(NewId, pGeom, pProperties));
}

//************************************************************************************
//************************************************************************************
void RotationalSprings::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void RotationalSprings::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void RotationalSprings::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                                      const ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag,
                                                      bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    unsigned int number_of_nodes = this->GetGeometry().size();
    unsigned int dim = 3;

    //resizing LHS and RHS where needed
    if (CalculateStiffnessMatrixFlag)
    {
        if(rLeftHandSideMatrix.size1() != number_of_nodes*dim || rLeftHandSideMatrix.size2() != number_of_nodes*dim)
            rLeftHandSideMatrix.resize(number_of_nodes*dim, number_of_nodes*dim,false);
        noalias(rLeftHandSideMatrix) = ZeroMatrixType(number_of_nodes*dim, number_of_nodes*dim);
    }

    if (CalculateResidualVectorFlag)
    {
        if(rRightHandSideVector.size() != number_of_nodes*dim)
            rRightHandSideVector.resize(number_of_nodes*dim, false);
        noalias(rRightHandSideVector) = ZeroVectorType(number_of_nodes*dim);
    }

    const auto& bedding = this->GetProperties()[ELASTIC_BEDDING_STIFFNESS];

    for (unsigned int i = 0; i < number_of_nodes; ++i )
    {
        const auto& rotation = this->GetGeometry()[i].GetSolutionStepValue(VARSEL(DataType, ROTATION));

        if (CalculateStiffnessMatrixFlag)
        {
            for( unsigned int j = 0; j < dim; j++ )
            {
                rLeftHandSideMatrix(dim*i+j, dim*i+j) = bedding[j];
            }
        }

        if (CalculateResidualVectorFlag)
        {
            for( unsigned int j = 0; j < dim; j++ )
            {
                rRightHandSideVector[dim*i+j] = -rotation[j]*bedding[j];
            }
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void RotationalSprings::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int number_of_nodes = this->GetGeometry().size();
    unsigned int index;
    unsigned int dim = 3;

    if(rResult.size() != number_of_nodes*dim)
        rResult.resize(number_of_nodes*dim);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        index = i*dim;
        rResult[index] = (this->GetGeometry()[i].GetDof(VARSELC(DataType, ROTATION, X)).EquationId());
        rResult[index+1] = (this->GetGeometry()[i].GetDof(VARSELC(DataType, ROTATION, Y)).EquationId());
        rResult[index+2] = (this->GetGeometry()[i].GetDof(VARSELC(DataType, ROTATION, Z)).EquationId());
    }
}

//************************************************************************************
//************************************************************************************
void RotationalSprings::GetDofList(DofsVectorType& rConditionalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    unsigned int number_of_nodes = this->GetGeometry().size();
    unsigned int index;
    unsigned int dim = 3;

    if(rConditionalDofList.size() != number_of_nodes*dim)
        rConditionalDofList.resize(number_of_nodes*dim);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        index = i*dim;
        rConditionalDofList[index] = (this->GetGeometry()[i].pGetDof(VARSELC(DataType, ROTATION, X)));
        rConditionalDofList[index+1] = (this->GetGeometry()[i].pGetDof(VARSELC(DataType, ROTATION, Y)));
        rConditionalDofList[index+2] = (this->GetGeometry()[i].pGetDof(VARSELC(DataType, ROTATION, Z)));
    }
}

} // Namespace Kratos
