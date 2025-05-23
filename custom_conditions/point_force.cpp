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
//   Project Name:        Kratos
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2007-12-12 14:51:07 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes


// External includes


// Project includes
#include "geometries/point_3d.h"
#include "custom_conditions/point_force.h"
#include "structural_application_variables.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
template<int TDim, typename TNodeType>
PointForce<TDim, TNodeType>::PointForce(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

//************************************************************************************
//************************************************************************************
template<int TDim, typename TNodeType>
PointForce<TDim, TNodeType>::PointForce(IndexType NewId, typename GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

//************************************************************************************
//************************************************************************************
template<int TDim, typename TNodeType>
PointForce<TDim, TNodeType>::PointForce(IndexType NewId, typename GeometryType::PointType::Pointer const& pNode, typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, typename GeometryType::Pointer(new Point3D<typename GeometryType::PointType>( pNode ) ), pProperties)
{
}

template<int TDim, typename TNodeType>
typename PointForce<TDim, TNodeType>::BaseType::Pointer PointForce<TDim, TNodeType>::Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                        typename PropertiesType::Pointer pProperties) const
{
    return typename BaseType::Pointer(new PointForce(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template<int TDim, typename TNodeType>
typename PointForce<TDim, TNodeType>::BaseType::Pointer PointForce<TDim, TNodeType>::Create(IndexType NewId, typename GeometryType::Pointer pGeom,
                                        typename PropertiesType::Pointer pProperties) const
{
    return typename BaseType::Pointer(new PointForce(NewId, pGeom, pProperties));
}

template<int TDim, typename TNodeType>
PointForce<TDim, TNodeType>::~PointForce()
{
}

//************************************************************************************
//************************************************************************************
template<int TDim, typename TNodeType>
void PointForce<TDim, TNodeType>::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rRightHandSideVector.size() != TDim)
        rRightHandSideVector.resize(TDim, false);

    if (!this->GetGeometry()[0].SolutionStepsDataHas(VARSEL(DataType, FORCE)))
        KRATOS_ERROR << "FORCE is not assigned to node " << this->GetGeometry()[0].Id();

    const auto& force = this->GetGeometry()[0].GetSolutionStepValue(VARSEL(DataType, FORCE));
    if constexpr (TDim > 0) rRightHandSideVector[0] = force[0];
    if constexpr (TDim > 1) rRightHandSideVector[1] = force[1];
    if constexpr (TDim > 2) rRightHandSideVector[2] = force[2];

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
template<int TDim, typename TNodeType>
void PointForce<TDim, TNodeType>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rLeftHandSideMatrix.size1() != TDim)
        rLeftHandSideMatrix.resize(TDim, TDim, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrixType(TDim, TDim);

    if(rRightHandSideVector.size() != TDim)
        rRightHandSideVector.resize(TDim, false);

    if (!this->GetGeometry()[0].SolutionStepsDataHas(VARSEL(DataType, FORCE)))
        KRATOS_ERROR << "FORCE is not assigned to node " << this->GetGeometry()[0].Id();

    const auto& force = this->GetGeometry()[0].GetSolutionStepValue(VARSEL(DataType, FORCE));
    if constexpr (TDim > 0) rRightHandSideVector[0] = force[0];
    if constexpr (TDim > 1) rRightHandSideVector[1] = force[1];
    if constexpr (TDim > 2) rRightHandSideVector[2] = force[2];

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
template<int TDim, typename TNodeType>
void PointForce<TDim, TNodeType>::CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    rDampMatrix.resize(0, 0, false);
}

//************************************************************************************
//************************************************************************************
template<int TDim, typename TNodeType>
void PointForce<TDim, TNodeType>::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    rMassMatrix.resize(0, 0, false);
}

//************************************************************************************
//************************************************************************************
template<int TDim, typename TNodeType>
void PointForce<TDim, TNodeType>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    unsigned int index;
    rResult.resize(number_of_nodes*TDim);
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        index = i*TDim;
        if constexpr (TDim > 0)
            rResult[index] = (this->GetGeometry()[i].GetDof(VARSELC(DataType, DISPLACEMENT, X)).EquationId());
        if constexpr (TDim > 1)
            rResult[index+1] = (this->GetGeometry()[i].GetDof(VARSELC(DataType, DISPLACEMENT, Y)).EquationId());
        if constexpr (TDim > 2)
            rResult[index+2] = (this->GetGeometry()[i].GetDof(VARSELC(DataType, DISPLACEMENT, Z)).EquationId());
    }
}

//************************************************************************************
//************************************************************************************
template<int TDim, typename TNodeType>
void PointForce<TDim, TNodeType>::GetDofList(DofsVectorType& ConditionalDofList, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int number_of_nodes = this->GetGeometry().size();
    unsigned int index;
    ConditionalDofList.resize(this->GetGeometry().size()*TDim);
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        index = i*TDim;
        if constexpr (TDim > 0)
            ConditionalDofList[index] = (this->GetGeometry()[i].pGetDof(VARSELC(DataType, DISPLACEMENT, X)));
        if constexpr (TDim > 1)
            ConditionalDofList[index+1] = (this->GetGeometry()[i].pGetDof(VARSELC(DataType, DISPLACEMENT, Y)));
        if constexpr (TDim > 2)
            ConditionalDofList[index+2] = (this->GetGeometry()[i].pGetDof(VARSELC(DataType, DISPLACEMENT, Z)));
    }
}

//************************************************************************************
//************************************************************************************
template class PointForce<2, RealNode>;
template class PointForce<2, ComplexNode>;
template class PointForce<2, GComplexNode>;
template class PointForce<3, RealNode>;
template class PointForce<3, ComplexNode>;
template class PointForce<3, GComplexNode>;

} // Namespace Kratos
