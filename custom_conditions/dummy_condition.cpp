//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 1 Mar 2017$
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_conditions/dummy_condition.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
DummyCondition::DummyCondition()
{
}

DummyCondition::DummyCondition( IndexType NewId,
                              GeometryType::Pointer pGeometry)
: Condition( NewId, pGeometry )
{
}

DummyCondition::DummyCondition( IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Condition( NewId, pGeometry, pProperties )
{
}

/**
 * Destructor. Never to be called manually
 */
DummyCondition::~DummyCondition()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Condition::Pointer DummyCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new DummyCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Condition::Pointer DummyCondition::Create(IndexType NewId, GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new DummyCondition(NewId, pGeom, pProperties));
}

void DummyCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void DummyCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType matrix = Matrix();
    CalculateAll( matrix, rRightHandSideVector,
                  rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag,
                  CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************
void DummyCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************    /
void DummyCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  const ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    rLeftHandSideMatrix.resize(0, 0, false);
    rRightHandSideVector.resize(0, false);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************    /
void DummyCondition::CalculateDampingMatrix( MatrixType& rDampingMatrix,
                                             const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************    /
void DummyCondition::CalculateMassMatrix( MatrixType& rMassMatrix,
                                          const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void DummyCondition::EquationIdVector( EquationIdVectorType& rResult,
                                       const ProcessInfo& CurrentProcessInfo) const
{
    rResult.resize(0);
}

//************************************************************************************
//************************************************************************************
void DummyCondition::GetDofList( DofsVectorType& ConditionalDofList,
                                 const ProcessInfo& CurrentProcessInfo) const
{
    ConditionalDofList.resize(0);
}

//************************************************************************************
//************************************************************************************
int DummyCondition::Check( const ProcessInfo& CurrentProcessInfo) const
{
    return 0;
}

} // Namespace Kratos

