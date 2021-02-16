/*
see license.txt
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 6 Jul 2016 $
//
//
// System includes

// External includes

// Project includes
#include "custom_conditions/embedded_node_penalty_tying_condition.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/kratos_flags.h"
#include "structural_application_variables.h"


namespace Kratos
{
    //************************************************************************************
    //************************************************************************************

    EmbeddedNodePenaltyTyingCondition::EmbeddedNodePenaltyTyingCondition()
    {
    }

    EmbeddedNodePenaltyTyingCondition::EmbeddedNodePenaltyTyingCondition( IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition( NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    EmbeddedNodePenaltyTyingCondition::EmbeddedNodePenaltyTyingCondition( IndexType NewId,
                GeometryType::Pointer pGeometry,
                NodeType::Pointer& pSlaveNode,
                Element::Pointer& pParentElement,
                PointType& rLocalPoint )
    : Condition( NewId, pGeometry ), mpSlaveNode(pSlaveNode), mpMasterElement(pParentElement), mLocalPoint(rLocalPoint)
    {
    }

    //********************************************************
    //**** Operations ****************************************
    //********************************************************
    Condition::Pointer EmbeddedNodePenaltyTyingCondition::Create( IndexType NewId, GeometryType::Pointer pGeometry,
                        NodeType::Pointer& pSlaveNode, Element::Pointer& pParentElement, PointType& rSolidLocalPoint ) const
    {
        return Condition::Pointer( new EmbeddedNodePenaltyTyingCondition(NewId, pGeometry,
                        pSlaveNode, pParentElement, rSolidLocalPoint));
    }

    /**
     * Destructor. Never to be called manually
     */
    EmbeddedNodePenaltyTyingCondition::~EmbeddedNodePenaltyTyingCondition()
    {
    }

    void EmbeddedNodePenaltyTyingCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void EmbeddedNodePenaltyTyingCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo)
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType dummy;
        CalculateAll( dummy, rRightHandSideVector,
                      rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag,
                      CalculateResidualVectorFlag);
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * calculates this contact element's local contributions
     */
    void EmbeddedNodePenaltyTyingCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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

    /**
     * This function calculates all system contributions due to the contact problem
     * with regard to the current master and slave partners.
     * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node
     */
    void EmbeddedNodePenaltyTyingCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

//        if( mpMasterElement->GetValue(IS_INACTIVE) || !mpMasterElement->Is(ACTIVE) )
//        {
//            rRightHandSideVector.resize(0, false);
//            rLeftHandSideMatrix.resize(0, 0, false);
//            return;
//        }

        if( !CalculateStiffnessMatrixFlag && !CalculateResidualVectorFlag )
            return;

        if( !this->Has(INITIAL_PENALTY) )
            KRATOS_THROW_ERROR(std::logic_error, "INITIAL_PENALTY is not set for embbedded tying link", Id())
        double Penalty = this->GetValue(INITIAL_PENALTY);

        //resizing the RHS & LHS
        unsigned int MasterSize = mpMasterElement->GetGeometry().size();
        unsigned int ndofs = 3*MasterSize + 3;
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if(rRightHandSideVector.size() != ndofs)
                rRightHandSideVector.resize(ndofs, false);
            noalias(rRightHandSideVector) = ZeroVector(ndofs); //resetting RHS    */
        }
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if(rLeftHandSideMatrix.size1() != ndofs || rLeftHandSideMatrix.size2() != ndofs)
                rLeftHandSideMatrix.resize(ndofs, ndofs, false);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(ndofs,ndofs);
        }

        // Calculate the relative displacement
        Vector uSlave = mpSlaveNode->GetSolutionStepValue(DISPLACEMENT);

        Vector uMaster = ZeroVector(3);
        Vector ShapeFunctionValuesOnMaster;
        ShapeFunctionValuesOnMaster = mpMasterElement->GetGeometry().ShapeFunctionsValues(ShapeFunctionValuesOnMaster, mLocalPoint);
        for( unsigned int node = 0; node < mpMasterElement->GetGeometry().size(); ++node )
            noalias(uMaster) += ShapeFunctionValuesOnMaster[node] * mpMasterElement->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

        Vector Rel_Disp = uSlave - uMaster;
//        KRATOS_WATCH(Rel_Disp)

        //----------------------------------------------
        // RHS
        for( unsigned int node = 0 ; node < MasterSize ; ++node )
        {
            rRightHandSideVector[3*node]     += Penalty * Rel_Disp[0];
            rRightHandSideVector[3*node + 1] += Penalty * Rel_Disp[1];
            rRightHandSideVector[3*node + 2] += Penalty * Rel_Disp[2];
        }
        rRightHandSideVector[3*MasterSize    ] -= Penalty * Rel_Disp[0];
        rRightHandSideVector[3*MasterSize + 1] -= Penalty * Rel_Disp[1];
        rRightHandSideVector[3*MasterSize + 2] -= Penalty * Rel_Disp[2];
//        KRATOS_WATCH(rRightHandSideVector)

        //----------------------------------------------
        // LHS
        unsigned int row, col;
        for( unsigned int node = 0; node < MasterSize; ++node )
        {
            for( unsigned int dim = 0; dim < 3; ++dim )
            {
                row = 3*node + dim;
                for( unsigned int node2 = 0; node2 < MasterSize; ++node2 )
                {
                    col = 3*node2 + dim;
                    rLeftHandSideMatrix(row, col) = Penalty * ShapeFunctionValuesOnMaster[node2];
                }
                col = 3*MasterSize + dim;
                rLeftHandSideMatrix(row, col) = -Penalty;
            }
        }
        for( unsigned int dim = 0; dim < 3; ++dim )
        {
            row = 3*MasterSize + dim;
            for( unsigned int node2 = 0; node2 < MasterSize; ++node2 )
            {
                col = 3*node2 + dim;
                rLeftHandSideMatrix(row, col) = -Penalty * ShapeFunctionValuesOnMaster[node2];
            }
            col = 3*MasterSize + dim;
            rLeftHandSideMatrix(row, col) = Penalty;
        }
//        KRATOS_WATCH(rLeftHandSideMatrix)

        KRATOS_CATCH("")
    } // END CalculateAll

    //************************************************************************************
    //************************************************************************************
    /**
    * Setting up the EquationIdVector for the current partners.
    * All conditions are assumed to be defined in 3D space with 3 DOFs per node.
    * All Equation IDs are given Master first, Slave second
    */
    void EmbeddedNodePenaltyTyingCondition::EquationIdVector( EquationIdVectorType& rResult,
                                          const ProcessInfo& CurrentProcessInfo) const
    {
//        if( mpMasterElement->GetValue(IS_INACTIVE) || !mpMasterElement->Is(ACTIVE) )
//        {
//            rResult.resize(0);
//            return;
//        }

        //determining size of DOF list (3D space)
        unsigned int MasterSize = mpMasterElement->GetGeometry().size();
        unsigned int ndofs = 3*MasterSize + 3;
        unsigned int index;

        if(rResult.size() != ndofs)
            rResult.resize(ndofs,false);

        index = 0;
        for( unsigned int node = 0; node < MasterSize ; ++node )
        {
            rResult[index++] = mpMasterElement->GetGeometry()[node].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index++] = mpMasterElement->GetGeometry()[node].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index++] = mpMasterElement->GetGeometry()[node].GetDof(DISPLACEMENT_Z).EquationId();
        }
        rResult[index++] = mpSlaveNode->GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = mpSlaveNode->GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = mpSlaveNode->GetDof(DISPLACEMENT_Z).EquationId();
    } // END EquationIdVector
    //************************************************************************************
    //************************************************************************************

    /**
     * Setting up the DOF list for the current partners.
     * All conditions are assumed to be defined in 3D space with 3 DOFs per Node.
     * All DOF are given Master first, Slave second
     */
    //************************************************************************************
    //************************************************************************************
    void EmbeddedNodePenaltyTyingCondition::GetDofList( DofsVectorType& ConditionalDofList, const ProcessInfo& CurrentProcessInfo) const
    {
        ConditionalDofList.resize(0);
//        if( mpMasterElement->GetValue(IS_INACTIVE) || !mpMasterElement->Is(ACTIVE) )
//        {
//            return;
//        }

        unsigned int MasterSize = mpMasterElement->GetGeometry().size();

        for( unsigned int node = 0; node < MasterSize ; ++node )
        {
            ConditionalDofList.push_back( mpMasterElement->GetGeometry()[node].pGetDof(DISPLACEMENT_X) );
            ConditionalDofList.push_back( mpMasterElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Y) );
            ConditionalDofList.push_back( mpMasterElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Z) );
        }
        ConditionalDofList.push_back( mpSlaveNode->pGetDof(DISPLACEMENT_X) );
        ConditionalDofList.push_back( mpSlaveNode->pGetDof(DISPLACEMENT_Y) );
        ConditionalDofList.push_back( mpSlaveNode->pGetDof(DISPLACEMENT_Z) );
    } // END GetDofList
} // Namespace Kratos

