/*
see license.txt
*/
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 6 Jul 2016 $
//   Last Modified by:    $Author: Marwan $
//   Date:                $Date: 23-11-2015 $
//   Last Modified by:    $Author: Gall $
//   Date:                $Date: 00-00-2015 $
//   Revision:            $Revision: 0.0 $
//
//
// System includes 

// External includes 

// Project includes 
#include "custom_conditions/embedded_node_lagrange_tying_condition.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/kratos_flags.h"
#include "structural_application.h"


namespace Kratos
{
    //************************************************************************************
    //************************************************************************************

    EmbeddedNodeLagrangeTyingCondition::EmbeddedNodeLagrangeTyingCondition()
    {
    }

    EmbeddedNodeLagrangeTyingCondition::EmbeddedNodeLagrangeTyingCondition( IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition( NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    EmbeddedNodeLagrangeTyingCondition::EmbeddedNodeLagrangeTyingCondition( IndexType NewId,
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
    Condition::Pointer EmbeddedNodeLagrangeTyingCondition::Create( IndexType NewId, GeometryType::Pointer pGeometry,
                        NodeType::Pointer& pSlaveNode, Element::Pointer& pParentElement, PointType& rSolidLocalPoint ) const
    {
        return Condition::Pointer( new EmbeddedNodeLagrangeTyingCondition(NewId, pGeometry,
                        pSlaveNode, pParentElement, rSolidLocalPoint));
    }

    /**
     * Destructor. Never to be called manually
     */
    EmbeddedNodeLagrangeTyingCondition::~EmbeddedNodeLagrangeTyingCondition()
    {
    }

    void EmbeddedNodeLagrangeTyingCondition::Initialize()
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    //************************************************************************************    
    //************************************************************************************
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void EmbeddedNodeLagrangeTyingCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, 
            ProcessInfo& rCurrentProcessInfo)
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
    void EmbeddedNodeLagrangeTyingCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                              VectorType& rRightHandSideVector, 
                                              ProcessInfo& rCurrentProcessInfo)
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * This function calculates all system contributions due to the contact problem
     * with regard to the current master and slave partners.
     * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node 
     */
    void EmbeddedNodeLagrangeTyingCondition::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                      VectorType& rRightHandSideVector, 
                                      ProcessInfo& rCurrentProcessInfo,
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

        //resizing the RHS & LHS
        unsigned int MasterSize = mpMasterElement->GetGeometry().size();
        unsigned int ndofs = 3*MasterSize + 6;
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

        if( !CalculateStiffnessMatrixFlag && !CalculateResidualVectorFlag )
            return;

        double Scale;
        if( !this->Has(LAGRANGE_SCALE) )
            Scale = 1.0e10;
        else
            Scale = this->GetValue(LAGRANGE_SCALE);

        // Calculate the relative displacement
        Vector uSlave = mpSlaveNode->GetSolutionStepValue(DISPLACEMENT);

        Vector uMaster = ZeroVector(3);
        Vector ShapeFunctionValuesOnMaster;
        ShapeFunctionValuesOnMaster = mpMasterElement->GetGeometry().ShapeFunctionsValues(ShapeFunctionValuesOnMaster, mLocalPoint);
        for( unsigned int node = 0; node < mpMasterElement->GetGeometry().size(); ++node )
            noalias(uMaster) += ShapeFunctionValuesOnMaster[node] * mpMasterElement->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

        Vector Rel_Disp = uSlave - uMaster;

        //----------------------------------------------
        // RHS
        for( unsigned int node = 0 ; node < MasterSize ; ++node )
        {
            rRightHandSideVector[3*node    ] -= ShapeFunctionValuesOnMaster[node] * mpSlaveNode->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_X);
            rRightHandSideVector[3*node + 1] -= ShapeFunctionValuesOnMaster[node] * mpSlaveNode->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Y);
            rRightHandSideVector[3*node + 2] -= ShapeFunctionValuesOnMaster[node] * mpSlaveNode->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Z);
        }
        rRightHandSideVector[3*MasterSize    ] += mpSlaveNode->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_X);
        rRightHandSideVector[3*MasterSize + 1] += mpSlaveNode->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Y);
        rRightHandSideVector[3*MasterSize + 2] += mpSlaveNode->GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Z);
        rRightHandSideVector[3*MasterSize + 3] += Scale * Rel_Disp[0];
        rRightHandSideVector[3*MasterSize + 4] += Scale * Rel_Disp[1];
        rRightHandSideVector[3*MasterSize + 5] += Scale * Rel_Disp[2];

        //----------------------------------------------
        // LHS
        for( unsigned int Dim = 0; Dim < 3 ; ++Dim )
        {
            for( unsigned int node = 0; node < MasterSize; ++node )
            {
                rLeftHandSideMatrix(3*node + Dim, 3*MasterSize + Dim + 3 ) += ShapeFunctionValuesOnMaster[node];
            }
            rLeftHandSideMatrix(3*MasterSize + Dim, 3*MasterSize + Dim + 3 ) -= 1.0;
            for( unsigned int node = 0; node < MasterSize; ++node )
            {
                rLeftHandSideMatrix(3*MasterSize + Dim + 3, 3*node + Dim) += Scale * ShapeFunctionValuesOnMaster[node];
            }
            rLeftHandSideMatrix(3*MasterSize + Dim + 3, 3*MasterSize + Dim) -= Scale;
        }
        KRATOS_CATCH("")
    } // END CalculateAll
    //************************************************************************************
    //************************************************************************************

    /**
    * Setting up the EquationIdVector for the current partners.    
    * All conditions are assumed to be defined in 3D space with 3 DOFs per node.
    * All Equation IDs are given Master first, Slave second
    */
    void EmbeddedNodeLagrangeTyingCondition::EquationIdVector( EquationIdVectorType& rResult, 
                                          ProcessInfo& CurrentProcessInfo)
    {
//        if( mpMasterElement->GetValue(IS_INACTIVE) || !mpMasterElement->Is(ACTIVE) )
//        {
//            rResult.resize(0);
//            return;
//        }

        //determining size of DOF list (3D space)
        unsigned int MasterSize = mpMasterElement->GetGeometry().size();
        unsigned int ndofs = 3*MasterSize + 6;
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
        rResult[index++] = mpSlaveNode->GetDof(LAGRANGE_DISPLACEMENT_X).EquationId();
        rResult[index++] = mpSlaveNode->GetDof(LAGRANGE_DISPLACEMENT_Y).EquationId();
        rResult[index++] = mpSlaveNode->GetDof(LAGRANGE_DISPLACEMENT_Z).EquationId();
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
    void EmbeddedNodeLagrangeTyingCondition::GetDofList( DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
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
        ConditionalDofList.push_back( mpSlaveNode->pGetDof(LAGRANGE_DISPLACEMENT_X) );
        ConditionalDofList.push_back( mpSlaveNode->pGetDof(LAGRANGE_DISPLACEMENT_Y) );
        ConditionalDofList.push_back( mpSlaveNode->pGetDof(LAGRANGE_DISPLACEMENT_Z) );
    } // END GetDofList
} // Namespace Kratos

