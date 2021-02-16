/*
==============================================================================
KratosR1StructuralApplication
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
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Jul 2016$
//   Last Modified by:    $Author: jelena $
//   Date:                $Date: 2015-03-25$
//   Revision:            $Revision: 1.2 $
//
//
// System includes

// External includes

// Project includes
#include "custom_conditions/embedded_point_lagrange_tying_condition.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/kratos_flags.h"
#include "structural_application_variables.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    EmbeddedPointLagrangeTyingCondition::EmbeddedPointLagrangeTyingCondition()
    {
    }

    EmbeddedPointLagrangeTyingCondition::EmbeddedPointLagrangeTyingCondition( IndexType NewId,
                                  GeometryType::Pointer pGeometry)
    : Condition( NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    EmbeddedPointLagrangeTyingCondition::EmbeddedPointLagrangeTyingCondition( IndexType NewId,
                                GeometryType::Pointer pGeometry,
                                Element::Pointer pMasterElement,
                                Element::Pointer pSlaveElement,
                                PointType& rMasterLocalPoint,
                                PointType& rSlaveLocalPoint )
    : Condition( NewId, pGeometry )
    {
        mMasterLocalPoint = rMasterLocalPoint;
        mSlaveLocalPoint = rSlaveLocalPoint;
        mpMasterElement = pMasterElement;
        mpSlaveElement = pSlaveElement;
    }

    //********************************************************
    //**** Operations ****************************************
    //********************************************************


    Condition::Pointer EmbeddedPointLagrangeTyingCondition::Create( IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       Element::Pointer pMasterElement,
                                       Element::Pointer pSlaveElement,
                                       PointType& rMasterLocalPoint,
                                       PointType& rSlaveLocalPoint ) const
    {
        return Condition::Pointer( new EmbeddedPointLagrangeTyingCondition(NewId, pGeometry, pMasterElement, pSlaveElement, rMasterLocalPoint, rSlaveLocalPoint) );
    }

    /**
     * Destructor. Never to be called manually
     */
    EmbeddedPointLagrangeTyingCondition::~EmbeddedPointLagrangeTyingCondition()
    {
    }


    void EmbeddedPointLagrangeTyingCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void EmbeddedPointLagrangeTyingCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
    /**
     * calculates this contact element's local contributions
     */
    void EmbeddedPointLagrangeTyingCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
    void EmbeddedPointLagrangeTyingCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

//        if(   ( mpMasterElement->GetValue(IS_INACTIVE) || !mpMasterElement->Is(ACTIVE) )
//           || ( mpSlaveElement->GetValue(IS_INACTIVE) || !mpSlaveElement->Is(ACTIVE) ) )
//        {
//            rRightHandSideVector.resize(0, false);
//            rLeftHandSideMatrix.resize(0, 0, false);
//            return;
//        }

        Vector uMaster = ZeroVector(3);
        Vector ShapeFunctionValuesOnMaster;
        ShapeFunctionValuesOnMaster = mpMasterElement->GetGeometry().ShapeFunctionsValues(ShapeFunctionValuesOnMaster, mMasterLocalPoint);
        for( unsigned int node = 0; node < mpMasterElement->GetGeometry().size(); node++ )
        {
            uMaster += ShapeFunctionValuesOnMaster[node] * mpMasterElement->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);
        }

        Vector uSlave = ZeroVector(3);
        Vector ShapeFunctionValuesOnSlave;
        ShapeFunctionValuesOnSlave = mpSlaveElement->GetGeometry().ShapeFunctionsValues(ShapeFunctionValuesOnSlave, mSlaveLocalPoint);
        for( unsigned int node = 0; node < mpSlaveElement->GetGeometry().size(); node++ )
        {
            uSlave += ShapeFunctionValuesOnSlave[node] * mpSlaveElement->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);
        }

        Vector relDisp = uSlave - uMaster;

        unsigned int MasterNN = mpMasterElement->GetGeometry().size();
        unsigned int SlaveNN = mpSlaveElement->GetGeometry().size();
        unsigned int MatSize = (MasterNN+SlaveNN)*3 + 3;
        KRATOS_WATCH(MatSize);

        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if(rRightHandSideVector.size() != MatSize)
                rRightHandSideVector.resize(MatSize);
            noalias(rRightHandSideVector) = ZeroVector(MatSize);
        }
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if(rLeftHandSideMatrix.size1() != MatSize || rLeftHandSideMatrix.size2() != MatSize)
                rLeftHandSideMatrix.resize(MatSize,MatSize);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize);
        }
        else return;

        for( unsigned int node = 0; node < MasterNN; node++ )
        {
            rRightHandSideVector[3*node    ] -= ShapeFunctionValuesOnMaster[node] * mpSlaveElement->GetGeometry()[0].GetSolutionStepValue(LAGRANGE_DISPLACEMENT_X);
            rRightHandSideVector[3*node + 1] -= ShapeFunctionValuesOnMaster[node] * mpSlaveElement->GetGeometry()[0].GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Y);
            rRightHandSideVector[3*node + 2] -= ShapeFunctionValuesOnMaster[node] * mpSlaveElement->GetGeometry()[0].GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Z);
        }
        for( unsigned int node = 0; node < SlaveNN; node++ )
        {
            rRightHandSideVector[3*MasterNN + 3*node    ] += ShapeFunctionValuesOnSlave[node] * mpSlaveElement->GetGeometry()[0].GetSolutionStepValue(LAGRANGE_DISPLACEMENT_X);
            rRightHandSideVector[3*MasterNN + 3*node + 1] += ShapeFunctionValuesOnSlave[node] * mpSlaveElement->GetGeometry()[0].GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Y);
            rRightHandSideVector[3*MasterNN + 3*node + 2] += ShapeFunctionValuesOnSlave[node] * mpSlaveElement->GetGeometry()[0].GetSolutionStepValue(LAGRANGE_DISPLACEMENT_Z);
        }
        rRightHandSideVector[3*MasterNN + 3*SlaveNN    ] += relDisp[0];
        rRightHandSideVector[3*MasterNN + 3*SlaveNN + 1] += relDisp[1];
        rRightHandSideVector[3*MasterNN + 3*SlaveNN + 2] += relDisp[2];

        for( unsigned int Dim = 0; Dim < 3 ; ++Dim )
        {
            for( unsigned int node = 0; node < MasterNN; ++node )
            {
                rLeftHandSideMatrix(3*node + Dim, 3*MasterNN + 3*SlaveNN + Dim ) += ShapeFunctionValuesOnMaster[node];
                rLeftHandSideMatrix(3*MasterNN + 3*SlaveNN + Dim, 3*node + Dim ) += ShapeFunctionValuesOnMaster[node];
            }
            for( unsigned int node = 0; node < SlaveNN; ++node )
            {
                rLeftHandSideMatrix(3*node + 3*MasterNN + Dim, 3*MasterNN + 3*SlaveNN + Dim ) -= ShapeFunctionValuesOnSlave[node];
                rLeftHandSideMatrix(3*MasterNN + 3*SlaveNN + Dim, 3*node + 3*MasterNN + Dim ) -= ShapeFunctionValuesOnSlave[node];
            }
        }

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************
    /**
    * Setting up the EquationIdVector for the current partners.
    * All conditions are assumed to be defined in 3D space with 3 DOFs per node.
    * All Equation IDs are given Master first, Slave second
    */
    void EmbeddedPointLagrangeTyingCondition::EquationIdVector( EquationIdVectorType& rResult,
                                          const ProcessInfo& CurrentProcessInfo) const
    {
//        if(   ( mpMasterElement->GetValue(IS_INACTIVE) || !mpMasterElement->Is(ACTIVE) )
//           || ( mpSlaveElement->GetValue(IS_INACTIVE) || !mpSlaveElement->Is(ACTIVE) ) )
//        {
//            rResult.resize(0);
//            return;
//        }

        //determining size of DOF list
        //dimension of space
        unsigned int MasterNN = mpMasterElement->GetGeometry().size();
        unsigned int SlaveNN = mpSlaveElement->GetGeometry().size();
        unsigned int ndofs = 3*(MasterNN+SlaveNN) + 3;
        unsigned int index;

        if(rResult.size() != ndofs)
            rResult.resize(ndofs, false);

        for( unsigned int node = 0; node < MasterNN; node++ )
        {
            index = node*3;
            rResult[index  ] = mpMasterElement->GetGeometry()[node].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index+1] = mpMasterElement->GetGeometry()[node].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index+2] = mpMasterElement->GetGeometry()[node].GetDof(DISPLACEMENT_Z).EquationId();
        }
        for( unsigned int node = 0; node < SlaveNN; node++ )
        {
            index = MasterNN*3 + node*3;
            rResult[index  ] = mpSlaveElement->GetGeometry()[node].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index+1] = mpSlaveElement->GetGeometry()[node].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index+2] = mpSlaveElement->GetGeometry()[node].GetDof(DISPLACEMENT_Z).EquationId();
        }

        index = MasterNN*3 + SlaveNN*3;
        rResult[index  ] = mpSlaveElement->GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_X).EquationId();
        rResult[index+1] = mpSlaveElement->GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Y).EquationId();
        rResult[index+2] = mpSlaveElement->GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Z).EquationId();
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * Setting up the DOF list for the current partners.
     * All conditions are assumed to be defined in 3D space with 3 DOFs per Node.
     * All DOF are given Master first, Slave second
     */
    //************************************************************************************
    //************************************************************************************
    void EmbeddedPointLagrangeTyingCondition::GetDofList( DofsVectorType& ConditionalDofList, const ProcessInfo& CurrentProcessInfo) const
    {
        ConditionalDofList.resize(0);
//        if(   ( mpMasterElement->GetValue(IS_INACTIVE) || !mpMasterElement->Is(ACTIVE) )
//           || ( mpSlaveElement->GetValue(IS_INACTIVE) || !mpSlaveElement->Is(ACTIVE) ) )
//        {
//            return;
//        }

        //determining size of DOF list
        //dimension of space
        unsigned int MasterNN = mpMasterElement->GetGeometry().size();
        unsigned int SlaveNN = mpSlaveElement->GetGeometry().size();

        for( unsigned int node = 0; node < MasterNN; ++node )
        {
            ConditionalDofList.push_back( mpMasterElement->GetGeometry()[node].pGetDof(DISPLACEMENT_X) );
            ConditionalDofList.push_back( mpMasterElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Y) );
            ConditionalDofList.push_back( mpMasterElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Z) );
        }

        for( unsigned int node = 0; node < SlaveNN; ++node )
        {
            ConditionalDofList.push_back( mpSlaveElement->GetGeometry()[node].pGetDof(DISPLACEMENT_X) );
            ConditionalDofList.push_back( mpSlaveElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Y) );
            ConditionalDofList.push_back( mpSlaveElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Z) );
        }

        // TODO be careful here, if we introduce LAGRANGE multiplier dof to node, that node must carry unique LAGRANGE dofs.
        // Because that L dofs are used to prescribe the contraint in this condition. Using L dofs at node may potentially
        // overlap with other master-slave couple
        ConditionalDofList.push_back( mpSlaveElement->GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_X) );
        ConditionalDofList.push_back( mpSlaveElement->GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Y) );
        ConditionalDofList.push_back( mpSlaveElement->GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Z) );
        //KRATOS_WATCH(ConditionalDofList[0]);
    }
} // Namespace Kratos

