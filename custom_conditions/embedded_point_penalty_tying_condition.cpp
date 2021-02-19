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
//   Date:                $Date: 22 Jan 2017$
//   Revision:            $Revision: 1.2 $
//
//
// System includes

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "custom_conditions/embedded_point_penalty_tying_condition.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    EmbeddedPointPenaltyTyingCondition::EmbeddedPointPenaltyTyingCondition()
    {
    }

    EmbeddedPointPenaltyTyingCondition::EmbeddedPointPenaltyTyingCondition( IndexType NewId,
                                  GeometryType::Pointer pGeometry)
    : Condition( NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    EmbeddedPointPenaltyTyingCondition::EmbeddedPointPenaltyTyingCondition( IndexType NewId,
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


    Condition::Pointer EmbeddedPointPenaltyTyingCondition::Create( IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       Element::Pointer pMasterElement,
                                       Element::Pointer pSlaveElement,
                                       PointType& rMasterLocalPoint,
                                       PointType& rSlaveLocalPoint ) const
    {
        return Condition::Pointer( new EmbeddedPointPenaltyTyingCondition(NewId, pGeometry, pMasterElement, pSlaveElement, rMasterLocalPoint, rSlaveLocalPoint) );
    }

    /**
     * Destructor. Never to be called manually
     */
    EmbeddedPointPenaltyTyingCondition::~EmbeddedPointPenaltyTyingCondition()
    {
    }


    void EmbeddedPointPenaltyTyingCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
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
    void EmbeddedPointPenaltyTyingCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
    void EmbeddedPointPenaltyTyingCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
    void EmbeddedPointPenaltyTyingCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        unsigned int MasterNN = mpMasterElement->GetGeometry().size();
        unsigned int SlaveNN = mpSlaveElement->GetGeometry().size();
        unsigned int MatSize = (MasterNN+SlaveNN)*3;

        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if(rRightHandSideVector.size() != MatSize)
                rRightHandSideVector.resize(MatSize, false);
            noalias(rRightHandSideVector) = ZeroVector(MatSize);
        }
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if(rLeftHandSideMatrix.size1() != MatSize || rLeftHandSideMatrix.size2() != MatSize)
                rLeftHandSideMatrix.resize(MatSize, MatSize, false);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize);
        }

        if(   ( mpMasterElement->GetValue(IS_INACTIVE) || !mpMasterElement->Is(ACTIVE) )
           || ( mpSlaveElement->GetValue(IS_INACTIVE) || !mpSlaveElement->Is(ACTIVE) ) )
        {
            return;
        }

        if( !CalculateStiffnessMatrixFlag && !CalculateResidualVectorFlag )
            return;

        if( !this->Has(INITIAL_PENALTY) )
            KRATOS_THROW_ERROR(std::logic_error, "INITIAL_PENALTY is not set for embbedded tying link", Id())
        double Penalty = this->GetValue(INITIAL_PENALTY);

        Vector uMaster = ZeroVector(3);
        Vector ShapeFunctionValuesOnMaster;
        ShapeFunctionValuesOnMaster = mpMasterElement->GetGeometry().ShapeFunctionsValues(ShapeFunctionValuesOnMaster, mMasterLocalPoint);
        for( unsigned int node = 0; node < mpMasterElement->GetGeometry().size(); ++node )
        {
            uMaster += ShapeFunctionValuesOnMaster[node] * mpMasterElement->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);
        }

        Vector uSlave = ZeroVector(3);
        Vector ShapeFunctionValuesOnSlave;
        ShapeFunctionValuesOnSlave = mpSlaveElement->GetGeometry().ShapeFunctionsValues(ShapeFunctionValuesOnSlave, mSlaveLocalPoint);
        for( unsigned int node = 0; node < mpSlaveElement->GetGeometry().size(); ++node )
        {
            uSlave += ShapeFunctionValuesOnSlave[node] * mpSlaveElement->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);
        }

        Vector relDisp = uSlave - uMaster;

        for( unsigned int node = 0; node < MasterNN; ++node )
        {
            for( unsigned int dim = 0; dim < 3; ++dim )
                rRightHandSideVector[3*node + dim] -= Penalty * relDisp[dim] * ShapeFunctionValuesOnMaster[node];
        }
        for( unsigned int node = 0; node < SlaveNN; ++node )
        {
            for( unsigned int dim = 0; dim < 3; ++dim )
                rRightHandSideVector[3*MasterNN + 3*node + dim] += Penalty * relDisp[dim] * ShapeFunctionValuesOnSlave[node];
        }

        for( unsigned int prim = 0; prim < MasterNN; ++prim )
        {
            for( unsigned int sec = 0; sec < MasterNN; ++sec )
            {
                for( unsigned int pdim = 0; pdim < 3; ++pdim )
                {
                    for( unsigned int sdim = 0; sdim < 3; ++sdim )
                    {
                        rLeftHandSideMatrix(3*prim + pdim, 3*sec+sdim)
                            += Penalty*ShapeFunctionValuesOnMaster[prim]*ShapeFunctionValuesOnMaster[sec];
                    }
                }
            }

            for( unsigned int sec = 0; sec < SlaveNN; ++sec )
            {
                for( unsigned int pdim = 0; pdim < 3; ++pdim )
                {
                    for( unsigned int sdim = 0; sdim < 3; ++sdim )
                    {
                        rLeftHandSideMatrix(3*prim + pdim, 3*MasterNN + 3*sec+sdim)
                            += Penalty*ShapeFunctionValuesOnMaster[prim]*ShapeFunctionValuesOnSlave[sec];
                    }
                }
            }
        }

        for( unsigned int prim = 0; prim < SlaveNN; ++prim )
        {
            for( unsigned int sec = 0; sec < MasterNN; ++sec )
            {
                for( unsigned int pdim = 0; pdim < 3; ++pdim )
                {
                    for( unsigned int sdim = 0; sdim < 3; ++sdim )
                    {
                        rLeftHandSideMatrix(3*MasterNN + 3*prim + pdim, 3*sec+sdim)
                            -= Penalty*ShapeFunctionValuesOnSlave[prim]*ShapeFunctionValuesOnMaster[sec];
                    }
                }
            }

            for( unsigned int sec = 0; sec < SlaveNN; ++sec )
            {
                for( unsigned int pdim = 0; pdim < 3; ++pdim )
                {
                    for( unsigned int sdim = 0; sdim < 3; ++sdim )
                    {
                        rLeftHandSideMatrix(3*MasterNN + 3*prim + pdim, 3*MasterNN + 3*sec+sdim)
                            += Penalty*ShapeFunctionValuesOnSlave[prim]*ShapeFunctionValuesOnSlave[sec];
                    }
                }
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
    void EmbeddedPointPenaltyTyingCondition::EquationIdVector( EquationIdVectorType& rResult,
                                          const ProcessInfo& CurrentProcessInfo) const
    {
        //determining size of DOF list
        //dimension of space
        unsigned int MasterNN = mpMasterElement->GetGeometry().size();
        unsigned int SlaveNN = mpSlaveElement->GetGeometry().size();
        unsigned int ndofs = 3*(MasterNN+SlaveNN);
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
    void EmbeddedPointPenaltyTyingCondition::GetDofList( DofsVectorType& ConditionalDofList, const ProcessInfo& CurrentProcessInfo) const
    {
        ConditionalDofList.resize(0);

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
    }
} // Namespace Kratos

