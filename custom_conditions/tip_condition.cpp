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
//   Last Modified by:    $Author: jelena $
//   Date:                $Date: 2015-03-25$
//   Revision:            $Revision: 1.2 $
//
//
// System includes

// External includes

// Project includes
#include "custom_utilities/sd_math_utils.h"
#include "custom_conditions/tip_condition.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    TipCondition::TipCondition( IndexType NewId,
                                  GeometryType::Pointer pGeometry)
    : Condition( NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    TipCondition::TipCondition( IndexType NewId, GeometryType::Pointer pGeometry,
                                  PropertiesType::Pointer pProperties
                                )
    : Condition( NewId, pGeometry, pProperties )
    {

    }

    TipCondition::TipCondition( IndexType NewId, GeometryType::Pointer pGeometry,
          PropertiesType::Pointer pProperties,
          Element::Pointer tip_soilElement, Element::Pointer tipElement,
          PointType& rTipSoilLocalPoint, PointType& rTipLocalPoint)
    : Condition( NewId, pGeometry, pProperties )
    {
        mTipSoilLocalPoint = rTipSoilLocalPoint;
        mTipLocalPoint = rTipLocalPoint;
        mpTipSoilElement = tip_soilElement;
        mpTipElement = tipElement;
    }

    //********************************************************
    //**** Operations ****************************************
    //********************************************************


    Condition::Pointer TipCondition::Create( IndexType NewId,
                                              NodesArrayType const& ThisNodes,
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new TipCondition(NewId, GetGeometry().Create(ThisNodes),
                                   pProperties));
    }
    /**
     * Destructor. Never to be called manually
     */
    TipCondition::~TipCondition()
    {
    }


    void TipCondition::Initialize(const ProcessInfo& CurrentProcessInfo)
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
    void TipCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
    void TipCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
    void TipCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        // double penalty = GetProperties()[INITIAL_PENALTY];

        Vector uTipSoil = ZeroVector(3);
        for( unsigned int node = 0; node < mpTipSoilElement->GetGeometry().size(); node++ )
        {
            uTipSoil+=mpTipSoilElement->GetGeometry().ShapeFunctionValue(node,mTipSoilLocalPoint)*mpTipSoilElement->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);
        }
        //KRATOS_WATCH(uTipSoil);
    //  KRATOS_WATCH(mTipSoilLocalPoint);
        Vector uTip = ZeroVector(3);
        for( unsigned int node = 0; node < mpTipElement->GetGeometry().size(); node++ )
        {
            uTip+=mpTipElement->GetGeometry().ShapeFunctionValue(node,mTipLocalPoint)*mpTipElement->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);
        }
    //  KRATOS_WATCH(uTip);
    //  KRATOS_WATCH(mTipSoilLocalPoint);

        Vector relDisp = uTipSoil-uTip;
        KRATOS_WATCH( relDisp );
        unsigned int tip_soilNN = mpTipSoilElement->GetGeometry().size();
        unsigned int tipNN = mpTipElement->GetGeometry().size();
        unsigned int dimension = 3;
        unsigned int MatSize = (tip_soilNN+tipNN)*dimension+3;
        KRATOS_WATCH(MatSize);
        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if(rRightHandSideVector.size() != MatSize)
                rRightHandSideVector.resize(MatSize);
            noalias(rRightHandSideVector) = ZeroVector(MatSize); //resetting RHS    */
        }
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if(rLeftHandSideMatrix.size1() != MatSize || rLeftHandSideMatrix.size2() != MatSize)
                rLeftHandSideMatrix.resize(MatSize,MatSize);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize);
        }

        if (!CalculateResidualVectorFlag && !CalculateStiffnessMatrixFlag)
            return;

        //subtracting relDisp*penalty from soil nodes' reaction vector
        for( unsigned int node=0; node < tip_soilNN; node++ )
        {
            rRightHandSideVector[3*node]   -= relDisp[0];
            rRightHandSideVector[3*node+1] -= relDisp[1];
            rRightHandSideVector[3*node+2] -= relDisp[2];
        }
        //adding relDisp*penalty to pile's reaction vector
        for( unsigned int node=0; node < tipNN; node++ )
        {
            rRightHandSideVector[3*tip_soilNN+3*node]   += relDisp[0];
            rRightHandSideVector[3*tip_soilNN+3*node+1] += relDisp[1];
            rRightHandSideVector[3*tip_soilNN+3*node+2] += relDisp[2];

        }
        //KRATOS_WATCH(rRightHandSideVector)


// std::cout<<"#################### THIS IS IN CALCULATE ALL ####################"<<std::endl;

        Vector tip_soil_local_shape = ZeroVector(tip_soilNN);
        Vector side_local_shape = ZeroVector(tipNN);

        for( IndexType PointNumber = 0; PointNumber < mpTipElement->GetGeometry().size(); PointNumber++ )
        {
            side_local_shape[PointNumber] = mpTipElement->GetGeometry().ShapeFunctionValue( PointNumber, mTipLocalPoint);
        }

        for( unsigned int tip_node = 0; tip_node < tipNN; tip_node++ )
        {
            for( unsigned int i=0; i<3; i++ )
            {
                for( unsigned int tip_soil_node = 0; tip_soil_node < tip_soilNN; tip_soil_node++ )
                {
                    //set connection tip_soil_node/lagrange -> -1
                    rLeftHandSideMatrix(tip_soilNN*3+tip_node*3+3+i,3*tip_soil_node+i) -=
                            1.0*(mpTipSoilElement->GetGeometry().ShapeFunctionValue(tip_soil_node,mTipSoilLocalPoint));
                    rLeftHandSideMatrix(3*tip_soil_node+i,tip_soilNN*3+tip_node*3+3+i) -=
                            1.0*(mpTipSoilElement->GetGeometry().ShapeFunctionValue(tip_soil_node,mTipSoilLocalPoint));
                }
                //set connection tip_node/lagrange -> +1
                rLeftHandSideMatrix(3*tip_soilNN+3*tip_node+3+i,3*tip_soilNN+3*tip_node+i) +=
                        1*(mpTipElement->GetGeometry().ShapeFunctionValue(tip_node,mTipLocalPoint));////
                rLeftHandSideMatrix(3*tip_soilNN+3*tip_node+i,3*tip_soilNN+3*tip_node+3+i) +=
                        1*(mpTipElement->GetGeometry().ShapeFunctionValue(tip_node,mTipLocalPoint));
            }
        }
//        KRATOS_WATCH(rLeftHandSideMatrix);
        KRATOS_CATCH("")
    }


//************************************************************************************
//************************************************************************************
/**
* Setting up the EquationIdVector for the current partners.
* All conditions are assumed to be defined in 3D space with 3 DOFs per node.
* All Equation IDs are given Master first, Slave second
*/
    void TipCondition::EquationIdVector( EquationIdVectorType& rResult,
                                          const ProcessInfo& CurrentProcessInfo) const
    {
        //determining size of DOF list
        //dimension of space
        unsigned int tip_soilNN = mpTipSoilElement->GetGeometry().size();
        unsigned int tipNN = mpTipElement->GetGeometry().size();
        unsigned int ndofs = 3*(tip_soilNN+tipNN)+3;
        unsigned int index;

        rResult.resize(ndofs,false);

        for( unsigned int node=0; node<tip_soilNN; node++ )
        {
            index = node*3;
            rResult[index]   = mpTipSoilElement->GetGeometry()[node].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index+1] = mpTipSoilElement->GetGeometry()[node].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index+2] = mpTipSoilElement->GetGeometry()[node].GetDof(DISPLACEMENT_Z).EquationId();
        }
        for( unsigned int node=0; node<tipNN; node++ )
        {
            index = tip_soilNN*3+node*3;
            rResult[index]   = mpTipElement->GetGeometry()[node].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index+1] = mpTipElement->GetGeometry()[node].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index+2] = mpTipElement->GetGeometry()[node].GetDof(DISPLACEMENT_Z).EquationId();
        }

        index = tip_soilNN*3+tipNN*3;
        rResult[index] = mpTipElement->GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_X).EquationId();
        rResult[index+1] = mpTipElement->GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Y).EquationId();
        rResult[index+2] = mpTipElement->GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Z).EquationId();
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
    void TipCondition::GetDofList( DofsVectorType& ConditionalDofList, const ProcessInfo& CurrentProcessInfo) const
    {
        //determining size of DOF list
        //dimension of space
        unsigned int tip_soilNN = mpTipSoilElement->GetGeometry().size();
        unsigned int tipNN = mpTipElement->GetGeometry().size();
        unsigned int index;
    unsigned int ndofs = 3*(tip_soilNN+tipNN)+3;
//  KRATOS_WATCH(ndofs);
        ConditionalDofList.resize(ndofs);
    for( unsigned int node=0; node<tip_soilNN; node++ )
        {
            index = node*3;
            ConditionalDofList[index]   = mpTipSoilElement->GetGeometry()[node].pGetDof(DISPLACEMENT_X);
            ConditionalDofList[index+1] = mpTipSoilElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Y);
            ConditionalDofList[index+2] = mpTipSoilElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Z);
        }
        for( unsigned int node=0; node<tipNN; node++ )
        {
            index = tip_soilNN*3+node*3;
            ConditionalDofList[index]   = mpTipElement->GetGeometry()[node].pGetDof(DISPLACEMENT_X);
            ConditionalDofList[index+1] = mpTipElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Y);
            ConditionalDofList[index+2] = mpTipElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Z);
        }
    index = tip_soilNN*3+tipNN*3;
    ConditionalDofList[index] = mpTipElement->GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_X);
    ConditionalDofList[index+1] = mpTipElement->GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Y);
    ConditionalDofList[index+2] = mpTipElement->GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Z);
    //KRATOS_WATCH(ConditionalDofList[0]);
    }
} // Namespace Kratos
