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
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-03-25 08:20:14 $
//   Revision:            $Revision: 1.5 $
//
//
// System includes

// External includes

// Project includes
#include "custom_conditions/foundation_condition.h"
#include "structural_application.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    FoundationCondition::FoundationCondition( IndexType NewId,
                                  GeometryType::Pointer pGeometry)
    : Condition( NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    FoundationCondition::FoundationCondition( IndexType NewId, GeometryType::Pointer pGeometry,
                                  PropertiesType::Pointer pProperties
                                )
    : Condition( NewId, pGeometry, pProperties )
    {

    }

    FoundationCondition::FoundationCondition( IndexType NewId, GeometryType::Pointer pGeometry,
                                  	PropertiesType::Pointer pProperties,
					Element::Pointer& soilElement,
					Element::Pointer& foundationElement,
					Point<3>& rSoilLocalPoint,
					Point<3>& rFoundationLocalPoint)
	: Condition( NewId, pGeometry, pProperties )
	{
       		mSoilLocalPoint = rSoilLocalPoint;
		mFoundationLocalPoint = rFoundationLocalPoint;
		mpSoilElement = soilElement;
		mpFoundationElement = foundationElement;
		//Test for calculating coordinates at time step midpoint
        	//mSoilGlobalPoint = GlobalCoordinates(mpSoilElement, mSoilGlobalPoint, mSoilLocalPoint );
		//Test for calculating coordinates at time step midpoint
        	//mTipGlobalPoint = GlobalCoordinates(mpFoundationElement, mTipGlobalPoint, mFoundationLocalPoint );
//
	}

    //********************************************************
    //**** Operations ****************************************
    //********************************************************


    Condition::Pointer FoundationCondition::Create( IndexType NewId,
                                              NodesArrayType const& ThisNodes,
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new FoundationCondition(NewId, GetGeometry().Create(ThisNodes),
                                   pProperties));
    }

    Condition::Pointer FoundationCondition::Create( IndexType NewId,
                                              GeometryType::Pointer pGeom,
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new FoundationCondition(NewId, pGeom,
                                   pProperties));
    }

    /**
     * Destructor. Never to be called manually
     */
    FoundationCondition::~FoundationCondition()
    {
    }


    void FoundationCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void FoundationCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo)
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
    void FoundationCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
    //************************************************************************************    /
    /**
     * This function calculates all system contributions due to the contact problem
     * with regard to the current master and slave partners.
     * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node
     */
    void FoundationCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {

	KRATOS_TRY

	Vector uSoil = ZeroVector(3);
	for( unsigned int node = 0; node < mpSoilElement->GetGeometry().size(); node++ )
	{
		uSoil+=mpSoilElement->GetGeometry().ShapeFunctionValue(node,mSoilLocalPoint)*mpSoilElement->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);
	}
	//KRATOS_WATCH(uSoil);
//	KRATOS_WATCH(mSoilLocalPoint);
	Vector uFoundation = ZeroVector(3);
	for( unsigned int node = 0; node < mpFoundationElement->GetGeometry().size(); node++ )
	{
		uFoundation+=mpFoundationElement->GetGeometry().ShapeFunctionValue(node,mFoundationLocalPoint)*mpFoundationElement->GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);
	}
//	KRATOS_WATCH(uTip);
//	KRATOS_WATCH(mSoilLocalPoint);

	Vector relDisp = uSoil-uFoundation;
	//KRATOS_WATCH( relDisp );
        unsigned int soilNN = mpSoilElement->GetGeometry().size();
        unsigned int foundationNN = mpFoundationElement->GetGeometry().size();
        unsigned int dimension = 3;
        unsigned int MatSize = (soilNN+foundationNN)*dimension+3;
	///	KRATOS_WATCH(foundationNN);
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
        else return;

	//subtracting relDisp*penalty from soil nodes' reaction vector

	for( unsigned int node=0; node < soilNN; node++ )
	{
		rRightHandSideVector[3*node]   -= relDisp[0];///(mpSoilElement->GetGeometry().ShapeFunctionValue(node,mSoilLocalPoint));
		rRightHandSideVector[3*node+1] -= relDisp[1];///(mpSoilElement->GetGeometry().ShapeFunctionValue(node,mSoilLocalPoint));
		rRightHandSideVector[3*node+2] -= relDisp[2];////(mpSoilElement->GetGeometry().ShapeFunctionValue(node,mSoilLocalPoint));
	}

	for( unsigned int node=0; node < foundationNN; node++ )
	{

		rRightHandSideVector[3*soilNN+3*node]   += relDisp[0];////(mpFoundationElement->GetGeometry().ShapeFunctionValue(node, mFoundationLocalPoint));
		rRightHandSideVector[3*soilNN+3*node+1] += relDisp[1];////(mpFoundationElement->GetGeometry().ShapeFunctionValue(node, mFoundationLocalPoint));
		rRightHandSideVector[3*soilNN+3*node+2] += relDisp[2];////(mpFoundationElement->GetGeometry().ShapeFunctionValue(node, mFoundationLocalPoint));
	}

//	KRATOS_WATCH(Id())
//	KRATOS_WATCH(rRightHandSideVector)
/*	for ( unsigned int i = (soilNN+foundationNN)*3 ; i < (soilNN+foundationNN+1)*3; i++ )
        {
               rRightHandSideVector[i] -= relDisp[i] / 3;
        }*/


// std::cout<<"#################### THIS IS IN CALCULATE ALL ####################"<<std::endl;

        Vector soil_local_shape = ZeroVector(soilNN);
    	Vector side_local_shape = ZeroVector(foundationNN);

//         interface_local_shape= mpInterfaceElement->GetShapeFunctionValues(mSurfaceIntegrationPointIndex, mLeftOrRight);

//		double soil_incremental_area= mpSoilElement->GetIncrementalArea(mSurfaceIntegrationPointIndex, mLeftOrRight);



        for( IndexType PointNumber = 0; PointNumber < mpFoundationElement->GetGeometry().size(); PointNumber++ )
        {
            side_local_shape[PointNumber] = mpFoundationElement->GetGeometry().ShapeFunctionValue( PointNumber, mFoundationLocalPoint);
	   // KRATOS_WATCH (side_local_shape);
        }
        for ( unsigned int i = 0; i < 3; i++ )
        {
                for ( unsigned int foundation_node = 0; foundation_node < foundationNN; foundation_node++ )
                {
                    for (  unsigned int soil_node = 0; soil_node < soilNN; soil_node++ )

                    {	//KRATOS_WATCH (mpSoilElement->GetGeometry().ShapeFunctionValue(soil_node,mSoilLocalPoint))
                        //rLeftHandSideMatrix( soil_node * 3 + i, soilNN * 3 + foundationNN * 3 + i ) -= 4.0*(mpSoilElement->GetGeometry().ShapeFunctionValue(soil_node,mSoilLocalPoint));
                       // rLeftHandSideMatrix( soilNN * 3 + foundationNN * 3 + i, soil_node * 3 + i ) -= 4.0*(mpSoilElement->GetGeometry().ShapeFunctionValue(soil_node,mSoilLocalPoint));
                       // rLeftHandSideMatrix( soilNN * 3 + foundation_node * 3 + i, soilNN * 3 + foundationNN * 3 + i ) += 0.5*(mpFoundationElement->GetGeometry().ShapeFunctionValue(foundation_node, mFoundationLocalPoint));
                       // rLeftHandSideMatrix( soilNN * 3 + foundationNN * 3 + i, soilNN * 3 + foundation_node * 3 + i ) += 0.5*(mpFoundationElement->GetGeometry().ShapeFunctionValue(foundation_node, mFoundationLocalPoint));
                        rLeftHandSideMatrix( soil_node * 3 + i, soilNN * 3 + foundationNN * 3 + i ) -= 1.0*(mpSoilElement->GetGeometry().ShapeFunctionValue(soil_node,mSoilLocalPoint));
                        rLeftHandSideMatrix( soilNN * 3 + foundationNN * 3 + i, soil_node * 3 + i ) -= 1.0*(mpSoilElement->GetGeometry().ShapeFunctionValue(soil_node,mSoilLocalPoint));
                        rLeftHandSideMatrix( soilNN * 3 + foundation_node * 3 + i, soilNN * 3 + foundationNN * 3 + i ) += 1.0*(mpFoundationElement->GetGeometry().ShapeFunctionValue(foundation_node, mFoundationLocalPoint));
                        rLeftHandSideMatrix( soilNN * 3 + foundationNN * 3 + i, soilNN * 3 + foundation_node * 3 + i ) += 1.0*(mpFoundationElement->GetGeometry().ShapeFunctionValue(foundation_node, mFoundationLocalPoint));
                     }
                 }
	}
//	KRATOS_WATCH (rLeftHandSideMatrix);

 /*       for( unsigned int foundation_node = 0; foundation_node < foundationNN; foundation_node++ )
        {
            for( unsigned int i=0; i<3; i++ )
            {
                for( unsigned int soil_node = 0; soil_node < soilNN; soil_node++ )
                {
                    //set connection soil_node/lagrange -> -1
                    rLeftHandSideMatrix(soilNN*3+foundation_node*3+3+i,3*soil_node+i) -=
                            1.0*(mpSoilElement->GetGeometry().ShapeFunctionValue(soil_node,mSoilLocalPoint));
                    rLeftHandSideMatrix(3*soil_node+i,soilNN*3+foundation_node*3+3+i) -=
                            1.0*(mpSoilElement->GetGeometry().ShapeFunctionValue(soil_node,mSoilLocalPoint));
                }
                //set connection foundation_node/lagrange -> +1
                rLeftHandSideMatrix(3*soilNN+3*foundation_node+3+i,3*soilNN+3*foundation_node+i) +=1*(mpFoundationElement->GetGeometry().ShapeFunctionValue(foundation_node, mFoundationLocalPoint));
                rLeftHandSideMatrix(3*soilNN+3*foundation_node+i,3*soilNN+3*foundation_node+3+i) +=1*(mpFoundationElement->GetGeometry().ShapeFunctionValue(foundation_node, mFoundationLocalPoint));
            }
        }*/
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
	void FoundationCondition::EquationIdVector( EquationIdVectorType& rResult,
                                          ProcessInfo& CurrentProcessInfo)
	{
        //determining size of DOF list
        //dimension of space
		unsigned int soilNN = mpSoilElement->GetGeometry().size();
	        unsigned int foundationNN = mpFoundationElement->GetGeometry().size();
		unsigned int ndofs = 3*(soilNN+foundationNN)+3;
	        unsigned int index;
		rResult.resize(ndofs,false);
 	        for( unsigned int node=0; node<soilNN; node++ )
 	        {
			index = node*3;
			rResult[index]   = mpSoilElement->GetGeometry()[node].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index+1] = mpSoilElement->GetGeometry()[node].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index+2] = mpSoilElement->GetGeometry()[node].GetDof(DISPLACEMENT_Z).EquationId();
 	        }
 	        for( unsigned int node=0; node<foundationNN; node++ )
 	        {
			index = soilNN*3+node*3;
			rResult[index]   = mpFoundationElement->GetGeometry()[node].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index+1] = mpFoundationElement->GetGeometry()[node].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index+2] = mpFoundationElement->GetGeometry()[node].GetDof(DISPLACEMENT_Z).EquationId();
//			rResult[index+3] = mpFoundationElement->GetGeometry()[node].GetDof(ROTATION_X).EquationId();
//			rResult[index+4] = mpFoundationElement->GetGeometry()[node].GetDof(ROTATION_Y).EquationId();
//			rResult[index+5] = mpFoundationElement->GetGeometry()[node].GetDof(ROTATION_Z).EquationId();
 	        }
	        index = soilNN*3+foundationNN*3;
		rResult[index] = mpFoundationElement->GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_X).EquationId();
		rResult[index+1] = mpFoundationElement->GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Y).EquationId();
		rResult[index+2] = mpFoundationElement->GetGeometry()[0].GetDof(LAGRANGE_DISPLACEMENT_Z).EquationId();
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
    void FoundationCondition::GetDofList( DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
    {
        //determining size of DOF list
        //dimension of space
        unsigned int soilNN = mpSoilElement->GetGeometry().size();
        unsigned int foundationNN = mpFoundationElement->GetGeometry().size();
        unsigned int index;
	unsigned int ndofs = 3*(soilNN+foundationNN)+3;
//	KRATOS_WATCH(ndofs);
        ConditionalDofList.resize(ndofs);
	for( unsigned int node=0; node<soilNN; node++ )
        {
            index = node*3;
            ConditionalDofList[index]   = mpSoilElement->GetGeometry()[node].pGetDof(DISPLACEMENT_X);
            ConditionalDofList[index+1] = mpSoilElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Y);
            ConditionalDofList[index+2] = mpSoilElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Z);
        }
        for( unsigned int node=0; node<foundationNN; node++ )
        {
            index = soilNN*3+node*3;
            ConditionalDofList[index]   = mpFoundationElement->GetGeometry()[node].pGetDof(DISPLACEMENT_X);
            ConditionalDofList[index+1] = mpFoundationElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Y);
            ConditionalDofList[index+2] = mpFoundationElement->GetGeometry()[node].pGetDof(DISPLACEMENT_Z);
 //           ConditionalDofList[index+3] = mpFoundationElement->GetGeometry()[node].pGetDof(ROTATION_X);
   //         ConditionalDofList[index+4] = mpFoundationElement->GetGeometry()[node].pGetDof(ROTATION_Y);
     //       ConditionalDofList[index+5] = mpFoundationElement->GetGeometry()[node].pGetDof(ROTATION_Z);
        }
	index = soilNN*3+foundationNN*3;
	ConditionalDofList[index] = mpFoundationElement->GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_X);
	ConditionalDofList[index+1] = mpFoundationElement->GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Y);
	ConditionalDofList[index+2] = mpFoundationElement->GetGeometry()[0].pGetDof(LAGRANGE_DISPLACEMENT_Z);
	//KRATOS_WATCH(ConditionalDofList[0]);
    }
} // Namespace Kratos
