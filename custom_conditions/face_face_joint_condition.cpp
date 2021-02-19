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
*   Date:                $Date: 6 Dec 2017 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "geometries/line_3d_2.h"
#include "utilities/math_utils.h"
#include "custom_conditions/face_face_joint_condition.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************

FaceFaceJointCondition::FaceFaceJointCondition( IndexType NewId, GeometryType::Pointer geom1, GeometryType::Pointer geom2,
        PropertiesType::Pointer pProperties )
    : Condition( NewId, GeometryType::Pointer( new Line3D2<Node<3> >( (*geom1)(0), (*geom2)(0) ) ), pProperties )
    , mpGeom1(geom1), mpGeom2(geom2)
{
    if (mpGeom1->GetGeometryType() != mpGeom2->GetGeometryType())
        KRATOS_THROW_ERROR(std::logic_error, "The geometry type is incompatible", "")
    mThisIntegrationMethod = mpGeom1->GetDefaultIntegrationMethod();//default method
}

FaceFaceJointCondition::FaceFaceJointCondition( IndexType NewId, Condition::Pointer cond1, Condition::Pointer cond2,
        PropertiesType::Pointer pProperties )
    : Condition( NewId, GeometryType::Pointer( new Line3D2<Node<3> >( (*cond1->pGetGeometry())(0), (*cond2->pGetGeometry())(0) ) ), pProperties )
    , mpGeom1(cond1->pGetGeometry()), mpGeom2(cond2->pGetGeometry())
{
    if (mpGeom1->GetGeometryType() != mpGeom2->GetGeometryType())
        KRATOS_THROW_ERROR(std::logic_error, "The geometry type is incompatible", "")
    mThisIntegrationMethod = mpGeom1->GetDefaultIntegrationMethod();//default method
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************

Condition::Pointer FaceFaceJointCondition::Create( IndexType NewId,
        GeometryType::Pointer geom1, GeometryType::Pointer geom2,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new FaceFaceJointCondition(NewId, geom1, geom2,
                                pProperties));
}

Condition::Pointer FaceFaceJointCondition::Create( IndexType NewId,
        Condition::Pointer cond1, Condition::Pointer cond2,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new FaceFaceJointCondition(NewId, cond1, cond2,
                                pProperties));
}



/**
 * Destructor. Never to be called manually
 */
FaceFaceJointCondition::~FaceFaceJointCondition()
{
}

void FaceFaceJointCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    // integration rule
    if(this->Has( INTEGRATION_ORDER ))
    {
        if(this->GetValue(INTEGRATION_ORDER) == 1)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 4)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "KinematicLinear element does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
    }
    else if(GetProperties().Has( INTEGRATION_ORDER ))
    {
        if(GetProperties()[INTEGRATION_ORDER] == 1)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 4)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "KinematicLinear element does not support for integration points", GetProperties()[INTEGRATION_ORDER])
    }
    else
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method
}

//************************************************************************************
//************************************************************************************

/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void FaceFaceJointCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
void FaceFaceJointCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
void FaceFaceJointCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag)
{
    unsigned int dim = 3;
    unsigned int MatSize = dim * (mpGeom1->size() + mpGeom2->size());

    //resize the LHS=StiffnessMatrix if its size is not correct
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize || rLeftHandSideMatrix.size2() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    // resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        //resize the RHS=force vector if its size is not correct
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    mpGeom1->Initialize(mThisIntegrationMethod);
//    mpGeom2->Initialize(mThisIntegrationMethod); // there is potential data race condition here. When different Bezier geometry may access the Bezier data at simultaneously. TODO: need check further
    #endif

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = mpGeom1->IntegrationPoints(mThisIntegrationMethod);
//    const GeometryType::ShapeFunctionsGradientsType& DN_De = mpGeom1->ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    Vector ShapeFunctionValues1;
    Vector ShapeFunctionValues2;

    //calculate the jacobian in reference configuration
    Matrix DeltaPosition(mpGeom1->size(), 3);
    for ( unsigned int i = 0; i < mpGeom1->size(); ++i )
    {
        noalias( row( DeltaPosition, i ) ) = (*mpGeom1)[i].GetSolutionStepValue(DISPLACEMENT);
    }

    GeometryType::JacobiansType J0;
    J0 = mpGeom1->Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

    const Matrix& joint_stiff = GetProperties()[ JOINT_STIFFNESS ]; // must be of size 6x6
    Matrix aux(6, 6);
    Vector displacements(6), forces(6);
    unsigned int offset_row, offset_col;

//    KRATOS_WATCH(mThisIntegrationMethod)
//    KRATOS_WATCH(integration_points.size())

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        double dA = sqrt(MathUtils<double>::Det(Matrix(prod(trans(J0[PointNumber]), J0[PointNumber]))));

        double IntegrationWeight = integration_points[PointNumber].Weight();

        noalias(aux) = joint_stiff * dA * IntegrationWeight;

        if ( CalculateResidualVectorFlag == true || CalculateStiffnessMatrixFlag == true )
        {
            ShapeFunctionValues1 = mpGeom1->ShapeFunctionsValues(ShapeFunctionValues1, integration_points[PointNumber]);
            ShapeFunctionValues2 = mpGeom2->ShapeFunctionsValues(ShapeFunctionValues2, integration_points[PointNumber]);

            noalias(displacements) = ZeroVector(6);
            for (unsigned int i = 0; i < mpGeom1->size(); ++i)
            {
                displacements[0] += ShapeFunctionValues1(i) * (*mpGeom1)[i].GetSolutionStepValue(DISPLACEMENT_X);
                displacements[1] += ShapeFunctionValues1(i) * (*mpGeom1)[i].GetSolutionStepValue(DISPLACEMENT_Y);
                displacements[2] += ShapeFunctionValues1(i) * (*mpGeom1)[i].GetSolutionStepValue(DISPLACEMENT_Z);
            }
            for (unsigned int i = 0; i < mpGeom2->size(); ++i)
            {
                displacements[3] += ShapeFunctionValues2(i) * (*mpGeom2)[i].GetSolutionStepValue(DISPLACEMENT_X);
                displacements[4] += ShapeFunctionValues2(i) * (*mpGeom2)[i].GetSolutionStepValue(DISPLACEMENT_Y);
                displacements[5] += ShapeFunctionValues2(i) * (*mpGeom2)[i].GetSolutionStepValue(DISPLACEMENT_Z);
            }

            noalias(forces) = prod(aux, displacements);

            offset_row = 0;
            for (unsigned int i = 0; i < mpGeom1->size(); ++i)
            {
                subrange(rRightHandSideVector, offset_row + 3*i, offset_row + 3*i+3) -= ShapeFunctionValues1(i) * subrange(forces, 0, 3);
            }

            if ( CalculateStiffnessMatrixFlag == true )
            {
                for (unsigned int i = 0; i < mpGeom1->size(); ++i)
                {
                    offset_col = 0;
                    for (unsigned int j = 0; j < mpGeom1->size(); ++j)
                        subrange(rLeftHandSideMatrix, offset_row + 3*i, offset_row + 3*i+3,
                            offset_col + 3*j, offset_col + 3*j+3) += ShapeFunctionValues1(i) * ShapeFunctionValues1(j) * subrange(aux, 0, 3, 0, 3);
                    offset_col = 3*mpGeom1->size();
                    for (unsigned int j = 0; j < mpGeom2->size(); ++j)
                        subrange(rLeftHandSideMatrix, offset_row + 3*i, offset_row + 3*i+3,
                            offset_col + 3*j, offset_col + 3*j+3) += ShapeFunctionValues1(i) * ShapeFunctionValues2(j) * subrange(aux, 0, 3, 3, 6);
                }
            }

            offset_row = 3*mpGeom1->size();
            for (unsigned int i = 0; i < mpGeom2->size(); ++i)
            {
                subrange(rRightHandSideVector, offset_row + 3*i, offset_row + 3*i+3) -= ShapeFunctionValues2(i) * subrange(forces, 3, 6);
            }

            if ( CalculateStiffnessMatrixFlag == true )
            {
                for (unsigned int i = 0; i < mpGeom2->size(); ++i)
                {
                    offset_col = 0;
                    for (unsigned int j = 0; j < mpGeom1->size(); ++j)
                        subrange(rLeftHandSideMatrix, offset_row + 3*i, offset_row + 3*i+3,
                            offset_col + 3*j, offset_col + 3*j+3) += ShapeFunctionValues2(i) * ShapeFunctionValues1(j) * subrange(aux, 3, 6, 0, 3);
                    offset_col = 3*mpGeom1->size();
                    for (unsigned int j = 0; j < mpGeom2->size(); ++j)
                        subrange(rLeftHandSideMatrix, offset_row + 3*i, offset_row + 3*i+3,
                            offset_col + 3*j, offset_col + 3*j+3) += ShapeFunctionValues2(i) * ShapeFunctionValues2(j) * subrange(aux, 3, 6, 3, 6);
                }
            }
        }
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    //clean the geometry
    mpGeom1->Clean();
    mpGeom2->Clean();
    #endif
}

//************************************************************************************
//************************************************************************************

void FaceFaceJointCondition::EquationIdVector( EquationIdVectorType& rResult,
        const ProcessInfo& CurrentProcessInfo ) const
{
    DofsVectorType ConditionalDofList;
    GetDofList(ConditionalDofList, CurrentProcessInfo);
    if (rResult.size() != ConditionalDofList.size())
        rResult.resize(ConditionalDofList.size());
    for ( unsigned int i = 0; i < ConditionalDofList.size(); ++i)
        rResult[i] = ConditionalDofList[i]->EquationId();
}

//************************************************************************************
//************************************************************************************
void FaceFaceJointCondition::GetDofList( DofsVectorType& ConditionalDofList,
                                        const ProcessInfo& CurrentProcessInfo) const
{
    //determining size of DOF list
    //dimension of space
    unsigned int dim = 3;
    ConditionalDofList.resize(dim*(mpGeom1->size() + mpGeom2->size()));
    std::size_t cnt = 0;
    for( unsigned int i = 0; i < mpGeom1->size(); ++i )
    {
        ConditionalDofList[cnt++] = (*mpGeom1)[i].pGetDof(DISPLACEMENT_X);
        ConditionalDofList[cnt++] = (*mpGeom1)[i].pGetDof(DISPLACEMENT_Y);
        ConditionalDofList[cnt++] = (*mpGeom1)[i].pGetDof(DISPLACEMENT_Z);
    }
    for( unsigned int i = 0; i < mpGeom2->size(); ++i )
    {
        ConditionalDofList[cnt++] = (*mpGeom2)[i].pGetDof(DISPLACEMENT_X);
        ConditionalDofList[cnt++] = (*mpGeom2)[i].pGetDof(DISPLACEMENT_Y);
        ConditionalDofList[cnt++] = (*mpGeom2)[i].pGetDof(DISPLACEMENT_Z);
    }
}

} // Namespace Kratos
