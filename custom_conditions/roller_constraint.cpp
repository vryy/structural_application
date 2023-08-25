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
*   Date:                $Date: 18 Aug 2023 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "utilities/math_utils.h"
#include "geometries/line_2d_3.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"
#include "custom_conditions/roller_constraint.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
RollerConstraint::RollerConstraint( IndexType NewId,
        GeometryType::Pointer pGeometry) :
    Condition( NewId, pGeometry )
{
}

RollerConstraint::RollerConstraint( IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) :
    Condition( NewId, pGeometry, pProperties )
{
}

RollerConstraint::RollerConstraint( IndexType NewId, Node<3>::Pointer const& node1, Node<3>::Pointer const& node2, Node<3>::Pointer const& node3,
        PropertiesType::Pointer pProperties )
    : Condition( NewId, GeometryType::Pointer( new Line2D3<Node<3> >( node1, node2, node3 ) ), pProperties )
{
}

RollerConstraint::RollerConstraint( IndexType NewId, NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties )
    : Condition( NewId, GeometryType::Pointer( new Line2D3<Node<3> >( ThisNodes ) ), pProperties )
{
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************

Condition::Pointer RollerConstraint::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new RollerConstraint(NewId, GetGeometry().Create(ThisNodes),
                                pProperties));
}

Condition::Pointer RollerConstraint::Create( IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new RollerConstraint(NewId, pGeom,
                                pProperties));
}

/**
 * Destructor. Never to be called manually
 */
RollerConstraint::~RollerConstraint()
{
}

//************************************************************************************
//************************************************************************************

void RollerConstraint::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY //EXCEPTION HANDLING (see corresponding KRATOS_CATCH("") )

    // integration rule
    if(this->Has( INTEGRATION_ORDER ))
    {
        if(this->GetValue(INTEGRATION_ORDER) == 1)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 2)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 3)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 4)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 5)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "ConvectiveHeatFluxSubCondition Condition does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
    }
    else if(GetProperties().Has( INTEGRATION_ORDER ))
    {
        if(GetProperties()[INTEGRATION_ORDER] == 1)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 2)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 3)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 4)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 5)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "ConvectiveHeatFluxSubCondition Condition does not support for integration points", GetProperties()[INTEGRATION_ORDER])
    }
    else
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void RollerConstraint::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
void RollerConstraint::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
void RollerConstraint::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag)
{
    unsigned int number_of_nodes_disp = GetGeometry().size();
    unsigned int number_of_nodes_lambda = GetGeometry().size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = dim*number_of_nodes_disp + number_of_nodes_lambda;

    //resize as needed the LHS
    if ( CalculateStiffnessMatrixFlag ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) //calculation of the matrix is required
    {
        //resize the RHS=force vector if its size is not correct
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    const Matrix& Ncontainer_disp = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
    const Matrix& Ncontainer_lambda = Ncontainer_disp;

    double IntegrationWeight;

    Vector t1(3), t2(3), v3(3);
    array_1d<double, 3> u;
    double lambda;
    Vector Ndisp(number_of_nodes_disp);
    const Vector& Nlambda = Ndisp;

    unsigned int offset_row, offset_column;

    // loop over integration points
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();
        if (dim == 2)
            IntegrationWeight *= GetProperties()[THICKNESS];

        noalias( Ndisp ) = row( Ncontainer_disp, PointNumber );

        // compute the normal vector
        t1.clear(); // first tangential vector

        if (dim == 2)
        {
            for ( unsigned int n = 0; n < number_of_nodes_disp; n++ )
            {
                t1[0] += GetGeometry().GetPoint( n ).X0() * DN_DeContainer[PointNumber]( n, 0 );
                t1[1] += GetGeometry().GetPoint( n ).Y0() * DN_DeContainer[PointNumber]( n, 0 );
            }

            // calculating normal
            v3[0] = -t1[1];
            v3[1] = t1[0];
            v3[2] = 0.0;
        }
        else if (dim == 3)
        {
            t2.clear(); // second tangential vector

            for ( unsigned int n = 0; n < number_of_nodes_disp; n++ )
            {
                t1[0] += GetGeometry().GetPoint( n ).X0() * DN_DeContainer[PointNumber]( n, 0 );
                t1[1] += GetGeometry().GetPoint( n ).Y0() * DN_DeContainer[PointNumber]( n, 0 );
                t1[2] += GetGeometry().GetPoint( n ).Z0() * DN_DeContainer[PointNumber]( n, 0 );
                t2[0] += GetGeometry().GetPoint( n ).X0() * DN_DeContainer[PointNumber]( n, 1 );
                t2[1] += GetGeometry().GetPoint( n ).Y0() * DN_DeContainer[PointNumber]( n, 1 );
                t2[2] += GetGeometry().GetPoint( n ).Z0() * DN_DeContainer[PointNumber]( n, 1 );
            }

            // calculating normal
            v3[0] = t1[1] * t2[2] - t1[2] * t2[1];

            v3[1] = t1[2] * t2[0] - t1[0] * t2[2];

            v3[2] = t1[0] * t2[1] - t1[1] * t2[0];

            // double dA = sqrt( v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2] );
        }
        else
            KRATOS_ERROR << "RollerConstraint condition is not applicable for dim = " << dim;

        // interpolating the displacement
        u.clear();
        for (unsigned int i = 0; i < number_of_nodes_disp; ++i)
            noalias(u) += Ndisp[i] * GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);

        // interpolating the multiplier
        lambda = 0.0;
        for (unsigned int i = 0; i < number_of_nodes_lambda; ++i)
            lambda += Nlambda[i] * GetGeometry()[i].GetSolutionStepValue(LAMBDA);

        double aux = inner_prod(u, v3);

        // contribution to residual
        if (CalculateResidualVectorFlag)
        {
            for (unsigned int i = 0; i < number_of_nodes_lambda; ++i)
                rRightHandSideVector( i ) -= Nlambda[i] * aux * IntegrationWeight;

            offset_row = number_of_nodes_lambda;
            for (unsigned int i = 0; i < number_of_nodes_disp; ++i)
                for (unsigned int k = 0; k < dim; ++k)
                    rRightHandSideVector( i*dim + k + offset_row ) -= lambda * Ndisp[i] * v3[k] * IntegrationWeight;
        }

        // contribution to stiffness
        if (CalculateStiffnessMatrixFlag)
        {
            for (unsigned int i = 0; i < number_of_nodes_lambda; ++i)
            {
                // contribution to lambda-lambda block is zero

                // contribution to lambda-u block
                offset_column = number_of_nodes_lambda;
                for (unsigned int j = 0; j < number_of_nodes_disp; ++j)
                    for (unsigned int k = 0; k < dim; ++k)
                        rLeftHandSideMatrix( i, j*dim+k + offset_column ) +=
                            Nlambda[i] * Ndisp[j] * v3[k] * IntegrationWeight;
            }

            offset_row = number_of_nodes_lambda;
            for (unsigned int i = 0; i < number_of_nodes_disp; ++i)
            {
                for (unsigned int k = 0; k < dim; ++k)
                {
                    // contribution to u-lambda block
                    for (unsigned int j = 0; j < number_of_nodes_lambda; ++j)
                        rLeftHandSideMatrix( i*dim + k + offset_row, j ) +=
                            Nlambda[j] * Ndisp[i] * v3[k] * IntegrationWeight;

                    // contribution to u-u block is zero
                }
            }
        }
    } // end loop over integration points

    #ifdef ENABLE_BEZIER_GEOMETRY
    //clean the internal data of the geometry
    GetGeometry().Clean();
    #endif
}

//***********************************************************************************
//***********************************************************************************
void RollerConstraint::CalculateMassMatrix( MatrixType& rMassMatrix,
                              const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rMassMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void RollerConstraint::CalculateDampingMatrix( MatrixType& rDampingMatrix,
                              const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rDampingMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void RollerConstraint::EquationIdVector( EquationIdVectorType& rResult,
        const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY

    unsigned int number_of_nodes_disp = GetGeometry().size();
    unsigned int number_of_nodes_lambda = GetGeometry().size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes_lambda + dim*number_of_nodes_disp;

    if ( rResult.size() != MatSize )
        rResult.resize( MatSize, false );

    unsigned int cnt = 0;

    for ( unsigned int i = 0; i < number_of_nodes_lambda; i++ )
        rResult[cnt++] = GetGeometry()[i].GetDof( LAMBDA ).EquationId();

    for ( unsigned int i = 0; i < number_of_nodes_disp; i++ )
    {
        rResult[cnt++] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[cnt++] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        if (dim > 2)
            rResult[cnt++] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void RollerConstraint::GetDofList( DofsVectorType& ConditionalDofList,
                                        const ProcessInfo& CurrentProcessInfo) const
{
    ConditionalDofList.resize( 0 );

    unsigned int number_of_nodes_disp = GetGeometry().size();
    unsigned int number_of_nodes_lambda = GetGeometry().size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < number_of_nodes_lambda; i++ )
        ConditionalDofList.push_back( GetGeometry()[i].pGetDof( LAMBDA ) );

    for ( unsigned int i = 0; i < number_of_nodes_disp; i++ )
    {
        ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        if (dim > 2)
            ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

} // Namespace Kratos
