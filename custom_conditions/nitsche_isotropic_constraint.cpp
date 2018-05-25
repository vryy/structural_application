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
*   Date:                $Date: 9 Feb 2018 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/legacy_structural_app_vars.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "geometries/line_2d_3.h"
#include "custom_conditions/nitsche_isotropic_constraint.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
NitscheIsotropicConstraint::NitscheIsotropicConstraint( IndexType NewId,
        GeometryType::Pointer pGeometry) :
    Condition( NewId, pGeometry )
{
}

NitscheIsotropicConstraint::NitscheIsotropicConstraint( IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) :
    Condition( NewId, pGeometry, pProperties )
{
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************

Condition::Pointer NitscheIsotropicConstraint::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new NitscheIsotropicConstraint(NewId, GetGeometry().Create(ThisNodes),
                                pProperties));
}

/**
 * Destructor. Never to be called manually
 */
NitscheIsotropicConstraint::~NitscheIsotropicConstraint()
{
}


//************************************************************************************
//************************************************************************************

/**
 * Initialization of the element, called at the begin of each simulation.
 * Member variables and the Material law are initialized here
 */
void NitscheIsotropicConstraint::Initialize()
{
    KRATOS_TRY//EXCEPTION HANDLING (see corresponding KRATOS_CATCH("") )

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
            KRATOS_THROW_ERROR(std::logic_error, "NitscheIsotropicConstraint element does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
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
            KRATOS_THROW_ERROR(std::logic_error, "NitscheIsotropicConstraint element does not support for integration points", GetProperties()[INTEGRATION_ORDER])
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
void NitscheIsotropicConstraint::CalculateRightHandSide( VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
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
void NitscheIsotropicConstraint::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
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
void NitscheIsotropicConstraint::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag)
{
    Element::Pointer pElement = this->GetValue(ASSOCIATED_ELEMENT);
//KRATOS_WATCH(pElement->Id())
    GeometryType& rGeometry = pElement->GetGeometry();

    unsigned int number_of_nodes = rGeometry.size();
    unsigned int dim = rGeometry.WorkingSpaceDimension();
    unsigned int strain_size = dim * ( dim + 1 ) / 2;

    //this is the size of the elements stiffness matrix/force vector
    unsigned int mat_size = rGeometry.size() * dim;

    //Initialize local variables
    Matrix B( strain_size, mat_size );
    Matrix TanC( strain_size, strain_size );
    Vector StrainVector( strain_size );
    Vector StressVector( strain_size );
    Matrix DN_De( number_of_nodes, dim );
    Matrix DN_DX( number_of_nodes, dim );
    Vector Ncontainer( number_of_nodes );
    Matrix Nmatrix( dim, mat_size );
    Matrix Noperator( dim, strain_size );
    Vector CurrentDisp( number_of_nodes * dim );
    Matrix J0( dim, dim );
    Matrix InvJ0( dim, dim );
    double DetJ0;
    Vector t1(dim), t2(dim), n(dim);
    Matrix Aux1, Aux2, Aux3;
    Vector aux;
    

    //resize the LHS=StiffnessMatrix if its size is not correct
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        //resize the RHS=force vector if its size is not correct
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    rGeometry.Initialize(mThisIntegrationMethod);
    #endif

    //initializing the Jacobian in the reference configuration
    Matrix DeltaPosition(rGeometry.size(), 3);
    for ( unsigned int node = 0; node < rGeometry.size(); ++node )
        noalias( row( DeltaPosition, node ) ) = rGeometry[node].Coordinates() - rGeometry[node].GetInitialPosition();

    //Current displacements
    for ( unsigned int node = 0; node < rGeometry.size(); ++node )
    {
        CurrentDisp( node*dim ) = rGeometry[node].GetSolutionStepValue( DISPLACEMENT_X );
        CurrentDisp( node*dim + 1 ) = rGeometry[node].GetSolutionStepValue( DISPLACEMENT_Y );
        if (dim == 3)
            CurrentDisp( node*dim + 2 ) = rGeometry[node].GetSolutionStepValue( DISPLACEMENT_Z );
    }

    //reading integration points and local gradients on the surface
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DN_De_surf = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    // material properties
    const Matrix& H = this->GetValue(CONSTRAINT_MATRIX);
    const Vector& g = this->GetValue(CONSTRAINT_VECTOR);
    const double& E = GetProperties()[YOUNG_MODULUS];
    const double& Nu = GetProperties()[POISSON_RATIO];
    const double& Thickness = GetProperties()[THICKNESS];
    const double& tau = GetProperties()[STABILISATION_FACTOR];

    if (dim == 2)
        this->CalculateElasticMatrix2DPlaneStrain( TanC, E, Nu );
    else if (dim == 3)
        this->CalculateElasticMatrix3D( TanC, E, Nu );

    CoordinatesArrayType CurrentGlobalCoords( ZeroVector( 3 ) );
    CoordinatesArrayType CurrentLocalCoords( ZeroVector( 3 ) );

    const unsigned int number_of_constraints = H.size1();

    Aux1.resize( mat_size, number_of_constraints, false );
    Aux2.resize( mat_size, mat_size, false );
    Aux3.resize( number_of_constraints, mat_size, false );
    aux.resize( number_of_constraints, false );
//KRATOS_WATCH(CurrentDisp)
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        CurrentGlobalCoords = GetGeometry().GlobalCoordinates( CurrentGlobalCoords, integration_points[PointNumber], DeltaPosition );
//KRATOS_WATCH(CurrentGlobalCoords)
        CurrentLocalCoords = rGeometry.PointLocalCoordinates( CurrentLocalCoords, CurrentGlobalCoords, DeltaPosition ); 
//KRATOS_WATCH(CurrentLocalCoords)

        Ncontainer = rGeometry.ShapeFunctionsValues( Ncontainer, CurrentLocalCoords );

        J0 = rGeometry.Jacobian( J0, CurrentLocalCoords, DeltaPosition );

        rGeometry.ShapeFunctionsLocalGradients( DN_De, CurrentLocalCoords );

        MathUtils<double>::InvertMatrix( J0, InvJ0, DetJ0 );
        noalias( DN_DX ) = prod( DN_De, InvJ0 );

        CalculateBoperator( B, DN_DX );

        if (dim == 2)
        {
            noalias( t1 ) = ZeroVector( 2 ); //tangential vector
            for ( unsigned int n = 0; n < GetGeometry().size(); ++n )
            {
                t1[0] += GetGeometry().GetPoint( n ).X0() * DN_De_surf[PointNumber]( n, 0 );
                t1[1] += GetGeometry().GetPoint( n ).Y0() * DN_De_surf[PointNumber]( n, 0 );
            }
            n( 0 ) = -t1( 1 );
            n( 1 ) = t1( 0 );
            n *= 1.0/norm_2(n);
        }
        else if (dim == 3)
        {
            noalias( t1 ) = ZeroVector( 3 ); //first tangential vector
            noalias( t2 ) = ZeroVector( 3 ); //second tangential vector
            for ( unsigned int n = 0; n < GetGeometry().size(); ++n )
            {
                t1[0] += GetGeometry().GetPoint( n ).X0() * DN_De_surf[PointNumber]( n, 0 );
                t1[1] += GetGeometry().GetPoint( n ).Y0() * DN_De_surf[PointNumber]( n, 0 );
                t1[2] += GetGeometry().GetPoint( n ).Z0() * DN_De_surf[PointNumber]( n, 0 );
                t2[0] += GetGeometry().GetPoint( n ).X0() * DN_De_surf[PointNumber]( n, 1 );
                t2[1] += GetGeometry().GetPoint( n ).Y0() * DN_De_surf[PointNumber]( n, 1 );
                t2[2] += GetGeometry().GetPoint( n ).Z0() * DN_De_surf[PointNumber]( n, 1 );
            }
            noalias( n ) = MathUtils<double>::CrossProduct( t1, t2 );
            n *= 1.0/norm_2(n);
        }

        CalculateNoperator( dim, Noperator, n );

        CalculateInterpolationOperator( dim, Nmatrix, Ncontainer );

        if ( CalculateResidualVectorFlag == true || CalculateStiffnessMatrixFlag == true )
        {
            double IntToReferenceWeight = integration_points[PointNumber].Weight();

            if ( dim == 2 ) IntToReferenceWeight *= Thickness;

            noalias( Aux1 ) = prod( trans(B), Matrix( prod( TanC, Matrix ( prod( trans(Noperator), trans(H) ) ) ) ) );
            noalias( Aux3 ) = prod( H, Nmatrix );
            noalias( Aux2 ) = prod( trans(Aux3), trans(Aux1) );
//KRATOS_WATCH(Aux1)
//KRATOS_WATCH(Aux2)
//KRATOS_WATCH(H)
//KRATOS_WATCH(Nmatrix)
//KRATOS_WATCH(Aux3)

            if ( CalculateStiffnessMatrixFlag )
            {
                noalias( rLeftHandSideMatrix ) += ( trans(Aux2) + Aux2 + tau*prod(trans(Aux3), Aux3) ) * DetJ0 * IntToReferenceWeight;
            }

            if ( CalculateResidualVectorFlag )
            {
                noalias( aux ) = prod(Aux3, CurrentDisp) - g;
                noalias( rRightHandSideVector ) -= (
                    prod( Aux1, aux )
                    + prod( Aux2, CurrentDisp )
                    + tau * prod( trans(Aux3), aux )
                    ) * DetJ0 * IntToReferenceWeight;

//                noalias( rRightHandSideVector ) -= prod( rLeftHandSideMatrix, CurrentDisp );
            }

            
        }
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    // clean the geometry
    rGeometry.Clean();
    #endif
//    KRATOS_WATCH(rRightHandSideVector)
//    KRATOS_WATCH(rLeftHandSideMatrix)
//    KRATOS_WATCH("----------------------------------")
}

//************************************************************************************
//************************************************************************************

void NitscheIsotropicConstraint::EquationIdVector( EquationIdVectorType& rResult,
        ProcessInfo& CurrentProcessInfo
                                            )
{
    //determining size of DOF list
    //dimension of space
    Element::Pointer pElement = this->GetValue(ASSOCIATED_ELEMENT);
    GeometryType& rGeometry = pElement->GetGeometry();
    unsigned int dim = rGeometry.WorkingSpaceDimension();
    rResult.resize(rGeometry.size()*dim, false);
    for( unsigned int i = 0; i < rGeometry.size(); ++i )
    {
        rResult[dim*i] = rGeometry[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[dim*i+1] = rGeometry[i].GetDof(DISPLACEMENT_Y).EquationId();
        if (dim == 3)
            rResult[dim*i+2] = rGeometry[i].GetDof(DISPLACEMENT_Z).EquationId();
    }
}

//************************************************************************************
//************************************************************************************

void NitscheIsotropicConstraint::GetDofList( DofsVectorType& ConditionalDofList,
                                        ProcessInfo& CurrentProcessInfo)
{
    //determining size of DOF list
    //dimension of space
    Element::Pointer pElement = this->GetValue(ASSOCIATED_ELEMENT);
    GeometryType& rGeometry = pElement->GetGeometry();
    unsigned int dim = rGeometry.WorkingSpaceDimension();
    ConditionalDofList.resize(rGeometry.size()*dim);
    for( unsigned int i = 0; i < rGeometry.size(); ++i )
    {
        ConditionalDofList[dim*i] = rGeometry[i].pGetDof(DISPLACEMENT_X);
        ConditionalDofList[dim*i+1] = rGeometry[i].pGetDof(DISPLACEMENT_Y);
        if (dim == 3)
            ConditionalDofList[dim*i+2] = rGeometry[i].pGetDof(DISPLACEMENT_Z);
    }
}

//************************************************************************************
//************************************************************************************

void NitscheIsotropicConstraint::CalculateInterpolationOperator( const unsigned int& Dim, Matrix& N_Operator, Vector& Ncontainer )
{
    const unsigned int number_of_nodes = Ncontainer.size();

    noalias( N_Operator ) = ZeroMatrix( Dim, Dim*number_of_nodes );

    for ( unsigned int node = 0; node < number_of_nodes; ++node )
    {
        for ( unsigned int dim = 0; dim < Dim; ++dim )
        {
            N_Operator( dim, Dim*node + dim ) = Ncontainer( node );
        }
    }
}

void NitscheIsotropicConstraint::CalculateNoperator( const unsigned int& Dim, Matrix& N_Operator, const Vector& n )
{
    if ( Dim == 2 )
    {
        noalias( N_Operator ) = ZeroMatrix( 2, 3 );

        N_Operator( 0, 0 ) = n( 0 );
        N_Operator( 1, 1 ) = n( 1 );
        N_Operator( 0, 2 ) = n( 1 );

        N_Operator( 1, 2 ) = n( 0 );
    }
    else if ( Dim == 3 )
    {
        noalias( N_Operator ) = ZeroMatrix( 3, 6 );

        N_Operator( 0, 0 ) = n( 0 );
        N_Operator( 1, 1 ) = n( 1 );
        N_Operator( 2, 2 ) = n( 2 );

        N_Operator( 0, 3 ) = n( 1 );
        N_Operator( 0, 5 ) = n( 2 );

        N_Operator( 1, 3 ) = n( 0 );
        N_Operator( 1, 4 ) = n( 2 );

        N_Operator( 2, 4 ) = n( 1 );
        N_Operator( 2, 5 ) = n( 0 );
    }
}

void NitscheIsotropicConstraint::CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX )
{
    KRATOS_TRY

    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int strain_size = dim * (dim + 1) / 2;
    const unsigned int number_of_nodes = DN_DX.size1();

    noalias( B_Operator ) = ZeroMatrix( strain_size, number_of_nodes * dim );

    if(dim == 2)
    {
        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            B_Operator( 0, i*2 ) = DN_DX( i, 0 );
            B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 );
            B_Operator( 2, i*2 ) = DN_DX( i, 1 );
            B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
        }
    }
    else if(dim == 3)
    {
        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            B_Operator( 0, i*3 ) = DN_DX( i, 0 );
            B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
            B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 );
            B_Operator( 3, i*3 ) = DN_DX( i, 1 );
            B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
            B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
            B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
            B_Operator( 5, i*3 ) = DN_DX( i, 2 );
            B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
        }
    }

    KRATOS_CATCH( "" )
}

void NitscheIsotropicConstraint::CalculateElasticMatrix2DPlaneStrain( Matrix& C, const double& E, const double& NU )
{
    double c1 = E * ( 1.00 - NU ) / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
    double c2 = E * NU / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
    double c3 = 0.5 * E / ( 1 + NU );

    C( 0, 0 ) = c1;
    C( 0, 1 ) = c2;
    C( 0, 2 ) = 0.0;
    C( 1, 0 ) = c2;
    C( 1, 1 ) = c1;
    C( 1, 2 ) = 0.0;
    C( 2, 0 ) = 0.0;
    C( 2, 1 ) = 0.0;
    C( 2, 2 ) = c3;
}

void NitscheIsotropicConstraint::CalculateElasticMatrix3D( Matrix& C, const double& E, const double& NU )
{
    //setting up material matrix
    double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
    double c2 = c1 * ( 1 - NU );
    double c3 = c1 * NU;
    double c4 = c1 * 0.5 * ( 1 - 2 * NU );
    //filling material matrix
    C( 0, 0 ) = c2;
    C( 0, 1 ) = c3;
    C( 0, 2 ) = c3;
    C( 0, 3 ) = 0.0;
    C( 0, 4 ) = 0.0;
    C( 0, 5 ) = 0.0;
    C( 1, 0 ) = c3;
    C( 1, 1 ) = c2;
    C( 1, 2 ) = c3;
    C( 1, 3 ) = 0.0;
    C( 1, 4 ) = 0.0;
    C( 1, 5 ) = 0.0;
    C( 2, 0 ) = c3;
    C( 2, 1 ) = c3;
    C( 2, 2 ) = c2;
    C( 2, 3 ) = 0.0;
    C( 2, 4 ) = 0.0;
    C( 2, 5 ) = 0.0;
    C( 3, 0 ) = 0.0;
    C( 3, 1 ) = 0.0;
    C( 3, 2 ) = 0.0;
    C( 3, 3 ) = c4;
    C( 3, 4 ) = 0.0;
    C( 3, 5 ) = 0.0;
    C( 4, 0 ) = 0.0;
    C( 4, 1 ) = 0.0;
    C( 4, 2 ) = 0.0;
    C( 4, 3 ) = 0.0;
    C( 4, 4 ) = c4;
    C( 4, 5 ) = 0.0;
    C( 5, 0 ) = 0.0;
    C( 5, 1 ) = 0.0;
    C( 5, 2 ) = 0.0;
    C( 5, 3 ) = 0.0;
    C( 5, 4 ) = 0.0;
    C( 5, 5 ) = c4;
}

int NitscheIsotropicConstraint::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int dim = this->GetGeometry().WorkingSpaceDimension();

    if ( this->Id() < 1 )
    {
        KRATOS_THROW_ERROR( std::logic_error, "Element found with Id 0 or negative, Id() =", Id() );
    }

    if ( this->GetProperties().Has( YOUNG_MODULUS ) == false )
    {
        KRATOS_THROW_ERROR( std::logic_error, "YOUNG_MODULUS not provided for property", this->GetProperties().Id() )
    }

    if ( this->GetProperties().Has( POISSON_RATIO ) == false )
    {
        KRATOS_THROW_ERROR( std::logic_error, "POISSON_RATIO not provided for property", this->GetProperties().Id() )
    }

    if ( dim == 2 )
        if ( this->GetProperties().Has( THICKNESS ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "POISSON_RATIO not provided for property", this->GetProperties().Id() )
        }

//    if ( this->Has( CONSTRAINT_MATRIX ) == false )
//    {
//        KRATOS_THROW_ERROR( std::logic_error, "CONSTRAINT_MATRIX not provided for condition", this->Id() )
//    }

//    if ( this->Has( CONSTRAINT_VECTOR ) == false )
//    {
//        KRATOS_THROW_ERROR( std::logic_error, "CONSTRAINT_VECTOR not provided for condition", this->Id() )
//    }

//    const Matrix& H = this->GetValue(CONSTRAINT_MATRIX);
//    const Vector& g = this->GetValue(CONSTRAINT_VECTOR);
//    const unsigned int number_of_constraints = H.size1();

//    if ( H.size2() < dim )
//        KRATOS_THROW_ERROR(std::logic_error, "The constraint matrix does not have correct dimension", "")

//    if ( g.size() < number_of_constraints )
//        KRATOS_THROW_ERROR(std::logic_error, "The constraint vector does not have correct size", "")

    return 0;

    KRATOS_CATCH( "" );
}

} // Namespace Kratos
