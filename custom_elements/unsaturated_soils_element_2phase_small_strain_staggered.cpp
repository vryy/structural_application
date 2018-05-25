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
*   Date:                $Date: 7 Jan 2016 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes
#include <typeinfo>


// External includes


// Project includes
#include "custom_elements/unsaturated_soils_element_2phase_small_strain_staggered.h"
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/prism_3d_6.h"
#include "structural_application.h"

#define ENABLE_DEBUG_CONSTITUTIVE_LAW

namespace Kratos
{

UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::UnsaturatedSoilsElement_2phase_SmallStrain_Staggered( IndexType NewId,
        GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
    mIsInitialized = false;
}

//************************************************************************************
//************************************************************************************
UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::UnsaturatedSoilsElement_2phase_SmallStrain_Staggered( IndexType NewId,
        GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
    mIsInitialized = false;
}

Element::Pointer UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UnsaturatedSoilsElement_2phase_SmallStrain_Staggered( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UnsaturatedSoilsElement_2phase_SmallStrain_Staggered( NewId, pGeom, pProperties ) );
}

UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::~UnsaturatedSoilsElement_2phase_SmallStrain_Staggered()
{
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::Initialize()
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if ( mIsInitialized )
    {
        //Set Up Initial displacement for StressFreeActivation of Elements
        mInitialDisp.resize( GetGeometry().size(), dim, false );

        for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
            for ( unsigned int i = 0; i < 3; i++ )
                mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

        return;
    }

    if(GetProperties().Has( INTEGRATION_ORDER ) == true)
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
            KRATOS_THROW_ERROR(std::logic_error, "Does not support for more integration points", *this)
    }

    // number of integration points used, mThisIntegrationMethod refers to the
    // integration method defined in the constructor
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    // initializing the Jacobian, the inverse Jacobian and Jacobians determinant in the reference configuration
    mInvJ0.resize( integration_points.size() );

    for ( unsigned int i = 0; i < integration_points.size(); i++ )
    {
        mInvJ0[i].resize( dim, dim, false );
        noalias( mInvJ0[i] ) = ZeroMatrix( dim, dim );
    }

    mDetJ0.resize( integration_points.size() );

    noalias( mDetJ0 ) = ZeroVector( integration_points.size() );

    mTotalDomainInitialSize = 0.0;
    
    GeometryType::JacobiansType J0( integration_points.size() );

    // calculating the Jacobian
    J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod );

    // calculating the inverse J0
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();
        //calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );
        //calculating the total area/volume
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }

    //Set Up Initial displacement for StressFreeActivation of Elements
    mInitialDisp.resize( GetGeometry().size(), dim, false );

    for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
        for ( unsigned int i = 0; i < 3; i++ )
            mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

    //Constitutive Law initialisation
    if ( mConstitutiveLawVector.size() != integration_points.size() )
    {
        mConstitutiveLawVector.resize( integration_points.size() );
        mReferencePressures.resize( integration_points.size() );
        InitializeMaterial();
    }

    mIsInitialized = true;

    KRATOS_CATCH( "" )
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::ResetConstitutiveLaw()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            Vector dummy;
//            dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i );
            mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(),  dummy);
        }
    }

    KRATOS_CATCH( "" )
}


/**
* THIS method is called from the scheme at the start of each solution step
* @param rCurrentProcessInfo
*/
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
    if(mConstitutiveLawVector[0]->Has(CURRENT_STRAIN_VECTOR))
    {
        std::vector<Vector> Values;
        this->GetValueOnIntegrationPoints(CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo);
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
        }
    }

    for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
    {
        Vector dummy;
//        dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point );
        mConstitutiveLawVector[Point]->InitializeSolutionStep( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
{
    if(mConstitutiveLawVector[0]->Has(CURRENT_STRAIN_VECTOR))
    {
        std::vector<Vector> Values;
        this->GetValueOnIntegrationPoints(CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo);
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
        }
    }

    for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
    {
        Vector dummy;
//        dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point );
        mConstitutiveLawVector[Point]->InitializeNonLinearIteration( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
    }
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double Weight;

    double DetJ;

    Vector N(number_of_nodes);

    #ifdef ENABLE_BEZIER_GEOMETRY
    // initialize the geometry
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    // extract the integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    // extract the shape function values
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

    // extract the shape function local gradients
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    // compute Jacobian
    GeometryType::JacobiansType J(integration_points.size());
    J = GetGeometry().Jacobian(J, mThisIntegrationMethod);

    if(rCurrentProcessInfo[FRACTIONAL_STEP] == 0) // solve for displacement
    {
        // variable declaration
        unsigned int mat_size = dim * number_of_nodes;

        unsigned int strain_size = dim * (dim + 1) / 2;

        Matrix B( strain_size, mat_size );

        Matrix TanC( strain_size, strain_size );

        Vector StrainVector( strain_size );

        Vector StressVector( strain_size );

        Matrix DN_DX( number_of_nodes, dim );

        Matrix CurrentDisp( number_of_nodes, dim );

        double capillaryPressure;

        double waterPressure;

        double porosity;

        double density;

        // Current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            CurrentDisp( node, 0 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_X );
            CurrentDisp( node, 1 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_Y );
            CurrentDisp( node, 2 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_Z );
        }

        // resize as needed the LHS
        if ( CalculateStiffnessMatrixFlag == true )
        {
            if ( rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size )
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
        }

        // resize as needed the RHS
        if ( CalculateResidualVectorFlag == true )
        {
            if ( rRightHandSideVector.size() != mat_size )
                rRightHandSideVector.resize( mat_size, false );

            noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
        }

        /////////////////////////////////////////////////////////////////////////
        //// Compute the B-dilatational operator
        //// Reference: Thomas Hughes, The Finite Element Method
        /////////////////////////////////////////////////////////////////////////
        Matrix Bdil_bar;
        if(GetProperties().Has(IS_BBAR))
        {
            if(GetProperties()[IS_BBAR] == true)
            {
                Bdil_bar.resize(number_of_nodes, dim);
                noalias(Bdil_bar) = ZeroMatrix(number_of_nodes, dim);
                for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
                {
                    Weight = integration_points[PointNumber].Weight();
                    DetJ   = mDetJ0[PointNumber];
                    noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );
    //                KRATOS_WATCH(DN_DX_DISP)
    //                KRATOS_WATCH(Weight)
    //                KRATOS_WATCH(DetJ)
                    noalias(Bdil_bar) += DN_DX * Weight * DetJ;
                }
                Bdil_bar /= mTotalDomainInitialSize;
            }
        }

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space sum_(beta=0)^(number of quadrature points) ////
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

            Weight = integration_points[PointNumber].Weight();

            DetJ = mDetJ0[PointNumber];

            // Shape Functions on current spatial quadrature point
            noalias( N ) = row( Ncontainer, PointNumber );

            // Initializing B_Operator at the current integration point
            if(GetProperties().Has(IS_BBAR))
            {
                if(GetProperties()[IS_BBAR] == true)
                    CalculateBBaroperator( B, DN_DX, Bdil_bar );
                else
                    CalculateBoperator( B, DN_DX );
            }
            else
                CalculateBoperator( B, DN_DX );

            // Calculate the current strain vector using the B-Operator
            CalculateStrain( B, CurrentDisp, StrainVector );

            noalias(StressVector) = ZeroVector(6);

            GetPressures( N, capillaryPressure, waterPressure );
            // REMARKS: comment this to maintain consistent linearisation
//            if( waterPressure < 0.0 ) // avoid going into saturated state
//            {
//                waterPressure = 0.0;
//            }

            porosity = GetPorosity( DN_DX );
            if ( porosity > 1.0 || porosity < 0.0 )
            {
                KRATOS_THROW_ERROR(std::logic_error, "porosity is ", porosity)
            }

            density = GetAveragedDensity( capillaryPressure, porosity );

            // Set suction in const. law
            mConstitutiveLawVector[PointNumber]->SetValue( SUCTION, -waterPressure,  rCurrentProcessInfo );

            // retrieve the material response
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(StrainVector,
                    ZeroMatrix( 1 ), StressVector, TanC, rCurrentProcessInfo,
                    GetProperties(), GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod), PointNumber ),
                    true, 1, true );

            // calculate effective stress
            double capillaryPressure = -waterPressure;

            double saturation = GetSaturation( capillaryPressure );

            for ( unsigned int i = 0; i < 3; ++i )
                StressVector( i ) -= ( saturation * waterPressure );

            if ( CalculateStiffnessMatrixFlag == true )
            {
                CalculateStiffnesMatrixUU( rLeftHandSideMatrix, TanC, B, DN_DX, N, density, capillaryPressure, Weight, DetJ );
            }

            if ( CalculateResidualVectorFlag == true )
            {
                //Calculation of spatial Loadvector
                AddBodyForcesToRHSVectorU( rRightHandSideVector, N, density, Weight, DetJ );
                AddInternalForcesToRHSU( rRightHandSideVector, B, StressVector, Weight, DetJ );
            }
        }
        ///////////////////////////////////////////////////////////////////////
        // END Integration in space sum_(beta=0)^(number of quadrature points)
        ///////////////////////////////////////////////////////////////////////
    }
    else if(rCurrentProcessInfo[FRACTIONAL_STEP] == 1) // solve for water pressure
    {
        unsigned int mat_size = number_of_nodes;

        Matrix TanW( dim, dim );

        Matrix DN_DX( number_of_nodes, dim );

        double capillaryPressure;

        double waterPressure;

        // resize as needed the LHS
        if ( CalculateStiffnessMatrixFlag == true )
        {
            if ( rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size )
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
        }

        // resize as needed the RHS
        if ( CalculateResidualVectorFlag == true )
        {
            if ( rRightHandSideVector.size() != mat_size )
                rRightHandSideVector.resize( mat_size, false );

            noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
        }

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space sum_(beta=0)^(number of quadrature points) ////
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

            Weight = integration_points[PointNumber].Weight();

            DetJ = mDetJ0[PointNumber];

            // Shape Functions on current spatial quadrature point
            noalias( N ) = row( Ncontainer, PointNumber );

            GetPressures( N, capillaryPressure, waterPressure );
            // REMARKS: comment this to maintain consistent linearisation
//            if( waterPressure < 0.0 ) // avoid going into saturated state
//            {
//                waterPressure = 0.0;
//            }

            if ( CalculateStiffnessMatrixFlag == true )
            {
                CalculateStiffnesMatrixWW( rLeftHandSideMatrix, DN_DX, DN_DX, N, capillaryPressure, Weight, DetJ );
            }

            if ( CalculateResidualVectorFlag == true )
            {
                AddInternalForcesToRHSW( rRightHandSideVector, DN_DX, DN_DX, N, capillaryPressure, Weight, DetJ );
            }
        }
        ///////////////////////////////////////////////////////////////////////
        // END Integration in space sum_(beta=0)^(number of quadrature points)
        ///////////////////////////////////////////////////////////////////////
    }
    else // solve for nothing (should not happen)
    {
        noalias(rLeftHandSideMatrix) = ZeroMatrix(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
        noalias(rRightHandSideVector) = ZeroVector(rRightHandSideVector.size());
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    // clean the geometry internal data
    GetGeometry().Clean();
    #endif

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateRightHandSide( VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

////************************************************************************************
////************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    DampMatrix(rMassMatrix, rCurrentProcessInfo);
}

////************************************************************************************
////************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateDampingMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if(rCurrentProcessInfo[FRACTIONAL_STEP] == 0) // solve for displacement
    {
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        //resizing as needed the LHS
        if(rDampMatrix.size1() != mat_size || rDampMatrix.size2() != mat_size)
            rDampMatrix.resize(mat_size, mat_size);
        rDampMatrix = ZeroMatrix(mat_size, mat_size);
    }
    if(rCurrentProcessInfo[FRACTIONAL_STEP] == 1) // solve for pressure
    {
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int mat_size = number_of_nodes;
        double Weight;
        double DetJ;

        #ifdef ENABLE_BEZIER_GEOMETRY
        // initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De =
            GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        double capillaryPressure;

        double waterPressure;

        Matrix DN_DX( number_of_nodes, 3 );

        Vector N( number_of_nodes );

        //resizing as needed the LHS
        if(rDampMatrix.size1() != mat_size || rDampMatrix.size2() != mat_size)
            rDampMatrix.resize(mat_size, mat_size);
        rDampMatrix = ZeroMatrix(mat_size, mat_size);

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space sum_(beta=0)^(number of quadrature points)
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

            Weight = integration_points[PointNumber].Weight();

            DetJ = mDetJ0[PointNumber];

            noalias( N ) = row( Ncontainer, PointNumber );

            GetPressures( N, capillaryPressure, waterPressure );

            CalculateDampingMatrixWW( rDampMatrix, DN_DX, N, capillaryPressure, Weight, DetJ );
        }
        ///////////////////////////////////////////////////////////////////////
        // END Integration in space sum_(beta=0)^(number of quadrature points)
        ///////////////////////////////////////////////////////////////////////

        #ifdef ENABLE_BEZIER_GEOMETRY
        // finalize the geometry
        GetGeometry().Clean();
        #endif
    }

    KRATOS_CATCH( "" )
}

////************************************************************************************
////************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::FinalizeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
{
    if(mConstitutiveLawVector[0]->Has(CURRENT_STRAIN_VECTOR))
    {
        std::vector<Vector> Values;
        this->GetValueOnIntegrationPoints(CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo);
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
        }
    }

    for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
    {
        Vector dummy;
//        dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point );
        mConstitutiveLawVector[Point]->FinalizeNonLinearIteration( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
    }
}

////************************************************************************************
////************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
    if(mConstitutiveLawVector[0]->Has(CURRENT_STRAIN_VECTOR))
    {
        std::vector<Vector> Values;
        this->GetValueOnIntegrationPoints(CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo);
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
        }
    }

    for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
    {
        Vector dummy;
//        dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point );
        mConstitutiveLawVector[Point]->FinalizeSolutionStep( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
    }

//    Vector Dummy_Vector( 9 );
// 
//    noalias( Dummy_Vector ) = mConstitutiveLawVector[0]->GetValue( INTERNAL_VARIABLES, Dummy_Vector );
// 
//    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
//    {
//        GetGeometry()[i].GetSolutionStepValue( MOMENTUM_X ) = Dummy_Vector( 0 );
//        GetGeometry()[i].GetSolutionStepValue( MOMENTUM_Y ) = Dummy_Vector( 1 );
//        GetGeometry()[i].GetSolutionStepValue( MOMENTUM_Z ) = Dummy_Vector( 2 );
//        GetGeometry()[i].GetSolutionStepValue( PRESSURE ) = Dummy_Vector( 3 );
//        GetGeometry()[i].GetSolutionStepValue( ERROR_RATIO ) = Dummy_Vector( 4 );
//        GetGeometry()[i].GetSolutionStepValue( EXCESS_PORE_WATER_PRESSURE ) = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) - 9.81 * 1000.0 * ( 20.0 - GetGeometry()[i].Z() );
//    }
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateOnIntegrationPoints( const Variable<double >& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();

    #ifdef ENABLE_BEZIER_GEOMETRY
    // initialize the geometry
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if ( Output.size() != integration_points.size() )
        Output.resize( integration_points.size() );

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

    Vector N( number_of_nodes );

    double capillaryPressure;

    double waterPressure;

    double saturation;

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        // Shape Functions on current spatial quadrature point
        if ( N.size() != number_of_nodes )
            N.resize( number_of_nodes );

        noalias( N ) = row( Ncontainer, PointNumber );
        GeometryType::CoordinatesArrayType gp_position;
        gp_position = GetGeometry().GlobalCoordinates( gp_position, GetGeometry().IntegrationPoints()[PointNumber] );

        GetPressures( N, capillaryPressure, waterPressure );

        saturation = GetSaturation( capillaryPressure );

        if ( rVariable == SATURATION )
        {
            Output[PointNumber] = saturation;
        }

        if ( rVariable == WATER_PRESSURE )
        {
            Output[PointNumber] = waterPressure;
        }

        if ( rVariable == EXCESS_PORE_WATER_PRESSURE )
        {
            Output[PointNumber] = waterPressure - mReferencePressures[PointNumber];
        }

        if ( rVariable == AIR_PRESSURE )
        {
            Output[PointNumber] = 0.0;
        }
    /////////////////////////////////////////////////////////////////////////
    //// End Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    // finalize the geometry
    GetGeometry().Clean();
    #endif

    KRATOS_CATCH( "" )
}

/**
* Calculate Vector Variables at each integration point, used for postprocessing etc.
* @param rVariable Global name of the variable to be calculated
* @param output Vector to store the values on the qudrature points, output of the method
* @param rCurrentProcessInfo
*/
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo )
{
    GetValueOnIntegrationPoints( rVariable, Output, rCurrentProcessInfo );
}


//************************************************************************************
//************************************************************************************

inline void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateAndAddExtForceContribution( const Vector& N, const ProcessInfo& CurrentProcessInfo, Vector& BodyForce, VectorType& rRightHandSideVector, double weight )
{
    KRATOS_TRY
    // TODO
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::EquationIdVector( EquationIdVectorType& rResult,
        ProcessInfo& CurrentProcessInfo )
{
    DofsVectorType ElementalDofList;
    this->GetDofList(ElementalDofList, CurrentProcessInfo);
    rResult.resize(ElementalDofList.size());
    for(unsigned int i = 0; i < ElementalDofList.size(); ++i)
        rResult[i] = ElementalDofList[i]->EquationId();
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetDofList( DofsVectorType& ElementalDofList,
        ProcessInfo& CurrentProcessInfo )
{
    ElementalDofList.resize( 0 );

    if(CurrentProcessInfo[FRACTIONAL_STEP] == 0) // solve for displacement
    {
        for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }
    else if(CurrentProcessInfo[FRACTIONAL_STEP] == 1) // solve for pressure
    {
        for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ) );
        }
    }
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetValuesVector( Vector& values, int Step )
{
    // TODO
}

//************************************************************************************
//CALCULATE EXTERNAL FORCEVECTORS DISPLACEMENT****************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::AddBodyForcesToRHSVectorU( Vector& R, Vector& N_DISP, double density, double Weight, double detJ )
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    Vector gravity( dim );

    double density = 0.0;

    if ( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
    {
        noalias( gravity ) = GetValue( GRAVITY );
        density = GetValue( DENSITY );
    }
    else
    {
        noalias( gravity ) = GetProperties()[GRAVITY];
        density = GetProperties()[DENSITY];
    }

    for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
    {
        for ( unsigned int i = 0; i < dim; ++i )
        {
            R( prim*dim + i ) +=
                N_DISP( prim ) * density * gravity( i ) * Weight * detJ;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//CALCULATE INTERNAL FORCEVECTORS DISPLACEMENT****************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::AddInternalForcesToRHSU( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
{
    KRATOS_TRY

    noalias( R ) -= detJ * Weight * prod( trans( B_Operator ), StressVector );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//CALCULATE STIFFNESS MATRICES DISPLACEMENT*******************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateStiffnesMatrixUU( Matrix& K, const
        Matrix& C, const Matrix& B_Operator, const Matrix& DN_DX_DISP, Vector& N_DISP,
        double density, double capillaryPressure, double Weight, double detJ )
{
    KRATOS_TRY
    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();

    Vector gravity( dim );

    double density_soil = 0.0;

    if ( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
    {
        noalias( gravity ) = GetValue( GRAVITY );
        density_soil = GetValue( DENSITY );
    }
    else
    {
        noalias( gravity ) = GetProperties()[GRAVITY];
        density_soil = GetProperties()[DENSITY];
    }
    double density_water = GetValue( DENSITY_WATER );
    double density_air = GetValue( DENSITY_AIR );

    double saturation = GetSaturation( capillaryPressure );

    double porosity_divu = 0.0;
//    double porosity_divu= GetDerivativeDPorosityDDivU(DN_DX_DISP);

    double DrhoDdivU =  porosity_divu * ( -density_soil + ( 1 - saturation ) * density_air + saturation * density_water );

    for ( unsigned int prim = 0; prim < number_of_nodes; ++prim )
    {
        for ( unsigned int i = 0; i < dim; ++i )
        {
            for ( unsigned int sec = 0; sec < number_of_nodes; ++sec )
            {
                for ( unsigned int j = 0; j < dim; ++j )
                {
                    K( prim*dim + i, sec*dim + j ) += N_DISP( prim ) * DrhoDdivU * gravity( i ) * DN_DX_DISP( sec, j ) * Weight * detJ;
                }
            }
        }
    }

    noalias( K ) += prod( trans( B_Operator ), ( Weight * detJ ) * Matrix( prod( C, B_Operator ) ) );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//CALCULATE FORCEVECTORS WATER********************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::AddInternalForcesToRHSW( Vector& Help_R_W, const
        Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double Weight, double  DetJ )
{

    unsigned int pressure_size = GetGeometry().size();

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double porosity = GetPorosity( DN_DX_DISP );

    double DS_Dpc = GetDerivativeDSaturationDpc( capillaryPressure );

    double saturation = GetSaturation( capillaryPressure );

    double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

    double div_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

    Vector flow_water( dim );
    noalias( flow_water ) = GetFlowWater( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

    for ( unsigned int prim = 0; prim < pressure_size; ++prim )
    {
        Help_R_W( prim ) -=
            N_PRESS( prim ) * porosity * DS_Dpc * Dpc_Dt * Weight * DetJ/* * GetProperties()[SCALE]*/
          + N_PRESS( prim ) * saturation * div_Dt * Weight * DetJ/* * GetProperties()[SCALE]*/ ;

        for ( unsigned int gamma = 0; gamma < dim; ++gamma )
        {
            Help_R_W( prim ) +=
                ( DN_DX_PRESS( prim, gamma ) * flow_water( gamma ) )
                * Weight * DetJ/* * GetProperties()[SCALE]*/ ;
        }
    }
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateStiffnesMatrixWW( Matrix& Help_K_WW,
        const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = GetGeometry().size();

    double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

    double porosity = GetPorosity( DN_DX_DISP );

    double D2S_Dpc2 = GetSecondDerivativeD2SaturationDpc2( capillaryPressure );

    double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

    Vector Dflow_waterDpw( dim );

    noalias( Dflow_waterDpw ) = GetDerivativeDWaterFlowDpw( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

    double Dflow_waterDgradpw = GetDerivativeDWaterFlowDGradpw( DN_DX_DISP, capillaryPressure );

    double Ddiv_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

    for ( unsigned int prim = 0; prim < pressure_size; ++prim )
    {
        for ( unsigned int sec = 0; sec < pressure_size; ++sec )
        {
            Help_K_WW( prim, sec ) -=
                N_PRESS( prim ) * porosity * D2S_Dpc2 * Dpc_Dt * N_PRESS( sec )
                * Weight * DetJ/* * GetProperties()[SCALE]*/;

            Help_K_WW( prim, sec ) -=
                N_PRESS( prim ) * DSDpc * Ddiv_Dt * N_PRESS( sec )
                * Weight * DetJ/* * GetProperties()[SCALE]*/;

            for ( unsigned int gamma = 0; gamma < dim; gamma++ )
            {
                Help_K_WW( prim, sec ) -=
                    DN_DX_PRESS( prim, gamma ) * Dflow_waterDpw( gamma )
                    * N_PRESS( sec )
                    * Weight * DetJ/* * GetProperties()[SCALE]*/;
                Help_K_WW( prim, sec ) -=
                    DN_DX_PRESS( prim, gamma ) * Dflow_waterDgradpw
                    * DN_DX_PRESS( sec, gamma ) * Weight * DetJ/* * GetProperties()[SCALE]*/;
            }
        }
    }
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateDampingMatrixWW( Matrix& Help_D_WW, const Matrix& DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
{
    unsigned int pressure_size = GetGeometry().size();
    double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
    double porosity = GetPorosity( DN_DX_DISP );

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < pressure_size; sec++ )
        {
            Help_D_WW( prim, sec ) -=
                N_PRESS( prim ) * porosity * DSDpc * N_PRESS( sec )
                * Weight * DetJ /* * GetProperties()[SCALE]*/;
        }
    }
}

//************************************************************************************
//************************************************************************************

Vector UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetGradientWater( const Matrix& DN_DX_PRESS )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    Vector result( dim );
    noalias( result ) = ZeroVector( dim );

    double presW_alpha;

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
    {
        presW_alpha = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );

        for ( unsigned int k = 0; k < 3; ++k )
        {
            result( k ) += presW_alpha * DN_DX_PRESS(( i ), k );
        }
    }

    return result;
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetPressures( const Vector& N_PRESS,
        double& capillaryPressure, double& waterPressure )
{
    capillaryPressure = 0.0;

    waterPressure = 0.0;

    double presW_alpha;

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
    {
        presW_alpha = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );

        capillaryPressure += ( -presW_alpha ) * N_PRESS( i );

        waterPressure += presW_alpha * N_PRESS( i );
    }
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetDerivativeDPressuresDt( const Vector& N_PRESS,
        double& capillaryPressure_Dt, double& waterPressure_Dt )
{
    capillaryPressure_Dt = 0.0;

    waterPressure_Dt = 0.0;

    double presW_alpha_Dt;

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
    {
        presW_alpha_Dt = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT );

        capillaryPressure_Dt += ( -presW_alpha_Dt ) * N_PRESS( i );

        waterPressure_Dt += presW_alpha_Dt * N_PRESS( i );
    }
}

//************************************************************************************
//************************************************************************************

double UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetDerivativeDCapillaryPressureDt( const Vector& N_PRESS )
{
    double capillaryPressure_Dt = 0.0;

    double presW_alpha_Dt;

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
    {
        presW_alpha_Dt = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT );

        capillaryPressure_Dt += ( -presW_alpha_Dt ) * N_PRESS( i );
    }

    return capillaryPressure_Dt;
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

//POROSITY AND ITS DERIVATIVES
//TODO check the implementation of GetPorosity, it should depend on the gradient of displacement
double UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetPorosity( const Matrix& DN_DX_DISP )
{
    //TODO implement the divergent porosity
    double porosity = GetValue(POROSITY);
    return porosity;
}

double UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetDerivativeDPorosityDDivU( const Matrix& DN_DX_DISP )
{
    //TODO implement the divergent porosity
    double porosity_divu = 0.0;
    return porosity_divu;
}


//************************************************************************************
//************************************************************************************
double UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetDivU( const Matrix& DN_DX_DISP )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double div = 0.0;

    Vector u_alpha( dim );

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
    {
        noalias( u_alpha ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT );

        for ( unsigned int k = 0; k < dim; ++k )
        {
            div += ( u_alpha( k ) ) * DN_DX_DISP( i, k );
        }
    }

    return div;
}

//************************************************************************************
//************************************************************************************
double UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetDerivativeDDivUDt( const Matrix& DN_DX_DISP )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double div = 0.0;

    Vector u_alpha_Dt( dim );

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
    {
        noalias( u_alpha_Dt ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT );

        for ( unsigned int k = 0; k < dim; ++k )
        {
            div += u_alpha_Dt( k ) * DN_DX_DISP( i, k );
        }
    }

    return div;
}

//AVERAGED DENSITY
//************************************************************************************
//************************************************************************************

double UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetAveragedDensity( double capillaryPressure, double porosity )
{
    double result = 0.0;

    double density_soil = 0.0;
    if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
    {
        density_soil = GetValue(DENSITY);
    }
    else
    {
        density_soil = GetProperties()[DENSITY];
    }
    double density_air = GetValue(DENSITY_AIR);
    double density_water = GetValue(DENSITY_WATER);
    double saturation = GetSaturation( capillaryPressure );

    result = ( 1 - porosity ) * density_soil +
             porosity * ( saturation * density_water + ( 1 - saturation ) * density_air );

    return result;
}


//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

//SATURATION AND ITS DERIVATIVES

double UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetSaturation( double capillaryPressure )
{
    double airEntryPressure = GetValue(AIR_ENTRY_VALUE);

    if ( airEntryPressure <= 0.0 )
        airEntryPressure = 1.0;

    double b = GetValue(FIRST_SATURATION_PARAM);

    double c = GetValue(SECOND_SATURATION_PARAM);

    double saturation = 0.0;

    if ( capillaryPressure < 0.0 )
        capillaryPressure = 0.0;

    saturation = pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c ) );

// For Liakopolous Benchmark
// saturation =  1.0-1.9722*1e-11*pow(capillaryPressure,2.4279);

    return saturation;
}

//************************************************************************************
//************************************************************************************

double UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetDerivativeDSaturationDpc( double capillaryPressure )
{
    double airEntryPressure = GetValue(AIR_ENTRY_VALUE);

    if ( airEntryPressure <= 0 )
        airEntryPressure = 1.0;

    double b = GetValue(FIRST_SATURATION_PARAM);

    double c = GetValue(SECOND_SATURATION_PARAM);

    double result = 0.0;

    if ( capillaryPressure < 0.0 )
    {
        capillaryPressure = 0.0;
    }

    result = ( -c ) * pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 1.0 ) ) * b *
             pow(( capillaryPressure / airEntryPressure ), ( b - 1 ) ) * 1.0 / airEntryPressure;

// For Liakopolous Benchmark
// result =  -1.9722*2.4279*1e-11*pow(capillaryPressure,1.4279);

    return result;
}

//************************************************************************************
//************************************************************************************

double UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetSecondDerivativeD2SaturationDpc2( double capillaryPressure )
{
    double airEntryPressure = GetValue(AIR_ENTRY_VALUE);

    if ( airEntryPressure <= 0 )
        airEntryPressure = 1.0;

    double b = GetValue(FIRST_SATURATION_PARAM);

    double c = GetValue(SECOND_SATURATION_PARAM);

    double result = 0.0;

    if ( capillaryPressure < 0.0 )
        capillaryPressure = 0.0;

    result = ( -c ) * b / airEntryPressure * (
                 ( -c - 1.0 ) * pow(( 1 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 2.0 ) )
                 * b / airEntryPressure * pow(( capillaryPressure / airEntryPressure ), ( 2.0 * ( b - 1.0 ) ) )
                 + pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 1.0 ) ) * ( b - 1.0 )
                 * pow(( capillaryPressure / airEntryPressure ), ( b - 2.0 ) ) * 1.0 / airEntryPressure );

//    double aux1 = capillaryPressure / airEntryPressure;
//    double aux2 = 1.0 + pow(aux1, b );
//    result = ( -c ) * b / airEntryPressure * (
//                   ( -c - 1.0 ) * pow( aux2 , -c - 2.0 ) * b / airEntryPressure * pow(aux1, 2.0 * ( b - 1.0 ) )
//                 + pow(aux2, -c - 1.0 ) * ( b - 1.0 ) * pow(aux1, b - 2.0 ) / airEntryPressure
//                 );

// For Liakopolous Benschmark
// result =  -1.9722*2.4279*1.4279*1e-11*pow(capillaryPressure,0.4279);

    return result;
}

//************************************************************************************
//************************************************************************************

Vector UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetFlowWater( const Matrix& DN_DX_PRESS, const Matrix& DN_DX_DISP, double capillaryPressure )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //Calculation of relative Permeability after Mualem
    double relPerm = GetSaturation( capillaryPressure );

    if ( relPerm <= 0.01 )
        relPerm = 0.01;

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

    Vector gravity( dim );

    if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
        noalias( gravity ) = GetValue(GRAVITY);
    else
        noalias( gravity ) = GetProperties()[GRAVITY];

    Vector result( dim );

    noalias( result ) = ZeroVector( dim );

    Vector grad_water( dim );

    noalias( grad_water ) = GetGradientWater( DN_DX_PRESS );
    
    for ( unsigned int i = 0; i < dim; ++i )
    {
        result( i ) = -relPerm * GetValue(PERMEABILITY_WATER) /
                      ( GetValue(DENSITY_WATER) * 9.81 )
                      * ( grad_water( i ) - GetValue(DENSITY_WATER)
                          * gravity( i ) );
    }

    return result;
}

//************************************************************************************
//************************************************************************************

Vector UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetDerivativeDWaterFlowDpw( const Matrix& DN_DX_PRESS, const Matrix& DN_DX_DISP, double capillaryPressure )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
    //Calculation of Derivative relative Permeability after Mualem
    double relPerm = GetSaturation( capillaryPressure );
    double relPerm_pw = ( -1 ) * DSDpc;

    if ( relPerm <= 0.0001 )
    {
        relPerm = 0.0001;
        relPerm_pw = 0.0;
    }

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm_pw= -2.207*1.0121*pow((1-saturation),0.0121)*(-DSDpc)*(-1);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

    Vector result( dim );

    noalias( result ) = ZeroVector( dim );

    Vector gravity( dim );

    if( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
        noalias( gravity ) = GetValue(GRAVITY);
    else
        noalias( gravity ) = GetProperties()[GRAVITY];

    Vector grad_water( dim );

    noalias( grad_water ) = GetGradientWater( DN_DX_PRESS );

    for ( unsigned int i = 0; i < dim; ++i )
    {
        result( i ) = -relPerm_pw * GetValue(PERMEABILITY_WATER) /
                      ( GetValue(DENSITY_WATER) * 9.81 )
                      * ( grad_water( i ) - GetValue(DENSITY_WATER)
                          * gravity( i ) );
    }

    return result;
}

//************************************************************************************
//************************************************************************************

double UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetDerivativeDWaterFlowDGradpw( const Matrix& DN_DX_DISP, double capillaryPressure )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //Calculation of Derivative relative Permeability after Mualem
    double relPerm = GetSaturation( capillaryPressure );

    if ( relPerm <= 0.0001 )
        relPerm = 0.0001;

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

    double result( dim );

    Vector gravity( dim );

    if( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
        noalias( gravity ) = GetValue(GRAVITY);
    else
        noalias( gravity ) = GetProperties()[GRAVITY];

    result = -relPerm * GetValue(PERMEABILITY_WATER) / ( GetValue(DENSITY_WATER) * norm_2(gravity) );

    return result;
}

/**
 * Computes the strain vector
 */
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector )
{
    KRATOS_TRY
    noalias( StrainVector ) = ZeroVector( 6 );

    for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
    {
        for ( unsigned int item = 0; item < 6; item++ )
            for ( unsigned int dim = 0; dim < 3; dim++ )
                StrainVector[item] += B( item, 3 * node + dim ) * ( Displacements( node, dim ) - mInitialDisp( node, dim ) );
    }

    KRATOS_CATCH( "" )
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX )
{
    unsigned int number_of_nodes_disp = GetGeometry().size();

    noalias( B_Operator ) = ZeroMatrix( 6, number_of_nodes_disp * 3 );

    for ( unsigned int i = 0; i < number_of_nodes_disp; ++i )
    {
        B_Operator( 0, i*3 )     = DN_DX( i, 0 );
        B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
        B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 );
        B_Operator( 3, i*3 )     = DN_DX( i, 1 );
        B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
        B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
        B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
        B_Operator( 5, i*3 )     = DN_DX( i, 2 );
        B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::CalculateBBaroperator( Matrix& B_Operator, const Matrix& DN_DX, const Matrix& Bdil_bar )
{
    unsigned int number_of_nodes_disp = GetGeometry().size();
    
    noalias( B_Operator ) = ZeroMatrix( 6, number_of_nodes_disp * 3 );

    double aux = (1.0 / 3);
    double tmp1;
    double tmp2;
    double tmp3;
    for ( unsigned int i = 0; i < number_of_nodes_disp; ++i )
    {
        tmp1 = aux * (Bdil_bar( i, 0 ) - DN_DX( i, 0 ));
        tmp2 = aux * (Bdil_bar( i, 1 ) - DN_DX( i, 1 ));
        tmp3 = aux * (Bdil_bar( i, 2 ) - DN_DX( i, 2 ));
        
        B_Operator( 0, i*3 ) = DN_DX( i, 0 ) + tmp1;
        B_Operator( 1, i*3 ) = tmp1;
        B_Operator( 2, i*3 ) = tmp1;

        B_Operator( 0, i*3 + 1)  = tmp2;
        B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 ) + tmp2;
        B_Operator( 2, i*3 + 1 ) = tmp2;
        
        B_Operator( 0, i*3 + 2)  = tmp3;
        B_Operator( 1, i*3 + 2 ) = tmp3;
        B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 ) + tmp3;
        
        B_Operator( 3, i*3 )     = DN_DX( i, 1 );
        B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
        B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
        B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
        B_Operator( 5, i*3 )     = DN_DX( i, 2 );
        B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
    }
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::InitializeMaterial
()
{
    KRATOS_TRY

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
    {
        mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[i]->SetValue( PARENT_ELEMENT_ID, this->Id(), *(ProcessInfo*)0);
        mConstitutiveLawVector[i]->SetValue( INTEGRATION_POINT_INDEX, i, *(ProcessInfo*)0);
        mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif

    KRATOS_CATCH( "" )
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
        rValues.resize( mConstitutiveLawVector.size() );

    if ( rVariable == ELASTIC_LEFT_CAUCHY_GREEN_OLD )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            if ( rValues[i].size1() != 3 || rValues[i].size2() != 3 )
                rValues[i].resize( 3, 3 );

            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( ELASTIC_LEFT_CAUCHY_GREEN_OLD, rValues[i] );
        }
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if ( rValues.size() != mConstitutiveLawVector.size() )
        rValues.resize( mConstitutiveLawVector.size() );

    if ( rVariable == MATERIAL_PARAMETERS )
    {
        if ( rValues.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
            rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
    }

    if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            if ( rValues[i].size() != 6 )
                rValues[i].resize( 6 );

            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( PRESTRESS, rValues[i] );
        }
    }
    
    if ( rVariable == PLASTIC_STRAIN_VECTOR )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            if ( rValues[i].size() != 6 )
                rValues[i].resize( 6 );

            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( PLASTIC_STRAIN_VECTOR, rValues[i] );
        }
    }

    //To Plot Internal variables
    if ( rVariable == INTERNAL_VARIABLES )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            if ( rValues[i].size() != 9 )
                rValues[i].resize( 9 );

            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( INTERNAL_VARIABLES, rValues[i] );
            
        }
    }

    //To Plot Stresses
    if ( rVariable == STRESSES )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            if ( rValues[i].size() != 6 )
                rValues[i].resize( 6 );

            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( STRESSES, rValues[i] );

        }
    }

    //To Plot Fluid Flows
    if ( rVariable == FLUID_FLOWS )
    {
        unsigned int number_of_nodes_press = GetGeometry().size();

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        Vector N( number_of_nodes_press );

        Matrix DN_DX( number_of_nodes_press, dim );

        double capillaryPressure;

        double waterPressure;

        double saturation;

        Vector waterFlow( 3 );

        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        const GeometryType::ShapeFunctionsGradientsType& DN_De =
            GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            if ( rValues[PointNumber].size() != 9 )
                rValues[PointNumber].resize( 9 );

            noalias( N ) = row( Ncontainer, PointNumber );

            GetPressures( N, capillaryPressure, waterPressure );

            saturation = GetSaturation( capillaryPressure );

            rValues[PointNumber][0] = waterPressure;

            rValues[PointNumber][1] = 0.0;

            rValues[PointNumber][2] = saturation;

            noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

            noalias( waterFlow ) = GetFlowWater( DN_DX, DN_DX, capillaryPressure );

            rValues[PointNumber][3] = waterFlow( 0 );

            rValues[PointNumber][4] = waterFlow( 1 );

            rValues[PointNumber][5] = waterFlow( 2 );

            rValues[PointNumber][6] = 0.0;

            rValues[PointNumber][7] = 0.0;

            rValues[PointNumber][8] = 0.0;
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Clean();
        #endif
    }

    //To Plot Coordinates of Integration Points
    if ( rVariable == COORDINATES )
    {
        for ( unsigned int i = 0; i < integration_points.size(); i++ )
        {
            if ( rValues[i].size() != 3 )
                rValues[i].resize( 3 );

            Geometry<Node<3> >::CoordinatesArrayType dummy;

            GetGeometry().GlobalCoordinates( dummy, integration_points[i] );

            noalias( rValues[i] ) = dummy;
        }
    }

    if ( rVariable == STRAIN || rVariable == CURRENT_STRAIN_VECTOR )
    {
        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        // calculate shape function values and local gradients
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int strain_size = dim * (dim + 1) / 2;
        unsigned int mat_size = dim * number_of_nodes;
        Matrix B(strain_size, mat_size);
        Vector StrainVector(strain_size);
        Matrix DN_DX(number_of_nodes, dim);
        Matrix CurrentDisp(number_of_nodes, dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
//        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        // extract current displacements
        for (unsigned int node = 0; node < GetGeometry().size(); ++node)
            noalias(row(CurrentDisp, node)) =
                GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        {
            if (rValues[i].size() != 6)
                rValues[i].resize(6);

            // compute B_Operator at the current integration point
            noalias(DN_DX) = prod(DN_De[i], mInvJ0[i]);
            CalculateBoperator(B, DN_DX);

            // compute the strain at integration point
            CalculateStrain(B, CurrentDisp, StrainVector);
            if(dim == 2)
            {
                rValues[i](0) = StrainVector(0);
                rValues[i](1) = StrainVector(1);
                rValues[i](2) = 0.0; // note: it's only correct for plane strain, TODO: we must make this available for plane stress constitutive law
                rValues[i](3) = StrainVector(2);
                rValues[i](4) = 0.0;
                rValues[i](5) = 0.0;
            }
            else if(dim == 3)
                noalias(rValues[i]) = StrainVector;
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Clean();
        #endif
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    if( rVariable == PLASTICITY_INDICATOR )
    {
        if ( rValues.size() != GetGeometry().IntegrationPoints().size() )
            rValues.resize( GetGeometry().IntegrationPoints().size() );

        //reading integration points and local gradients
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
        {
            rValues[Point] = mConstitutiveLawVector[Point]->GetValue( rVariable, rValues[Point] );
        }
        return;
    }
    
    if( rVariable == K0 )
    {
        if ( rValues.size() != GetGeometry().IntegrationPoints().size() )
            rValues.resize( GetGeometry().IntegrationPoints().size() );
        
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
        {
            rValues[Point] = GetValue( K0 );
        }
        
        return;
    }

    unsigned int number_of_nodes_press = GetGeometry().size();

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if ( rValues.size() != integration_points.size() )
        rValues.resize( integration_points.size() );

    const Matrix& Ncontainer_Pressure = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

    Vector N_PRESS( number_of_nodes_press );

    double capillaryPressure;

    double waterPressure;

    double saturation;


    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        GeometryType::CoordinatesArrayType gp_position;
        gp_position = GetGeometry().GlobalCoordinates( gp_position, GetGeometry().IntegrationPoints()[PointNumber] );

        // Shape Functions on current spatial quadrature point
        if ( N_PRESS.size() != number_of_nodes_press )
            N_PRESS.resize( number_of_nodes_press );

        noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

        GetPressures( N_PRESS, capillaryPressure, waterPressure );

        saturation = GetSaturation( capillaryPressure );

        if ( rVariable == SATURATION )
        {
            rValues[PointNumber] = saturation;
        }

        if ( rVariable == WATER_PRESSURE )
        {
            rValues[PointNumber] = waterPressure;
//            KRATOS_WATCH(waterPressure)
        }

        if ( rVariable == EXCESS_PORE_WATER_PRESSURE )
        {
            rValues[PointNumber] = waterPressure - mReferencePressures[PointNumber];
        }
    }
    /////////////////////////////////////////////////////////////////////////
    //// End Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif

    KRATOS_CATCH( "" )
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
    {
        std::cout << "wrong size: " << rValues.size() << "!=" << mConstitutiveLawVector.size() << std::endl;
        return;
    }

    if ( rVariable == ELASTIC_LEFT_CAUCHY_GREEN_OLD )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i]->SetValue( ELASTIC_LEFT_CAUCHY_GREEN_OLD, rValues[i], rCurrentProcessInfo );
        }
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
        return;

    if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->SetValue( PRESTRESS, rValues[i], rCurrentProcessInfo );

    if ( rVariable == INTERNAL_VARIABLES )
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->SetValue( INTERNAL_VARIABLES, rValues[i], rCurrentProcessInfo );

    if ( rVariable == MATERIAL_PARAMETERS )
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::SetValueOnIntegrationPoints( const Kratos::Variable< ConstitutiveLaw::Pointer >& rVariable, std::vector< ConstitutiveLaw::Pointer >& rValues, const Kratos::ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == CONSTITUTIVE_LAW )
    {
        for ( unsigned int i = 0; i < rValues.size(); i++ )
        {
            mConstitutiveLawVector[i] = rValues[i];
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::SetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if( rVariable == K0 )
    {
        SetValue( K0, rValues[0] );
    }
    
    else if( rVariable == REFERENCE_WATER_PRESSURE )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mReferencePressures[i] = rValues[i];
        }
    }
    else
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
        }
    }

}

UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::IntegrationMethod UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

int UnsaturatedSoilsElement_2phase_SmallStrain_Staggered::Check(const Kratos::ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (this->Id() < 1)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Element found with Id 0 or negative", __FUNCTION__);
    }

    if (mTotalDomainInitialSize < 0)
    {
        std::cout << "error on element -> " << this->Id() << std::endl;
        KRATOS_THROW_ERROR(std::logic_error, "Domain size can not be less than 0. Please check Jacobian.", __FUNCTION__);
    }

    //verify that the constitutive law exists
    if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id());
    }

    KRATOS_CATCH( "" )

    return 1;
}

} // Namespace Kratos

#undef ENABLE_DEBUG_CONSTITUTIVE_LAW
