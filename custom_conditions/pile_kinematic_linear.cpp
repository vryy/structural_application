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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jelena $
//   Date:                $Date: 2013-03-20 $
//   Revision:            $Revision: 1.2 $
//
//
// System includes

// External includes


// Project includes
#include "pile_kinematic_linear.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
Pile_Kinematic_Linear::Pile_Kinematic_Linear( IndexType NewId, GeometryType::Pointer pGeometry )
: Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************
Pile_Kinematic_Linear::Pile_Kinematic_Linear( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
: Condition( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

Pile_Kinematic_Linear::Pile_Kinematic_Linear( IndexType NewId,
                                              GeometryType::Pointer pGeometry,
                                              PropertiesType::Pointer pProperties,
                                              Element::Pointer soilElement,
                                              Element::Pointer pileElement,
                                              Point<3>& rSoilLocalPoint,
                                              Point<3>& rPileLocalPoint,
                                              int PileIntegrationPointIndex )
: Condition( NewId, pGeometry, pProperties )
{
    mSoilLocalPoint = rSoilLocalPoint;
    mPileLocalPoint = rPileLocalPoint;
    mpSoilElement = soilElement;
    mpPileElement = pileElement;

    mPileIntegrationPointIndex = PileIntegrationPointIndex;
    //Test for calculating coordinates at time step midpoint
    mSoilGlobalPoint = mpSoilElement->GetGeometry().GlobalCoordinates( mSoilGlobalPoint, mSoilLocalPoint );
    mPileGlobalPoint = mpPileElement->GetGeometry().GlobalCoordinates( mPileGlobalPoint, mPileLocalPoint );

    mvPile.resize( 3 );
    mTPile.resize( 2, 3 );
}

//********************************************************
//**** Operations ****************************************
//********************************************************

/**
* Destructor. Never to be called manually
*/
Pile_Kinematic_Linear::~Pile_Kinematic_Linear()
{
}

void Pile_Kinematic_Linear::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_CATCH( "" )
}

void Pile_Kinematic_Linear::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
}

Matrix Pile_Kinematic_Linear::TangentialVectors( Element::Pointer rElement,
        const GeometryType::CoordinatesArrayType& rPoint )
{
    //setting up result matrix
    Matrix T( 2, 3 );
    noalias( T ) = ZeroMatrix( 2, 3 );
    //shape function gradients
    Matrix DN = ZeroMatrix( rElement->GetGeometry().size(), 1 );

     rElement->GetGeometry().ShapeFunctionsLocalGradients( DN, rPoint );
    //calculating tangential vectors
    for ( unsigned int n = 0; n <  rElement->GetGeometry().PointsNumber(); n++ )
    {
        T( 0, 0 ) += ( rElement->GetGeometry().GetPoint( n ).X0()
                       + rElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                    * DN( n, 0 );
        T( 0, 1 ) += ( rElement->GetGeometry().GetPoint( n ).Y0()
                       + rElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                    * DN( n, 0 );
        T( 0, 2 ) += ( rElement->GetGeometry().GetPoint( n ).Z0()
                       + rElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                     * DN( n, 0 );
        T( 1, 0 ) += 0;
        T( 1, 1 ) += 0;
        T( 1, 2 ) += 0;
    }
    //KRATOS_WATCH (T);

    return( T );
}

/**
 * returns the tangential vectors of the current surface in an arbitrary point due to reference configuration
 * @param Surface Surface Condition for that the Tangential Vector should be calculated
 * @param rPoint local coordinates of the point in Surface for that the Tangential Vector should be calculated
 * @return Matrix(2,3) of the two tangential vectors
 */
Matrix Pile_Kinematic_Linear::TangentialVectors_inOrigin( Element::Pointer rElement,
        const GeometryType::CoordinatesArrayType& rPoint )
{
    //setting up result matrix
    Matrix T( 2, 3 );
    noalias( T ) = ZeroMatrix( 2, 3 );
    //shape function gradients
    Matrix DN = ZeroMatrix( rElement->GetGeometry().PointsNumber(), 1 );
    rElement->GetGeometry().ShapeFunctionsLocalGradients( DN, rPoint );
    //calculating tangential vectors

    for ( unsigned int n = 0; n < rElement->GetGeometry().PointsNumber(); n++ )
    {
        T( 0, 0 ) += rElement->GetGeometry().GetPoint( n ).X0() * DN( n, 0 );
        T( 0, 1 ) += rElement->GetGeometry().GetPoint( n ).Y0() * DN( n, 0 );
        T( 0, 2 ) += rElement->GetGeometry().GetPoint( n ).Z0() * DN( n, 0 );
        T( 1, 0 ) += rElement->GetGeometry().GetPoint( n ).X0() * DN( n, 1 );
        T( 1, 1 ) += rElement->GetGeometry().GetPoint( n ).Y0() * DN( n, 1 );
        T( 1, 2 ) += rElement->GetGeometry().GetPoint( n ).Z0() * DN( n, 1 );
    }

    return( T );
}

/**
 * returns the normalized tangential vectors of the current surface in an arbitrary point due to current configuration
 * @param Surface Surface Condition for that the Tangential Vector should be calculated
 * @param rPoint local coordinates of the point in Surface for that the Tangential Vector should be calculated
 * @return Matrix(2,3) of the two tangential vectors
 */
Matrix Pile_Kinematic_Linear::TangentialVectorsGlobal( Element::Pointer rElement,
        const GeometryType::CoordinatesArrayType& rPoint )
{
    //setting up result matrix
    Matrix T( 2, 3 );
    noalias( T ) = TangentialVectors( rElement, rPoint );

    //shape function gradients
    double normT1 = sqrt( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) );
    double normT2 = sqrt( T( 1, 0 ) * T( 1, 0 ) + T( 1, 1 ) * T( 1, 1 ) + T( 1, 2 ) * T( 1, 2 ) );
    T( 0, 0 ) = T( 0, 0 ) / normT1;
    T( 0, 1 ) = T( 0, 0 ) / normT1;
    T( 0, 2 ) = T( 0, 0 ) / normT1;
    T( 1, 0 ) = T( 0, 0 ) / normT2;
    T( 1, 1 ) = T( 0, 0 ) / normT2;
    T( 1, 2 ) = T( 0, 0 ) / normT2;
    return( T );
}

/**
 * returns the normal vector in arbitrary point
 * calculates the normalized vector orthogonal to the current surface in given point
 * @param rPoint the given point in local coordinates
 * @return the normal vector
 */
Vector Pile_Kinematic_Linear::NormalVector( Element::Pointer rElement,
        const GeometryType::CoordinatesArrayType& rPoint )
{
    Vector ShapeFunctionValues = ZeroVector( rElement->GetGeometry().PointsNumber() );

    for ( IndexType i = 0; i <rElement->GetGeometry().PointsNumber(); i++ )
        ShapeFunctionValues[i] = rElement->GetGeometry().ShapeFunctionValue( i, rPoint );

    Vector Result( 3 );
    noalias( Result ) = ZeroVector( 3 );

    //getting tangential vectors
    Matrix T( 2, 3 );
    noalias( T ) = TangentialVectors( rElement, rPoint );

    //calculating normal vector
    Vector vPC( 3 );
    for ( unsigned int n = 0; n < rElement->GetGeometry().PointsNumber(); n++ )
    {
        vPC[0] += ( rElement->GetGeometry()[n].X()
                             + mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_X,0 ) )
                         *  ShapeFunctionValues[n];
        vPC[1] += ( rElement->GetGeometry()[n].Y()
                                 + mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_Y, 0 ) )
                         *  ShapeFunctionValues[n];
        vPC[2] += ( rElement->GetGeometry()[n].Z()
                                + mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_Z, 0 ) )
                         *  ShapeFunctionValues[n];
    }

    Vector vP0( 3 );
    for ( unsigned int n = 0; n < rElement->GetGeometry().PointsNumber(); n++ )
    {
        vP0[0] += ( rElement->GetGeometry()[n].X())*  ShapeFunctionValues[n];
        vP0[1] += ( rElement->GetGeometry()[n].Y() )*  ShapeFunctionValues[n];
        vP0[2] += ( rElement->GetGeometry()[n].Z() )*  ShapeFunctionValues[n];
    }

    Vector vP0Proj( 3 );
    {
        double T0 = ( vP0(0)*T(0, 0) + vP0(1)*T(0,1) + vP0(2) * T(0,2) - ( T(0,0)*vPC[0] + T(0,1)* vPC[1] + T(0,2)* vPC[2])) / (T(0, 0) * T (0, 0) + T(0, 1) * T(0, 1) + T(0, 2) * T(0, 2) ) ;
        noalias( vP0Proj ) = ZeroVector( 3 );

        vP0Proj( 0 ) += ( vP0( 0 ) - T( 0, 0 ) * T0 );
        vP0Proj( 1 ) += ( vP0( 1 ) - T( 0, 1 ) * T0 );
        vP0Proj( 2 ) += ( vP0( 2 ) - T( 0, 2 ) * T0 );
    }

    Result[0] = vPC (0) - vP0Proj(0);
    Result[1] = vPC (0) - vP0Proj(0);
    Result[2] = vPC (0) - vP0Proj(0);
    //KRATOS_WATCH (Result);

    SD_MathUtils<double>::Normalize( Result );

    return( Result );
}


Matrix Pile_Kinematic_Linear::TangentialVectorsTotal( Element::Pointer rElement,
        const GeometryType::CoordinatesArrayType& rPoint )
{
    //calculating tangential vectors
    Matrix T( 2, 3 );
    noalias( T ) = TangentialVectors( rElement, rPoint );

    //calculating normal vectors
    Vector N( 1, 3 );
    noalias( N ) = NormalVector(  rElement, rPoint );

    T( 0, 0 ) = T( 0, 0 );
    T( 0, 1 ) = T( 0, 1 );
    T( 0, 2 ) = T( 0, 2 );

    T( 1, 0 ) = T( 0, 1 ) * N(2) - T( 0, 2 ) * N(1);
    T( 1, 1 ) = T( 0, 2 ) * N(0) - T( 0, 0 ) * N(2);
    T( 1, 2 ) = T( 0, 0 ) * N(1) - T( 0, 1 ) * N(0);

    return( T );
}

void Pile_Kinematic_Linear::CalculateRightHandSide( VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType matrix = Matrix();
    CalculateAll( matrix, rRightHandSideVector,
                  rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag,
                  CalculateResidualVectorFlag );
}

//************************************************************************************
//************************************************************************************
/**
 * calculates this contact element's local contributions
 */
void Pile_Kinematic_Linear::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

/**
 * This function calculates all system contributions due to the contact problem
 * with regard to the current Pile and Soil partners.
 * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node
 */
void Pile_Kinematic_Linear::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    //**********************************************
    //setting up the dimensions of the contributions
    //**********************************************
    unsigned int pileNN = mpPileElement->GetGeometry().size();
    unsigned int soilNN = mpSoilElement->GetGeometry().size();
    unsigned int dimension = 3;

    //resizing as needed the LHS
    int MatSize = ( pileNN + soilNN ) * dimension;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        rLeftHandSideMatrix.resize( MatSize, MatSize, false );
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        rRightHandSideVector.resize( MatSize, false );
        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    //*******************************
    //Calculating general information
    //*******************************

    //calculating shape function values for current soil element
    Vector SoilShapeFunctionValues = ZeroVector( soilNN );

    for ( IndexType i = 0; i < soilNN; i++ )
        SoilShapeFunctionValues[i] = mpSoilElement->GetGeometry().ShapeFunctionValue( i, mSoilLocalPoint );

    //calculating shape function gradients for current soil element
    Matrix SoilDN = ZeroMatrix( soilNN, dimension );

    noalias( SoilDN ) = mpSoilElement->GetGeometry().ShapeFunctionsLocalGradients( SoilDN, mSoilLocalPoint );

    //calculating shape function values for current pile element
    Vector PileShapeFunctionValues = ZeroVector( pileNN );

    for ( IndexType i = 0; i < pileNN; i++ )
        PileShapeFunctionValues[i] = mpPileElement->GetGeometry().ShapeFunctionValue( i, mPileLocalPoint );


    //calculating shape function gradients for current pile element
    Matrix PileDN = ZeroMatrix( pileNN, ( dimension - 2 ) ); ////??????????????????????????????????????????

    noalias( PileDN ) = mpPileElement->GetGeometry().ShapeFunctionsLocalGradients( PileDN, mPileLocalPoint );
      Matrix TSoil( 2, 3 );


    //********************************
    //Calculating contact area
    //********************************

    Matrix& inertia = mpPileElement->GetProperties()[INERTIA];
    double mArea = mpPileElement->GetProperties()[CROSS_AREA];
    double mInertia_x = inertia(0,0);
    double mInertia_y = inertia(1,1);
    double mInertia_Polar = inertia(2,2);

    double circumference_pile = sqrt ((12/mArea)*(mInertia_x*mInertia_y)+2*mArea) ;

    double DetJPile = mpPileElement->GetGeometry().DeterminantOfJacobian( mPileIntegrationPointIndex );

    double influence_area = mpPileElement->GetGeometry().IntegrationPoints()[mPileIntegrationPointIndex].Weight()
                            * DetJPile * circumference_pile*0.5;

    //KRATOS_WATCH( influence_area );
    for ( unsigned int n = 0; n < mpPileElement->GetGeometry().size(); n++ )
    {
        double DispX = mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_X, 0 );
    }


    //********************************
    //Calculating normals and tangents
    //********************************

    //calculating normal vector on Pile element

    Vector vPile( 3 );

    noalias( vPile ) = ZeroVector( 3 );

    //calculating tangential vectors
    Matrix TPileN( 2, 3 );

    noalias( TPileN ) = ZeroMatrix( 2, 3 );

    for ( unsigned int n = 0; n < mpPileElement->GetGeometry().size(); n++ )
    {
        TPileN( 0, 0 ) += ( mpPileElement->GetGeometry()[n].X())* PileDN( n, 0 );
        TPileN( 0, 1 ) += ( mpPileElement->GetGeometry()[n].Y())* PileDN( n, 0 );
        TPileN( 0, 2 ) += ( mpPileElement->GetGeometry()[n].Z() ) * PileDN( n, 0 );
        TPileN( 1, 0 ) +=0;
        TPileN( 1, 1 ) +=0;
        TPileN( 1, 2 ) +=0;
    }

    //calculating normal vetors vectors
    Vector vCurentPile( 3 );

   //calculating current position of target point on pile

    noalias( vCurentPile ) = ZeroVector( 3 );

    for ( unsigned int n = 0; n < pileNN; n++ )
    {
        vCurentPile[0] += ( mpPileElement->GetGeometry()[n].X() )*  PileShapeFunctionValues[n];
        vCurentPile[1] += ( mpPileElement->GetGeometry()[n].Y())*  PileShapeFunctionValues[n];
        vCurentPile[2] += ( mpPileElement->GetGeometry()[n].Z())*  PileShapeFunctionValues[n];

    }
   //calculating initial position of target point on pile

    Vector vOrigintPile = ZeroVector( 3 );

    for ( unsigned int n = 0; n < pileNN; n++ )
    {
        vOrigintPile[0] += ( mpPileElement->GetGeometry()[n].X()-mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_X,0 ) )
                           *  PileShapeFunctionValues[n];
        vOrigintPile[1] += ( mpPileElement->GetGeometry()[n].Y()-mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_Y,0 ) )
                           *  PileShapeFunctionValues[n];
        vOrigintPile[2] += ( mpPileElement->GetGeometry()[n].Z()-mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_Z,0 ))
                           *  PileShapeFunctionValues[n];

    }

    //calculating projection of curret point on nomal plane
    Vector vOriginPileProj( 3 );
    {
        double T0 = ( vOrigintPile( 0 ) * TPileN( 0, 0 ) +  vOrigintPile( 1 ) * TPileN( 0, 1 ) +  vOrigintPile( 2 ) * TPileN( 0, 2 ) - (TPileN( 0, 0 )* vCurentPile[0] + TPileN( 0, 1 )* vCurentPile[1] + TPileN( 0, 2 )* vCurentPile[2])   ) / ( TPileN( 0, 0 ) * TPileN( 0, 0 ) + TPileN( 0, 1 ) * TPileN( 0, 1 ) + TPileN( 0, 2 ) * TPileN( 0, 2 ) ) ;

        noalias( vOriginPileProj ) = ZeroVector( 3 );

        vOriginPileProj( 0 ) += ( vOrigintPile( 0 ) - TPileN( 0, 0 ) * T0 );
        vOriginPileProj( 1 ) += ( vOrigintPile( 1 ) - TPileN( 0, 1 ) * T0 );
        vOriginPileProj( 2 ) += ( vOrigintPile( 2 ) - TPileN( 0, 2 ) * T0 );

    }

    //calculating normal vector
    vPile = vCurentPile - vOriginPileProj;

    noalias( mPileGlobalPoint ) = GetGlobalCoordinates( mpPileElement, mPileGlobalPoint, mPileLocalPoint );

    noalias( mSoilGlobalPoint ) = GetGlobalCoordinates( mpSoilElement, mSoilGlobalPoint, mSoilLocalPoint );
    double norm_vPile;
    norm_vPile = sqrt( vPile( 0 ) * vPile( 0 ) + vPile( 1 ) * vPile( 1 ) + vPile(2 ) * vPile( 2 ) );
//    KRATOS_WATCH( vPile );
    //KRATOS_WATCH (norm_vPile);
     Vector relDisp( 3 );

    noalias( relDisp ) = ( mPileGlobalPoint - mSoilGlobalPoint );
    double NormrelDisp = MathUtils<double>::Norm3( relDisp );

    SD_MathUtils<double>::Normalize( vPile );

    //calculating cecond tangen as pordact of n*t1
    if (norm_vPile !=0)
    {
        TPileN( 1, 0 ) = TPileN( 0, 1 ) * vPile[2] - TPileN( 0, 2 ) * vPile[1];

        TPileN( 1, 1 ) = TPileN( 0, 2 ) * vPile[0] - TPileN( 0, 0 ) * vPile[2];

        TPileN( 1, 2 ) = TPileN( 0, 0 ) * vPile[1] - TPileN( 0, 1 ) * vPile[0];
    }
    else
    {
 //     noalias( TPileN ) = TangentialVectors( mpPileElement, mPileLocalPoint);
        TPileN( 1, 0 ) ==0;
        TPileN( 1, 1 ) ==0;
        TPileN( 1, 2 ) ==0;
    }

    //KRATOS_WATCH( TPileN );
    double Gap;
    if (norm_vPile !=0)
    {
        Gap = inner_prod( vPile, relDisp );
    }
    else
    {
        Gap == 0;
    }

    double penalty =1.0e+10;
    double penalty_T =1.0e+10;
    double friction_coeff = 0.2;

    //calculating normal contact stress
    double normalStress = abs( penalty * Gap );
   // double normalStress = penalty * Gap;

    //KRATOS_WATCH( normalStress );

    //tangential vector on Pile surface
    Matrix T( 2, 3 );

    noalias( T ) = ZeroMatrix( 2, 3 );

    //calculating tangential vectors

    for ( unsigned int n = 0; n < mpPileElement->GetGeometry().size(); n++ )
    {
        T( 0, 0 ) += ( mpPileElement->GetGeometry()[n].X()) * PileDN( n, 0 );
        T( 0, 1 ) += ( mpPileElement->GetGeometry()[n].Y())* PileDN( n, 0 );
        T( 0, 2 ) += ( mpPileElement->GetGeometry()[n].Z())* PileDN( n, 0 );
    }

    Vector vCurentPileB( 3 );

    noalias( vCurentPileB ) = ZeroVector( 3 );

    for ( unsigned int n = 0; n < pileNN; n++ )
    {
        vCurentPileB( 0 ) += ( mpPileElement->GetGeometry()[n].X())*  PileShapeFunctionValues[n];
        vCurentPileB( 1 ) += ( mpPileElement->GetGeometry()[n].Y())*  PileShapeFunctionValues[n];
        vCurentPileB( 2 ) += ( mpPileElement->GetGeometry()[n].Z() )*  PileShapeFunctionValues[n];
    }

    Vector vOrigintPileB( 3 );

    noalias( vOrigintPileB ) = ZeroVector( 3 );

    for ( unsigned int n = 0; n < pileNN; n++ )
    {
        vOrigintPileB( 0 ) += (( mpPileElement->GetGeometry()[n].X() - mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_X, 0 ) )
                               *  PileShapeFunctionValues[n]);
        vOrigintPileB( 1 ) += (( mpPileElement->GetGeometry()[n].Y() -mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_Y, 0 ) )
                               *  PileShapeFunctionValues[n]);
        vOrigintPileB( 2 ) += (( mpPileElement->GetGeometry()[n].Z() -mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_Z, 0 ))
                               *  PileShapeFunctionValues[n]);
    }

    Vector vOriginPileProjB( 3 );
    {
        double T0B = ( vOrigintPileB( 0 ) * T( 0, 0 ) +  vOrigintPileB( 1 ) * T( 0, 1 ) +  vOrigintPileB( 2 ) * T( 0, 2 ) - (T( 0, 0 )* vCurentPileB[0] + T( 0, 1 )* vCurentPileB[1] + T( 0, 2 )* vCurentPileB[2])   ) / ( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) ) ;

        double T0BX = ( vOrigintPileB( 0 ) * T( 0, 0 ) +  vOrigintPileB( 1 ) * T( 0, 1 ) +  vOrigintPileB( 2 ) * T( 0, 2 )   ) / ( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) ) ;

        noalias( vOriginPileProjB ) = ZeroVector( 3 );

        vOriginPileProjB[0] = ( vOrigintPileB[0] - T( 0, 0 ) * T0B );
        vOriginPileProjB[1] = ( vOrigintPileB[1] - T( 0, 1 ) * T0B  );
        vOriginPileProjB[2] = ( vOrigintPileB[2] - T( 0, 2 ) * T0B );
    }

    Vector vPileNonNormalized = vCurentPileB - vOriginPileProjB;

    Vector norm_T( 2 );

    norm_T[0] = sqrt( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) );
    norm_T[1] = sqrt( T( 1, 0 ) * T( 1, 0 ) + T( 1, 1 ) * T( 1, 1 ) + T( 1, 2 ) * T( 1, 2 ) );

    double NormvPile = MathUtils<double>::Norm3( vPileNonNormalized );
    if (NormvPile !=0)
    {
      T( 1, 0 ) = T( 0, 1 ) * vPileNonNormalized[2] - T( 0, 2 ) * vPileNonNormalized[1];
      T( 1, 1 ) = T( 0, 2 ) * vPileNonNormalized[0] - T( 0, 0 ) * vPileNonNormalized[2];
      T( 1, 2 ) = T( 0, 0 ) * vPileNonNormalized[1] - T( 0, 1 ) * vPileNonNormalized[0];
    }
    else
    {
        T( 1, 0 ) ==0;
        T( 1, 1 ) ==0;
        T( 1, 2 ) ==0;
    }

    Matrix m( 2, 2 );
    noalias( m ) = ZeroMatrix( 2, 2 );//m[alpha][beta]=T[alpha]*T[beta]

    for ( int i = 0; i < 2; i++ )
    {
        for ( int j = 0; j < 2; j++ )
        {
            m( i, j ) = 1.0;
        }
    }

    Vector tangentialStresses = ZeroVector( 2 );

    Vector tangentialStresses_trial = ZeroVector( 2 );

    Vector relativTangentialVelocity = ZeroVector( 2 );

    Vector relativVelocity = ZeroVector( 3 );

    double normTangentialStresses_trial = 0.0;

    bool Stick = true;

    if ( friction_coeff > 0.0 )
    {
        //TEST Calculation of relative velocity between Pile and Soil surface
//  noalias(relativTangentialVelocity)= GetRelativTangentialVelocity(T);
        Vector relativTangentialVelocity( 2 );
        noalias( relativTangentialVelocity ) = ZeroVector( 2 );
        if (NormrelDisp != 0 && (NormrelDisp != Gap))
        {
            Vector soil_velo = ZeroVector( 3 );
            Vector pile_velo = ZeroVector( 3 );

            for ( IndexType i = 0 ; i < mpSoilElement->GetGeometry().size() ; i++ )
            {
                double shape_func = mpSoilElement->GetGeometry().ShapeFunctionValue( i, mSoilLocalPoint );
                soil_velo += (( mpSoilElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) ) *
                             shape_func;
            }

            for ( IndexType i = 0 ; i < mpPileElement->GetGeometry().size() ; i++ )
            {
                double shape_func = mpPileElement->GetGeometry().ShapeFunctionValue( i, mPileLocalPoint );
                pile_velo += (( mpPileElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) ) *
                             shape_func;
            }

            Vector norm_T( 2 );

            norm_T( 0 ) = sqrt( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) );

            norm_T( 1 ) = sqrt( T( 1, 0 ) * T( 1, 0 ) + T( 1, 1 ) * T( 1, 1 ) + T( 1, 2 ) * T( 1, 2 ) );

        if ((norm_T (0)) != 0)
        {
          relativTangentialVelocity( 0 ) = (( soil_velo( 0 ) - pile_velo( 0 ) ) * T( 0, 0 ) + ( soil_velo( 1 ) - pile_velo( 1 ) ) * T( 0, 1 )
                                              + ( soil_velo( 2 ) - pile_velo( 2 ) ) * T( 0, 2 ) ) / norm_T( 0 );
        } else
        {
         relativTangentialVelocity( 0 ) = 0;
        }

        if ((norm_T (1)) != 0)
        {
        relativTangentialVelocity( 1 ) = (( soil_velo( 0 ) - pile_velo( 0 ) ) * T( 1, 0 ) + ( soil_velo( 1 ) - pile_velo( 1 ) ) * T( 1, 1 )
                                              + ( soil_velo( 2 ) - pile_velo( 2 ) ) * T( 1, 2 ) ) / norm_T( 1 );
        } else
        {
          relativTangentialVelocity( 1 ) = 0;
        }






        tangentialStresses_trial[0] +=  relativTangentialVelocity( 0 ) * penalty_T;

        tangentialStresses_trial[1] +=  relativTangentialVelocity( 1 ) * penalty_T;



        normTangentialStresses_trial =  sqrt( tangentialStresses_trial[0] * m( 0, 0 ) * tangentialStresses_trial[0]
                                              + tangentialStresses_trial[0] * m( 0, 1 ) * tangentialStresses_trial[1]
                                              + tangentialStresses_trial[1] * m( 1, 0 ) * tangentialStresses_trial[0]
                                              + tangentialStresses_trial[1] * m( 1, 1 ) * tangentialStresses_trial[1] );




      if ( normTangentialStresses_trial > friction_coeff*normalStress )//Slip
      {
            tangentialStresses[0] = friction_coeff * normalStress * tangentialStresses_trial[0] / normTangentialStresses_trial;//
            tangentialStresses[1] = friction_coeff * normalStress * tangentialStresses_trial[1] / normTangentialStresses_trial;//
            Stick = false;


      }
      else //Stick
      {

            tangentialStresses[0] = tangentialStresses_trial[0];
            tangentialStresses[1] = tangentialStresses_trial[1];
            Stick = true;
      }

    }
    else
      return;
    }



    //Calculating Soil element's current integration weight
    double SoilIntegrationWeight = mpSoilElement->GetGeometry().IntegrationPoints()[mPileIntegrationPointIndex].Weight();


    if ( CalculateResidualVectorFlag == true )
    {
        if ( normalStress > 0.0 )
        {
            CalculateAndAdd_RHS( rRightHandSideVector,
                                 PileShapeFunctionValues,
                                 SoilShapeFunctionValues,
                                 vPile,
                                 T,
                                 tangentialStresses,
                                 Gap,
                                 normalStress,
                                 SoilIntegrationWeight,
                                 influence_area,
                                 friction_coeff );
        }
    }

    //END OF ADDING RESIDUAL CONTRIBUTIONS
    //************************************

    if ( CalculateStiffnessMatrixFlag == false )
    {
        return;
    }

    //*******************************************************************************
    //****************** adding all contributions ***********************************
    //*******************************************************************************

    //*************************************************************
    //BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF CONTACT STRESS
    //*************************************************************

    //BEGIN OF CONTRIBUTION DUE TO DISPLACEMENTS ON CONTACT SURFACES

    //BEGIN OF CONTRIBUTION: MASTER-MASTER


    if (Gap != 0)
    {
          Vector cont1 (2);
    Vector cont2 (2);
    Vector cont3 (2);
    Vector cont4 (2);
    noalias (cont1) = ZeroVector (2);
        for ( unsigned int prim = 0; prim < pileNN ;  prim++ )
        {
            for ( unsigned int i = 0; i < dimension; i++ )
            {
                for ( unsigned int sec = 0; sec < pileNN; sec++ )
                {
                    for ( unsigned int j = 0; j < dimension; j++ )
                    {
                        rLeftHandSideMatrix( prim*dimension + i, sec*dimension + j )
                        += PileShapeFunctionValues[prim] * vPile[i]
                           * PileShapeFunctionValues[sec] * vPile[j]
                           * penalty * SoilIntegrationWeight * influence_area;

                cont1 (0) =  prim*dimension + i;
                cont1 (1) = sec*dimension + j;
            //    KRATOS_WATCH (cont1)
                    }
                }
            }
        }



        //END OF CONTRIBUTION: MASTER-MASTER
        //BEGIN OF CONTRIBUTION: MASTER-SLAVE
        for ( unsigned int prim = 0; prim < pileNN; prim++ )
        {
            for ( unsigned int i = 0; i < dimension; i++ )
            {
                for ( unsigned int sec = 0; sec < soilNN; sec++ )
                {
                    for ( unsigned int j = 0; j < dimension; j++ )
                    {
                        rLeftHandSideMatrix( prim*dimension + i, sec*dimension + j + pileNN*dimension )
                        -= PileShapeFunctionValues[prim] * vPile[i]
                           * SoilShapeFunctionValues[sec] * vPile[j] * penalty
                           * SoilIntegrationWeight * influence_area;

            //   //KRATOS_WATCH (cont2);

                    }
                }
            }
        }

        //END OF CONTRIBUTION: MASTER-SLAVE
        //BEGIN OF CONTRIBUTION: SLAVE-MASTER
        for ( unsigned int prim = 0; prim < soilNN; prim++ )
        {
            for ( unsigned int i = 0; i < dimension; i++ )
            {
                for ( unsigned int sec = 0; sec < pileNN; sec++ )
                {
                    for ( unsigned int j = 0; j < dimension; j++ )
                    {
                        rLeftHandSideMatrix( prim*dimension + i + pileNN*dimension, sec*dimension + j )
                        -= SoilShapeFunctionValues[prim] * vPile[i]
                           * PileShapeFunctionValues[sec] * vPile[j] * penalty
                           * SoilIntegrationWeight * influence_area;

                    }
                }
            }
//       //KRATOS_WATCH (rLeftHandSideMatrix);
        }

        //END OF CONTRIBUTION: SLAVE-MASTER
        //BEGIN OF CONTRIBUTION: SLAVE-SLAVE
        for ( unsigned int prim = 0; prim < soilNN; prim++ )
        {
            for ( unsigned int i = 0; i < dimension; i++ )
            {
                for ( unsigned int sec = 0;  sec < soilNN;  sec++ )
                {
                    for ( unsigned int j = 0; j < dimension; j++ )
                    {
                        rLeftHandSideMatrix( prim*dimension + i + pileNN*dimension, sec*dimension + j + pileNN*dimension )
                        += SoilShapeFunctionValues[prim] * vPile[i]
                           * SoilShapeFunctionValues[sec] * vPile[j] * penalty
                           * SoilIntegrationWeight * influence_area;

                    }
                }
            }
        }
    }

    else
        return;

    //END OF CONTRIBUTION: SLAVE-SLAVE
    //END OF CONTRIBUTION DUE TO DISPLACEMENTS ON CONTACT SURFACES
    if ( friction_coeff < 0.0 )
        return;

    if ( Stick )  return;

//             //*****************************************************************
//             //BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF FRICTIONAL STRESS
//              // ONLY IF STICK IS FALSE
//             //*****************************************************************
//BEGIN OF CONTRIBUTION MASTER-MASTER
    if (normTangentialStresses_trial != 0)
    {
      for ( unsigned int prim = 0; prim < pileNN; prim++ )
      {
        for ( unsigned int i = 0; i < dimension; i++ )
        {
            Vector XiPrim = ZeroVector( 2 );

        if ((norm_T (0)) != 0)
        {
          XiPrim[0] = -PileShapeFunctionValues[prim] * T( 0, i ) / norm_T( 0 );
        }
        else
        {
         XiPrim[0]=0;
        }

        if ((norm_T (1)) != 0)
        {
          XiPrim[1] = -PileShapeFunctionValues[prim] * T( 1, i ) / norm_T( 1 );
        }
        else
        {
          XiPrim[1]= 0;
        }

            for ( unsigned int sec = 0; sec < pileNN; sec++ )
            {
                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    double normalStress_DU = PileShapeFunctionValues[sec]
                                             * vPile[j] * penalty;

                    Vector tangentialStresses_DU( 2 );
            if (normTangentialStresses_trial != 0)
            {
              tangentialStresses_DU( 0 ) = friction_coeff
                                                 * normalStress_DU * tangentialStresses_trial[0] / normTangentialStresses_trial;

              tangentialStresses_DU( 1 ) = friction_coeff
                                                 * normalStress_DU * tangentialStresses_trial[1] / normTangentialStresses_trial;
            } else
            {
              tangentialStresses_DU( 0 ) = 0;
              tangentialStresses_DU( 1 ) = 0;
            }

                    rLeftHandSideMatrix( prim*dimension + i, sec*dimension + j) +=
                        ( tangentialStresses_DU[0] * XiPrim[0] +
                          tangentialStresses_DU[1] * XiPrim[1] )
                        * SoilIntegrationWeight * influence_area;
                }
            }
        }
      }


    //END OF CONTRIBUTION: MASTER-MASTER

    //BEGIN OF CONTRIBUTION MASTER-SLAVE
    for ( unsigned int prim = 0; prim < pileNN; prim++ )
    {
        for ( unsigned int i = 0; i < dimension; i++ )
        {
            Vector XiPrim = ZeroVector( 2 );

        if ((norm_T (0)) != 0)
        {
          XiPrim[0] = -PileShapeFunctionValues[prim] * T( 0, i ) / norm_T( 0 );
        }
        else
        {
         XiPrim[0]=0;
        }

        if ((norm_T (1)) != 0)
        {
          XiPrim[1] = -PileShapeFunctionValues[prim] * T( 1, i ) / norm_T( 1 );
        }
        else
        {
          XiPrim[1]= 0;
        }

            for ( unsigned int sec = 0; sec < soilNN; sec++ )
            {
                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    double normalStress_DU = ( -1 ) * SoilShapeFunctionValues[sec]
                                             * vPile[j] * penalty;

                    Vector tangentialStresses_DU( 2 );
            if (normTangentialStresses_trial != 0)
            {
              tangentialStresses_DU( 0 ) = friction_coeff
                                                 * normalStress_DU * tangentialStresses_trial[0] / normTangentialStresses_trial;

              tangentialStresses_DU( 1 ) = friction_coeff
                                                 * normalStress_DU * tangentialStresses_trial[1] / normTangentialStresses_trial;
            } else
            {
              tangentialStresses_DU( 0 ) = 0;
              tangentialStresses_DU( 1 ) = 0;
            }
                    rLeftHandSideMatrix( prim*dimension+ i, pileNN*dimension+ sec*dimension+ j ) +=
                        ( tangentialStresses_DU[0] * XiPrim[0] +
                          tangentialStresses_DU[1] * XiPrim[1] )
                        * SoilIntegrationWeight * influence_area;
                }
            }
        }
    }

    //END OF CONTRIBUTION: MASTER-SLAVE
    //BEGIN OF CONTRIBUTION: SLAVE -MASTER
    for ( unsigned int prim = 0; prim < soilNN; prim++ )
    {
        for ( unsigned int i = 0; i < dimension; i++ )
        {
            Vector XiPrim = ZeroVector( 2 );
        if ((norm_T (0)) != 0)
        {

          XiPrim[0] = SoilShapeFunctionValues[prim] * T( 0, i ) / norm_T( 0 );
        }
        else
        {
         XiPrim[0]=0;
        }

        if ((norm_T (1)) != 0)
        {

            XiPrim[1] = SoilShapeFunctionValues[prim] * T( 1, i ) / norm_T( 1 );
        }
        else
        {
          XiPrim[1]= 0;
        }
            for ( unsigned int sec = 0; sec < pileNN; sec++ )
            {
                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    double normalStress_DU = PileShapeFunctionValues[sec]
                                             * vPile[j] * penalty;

                    Vector tangentialStresses_DU( 2 );
            if (normTangentialStresses_trial != 0)
            {
              tangentialStresses_DU( 0 ) = friction_coeff
                                                 * normalStress_DU * tangentialStresses_trial[0] / normTangentialStresses_trial;

              tangentialStresses_DU( 1 ) = friction_coeff
                                                 * normalStress_DU * tangentialStresses_trial[1] / normTangentialStresses_trial;
            } else
            {
              tangentialStresses_DU( 0 ) = 0;
              tangentialStresses_DU( 1 ) = 0;
            }
                    rLeftHandSideMatrix( prim*dimension+ i + pileNN*dimension, sec*dimension+ j )     +=
                        ( tangentialStresses_DU[0] * XiPrim[0] +
                          tangentialStresses_DU[1] * XiPrim[1] )
                        * SoilIntegrationWeight * influence_area;
                }
            }
        }
    }

    //BEGIN OF CONTRIBUTION: SLAVE -MASTER
    //BEGIN OF CONTRIBUTION: SLAVE -SLAVE
    for ( unsigned int prim = 0; prim < soilNN; prim++ )
    {
        for ( unsigned int i = 0; i < dimension; i++ )
        {
            Vector XiPrim = ZeroVector( 2 );

        if ((norm_T (0)) != 0)
        {

          XiPrim[0] = SoilShapeFunctionValues[prim] * T( 0, i ) / norm_T( 0 );
        }
        else
        {
         XiPrim[0]=0;
        }

        if ((norm_T (1)) != 0)
        {

            XiPrim[1] = SoilShapeFunctionValues[prim] * T( 1, i ) / norm_T( 1 );
        }
        else
        {
          XiPrim[1]= 0;
        }

            for ( unsigned int sec = 0; sec < soilNN; sec++ )
            {
                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    double normalStress_DU = ( -1 ) * SoilShapeFunctionValues[sec]
                                             * vPile[j] * penalty;

                    Vector tangentialStresses_DU( 2 );
            if (normTangentialStresses_trial != 0)
            {
                    tangentialStresses_DU( 0 ) = friction_coeff
                                                 * normalStress_DU * tangentialStresses_trial[0] / normTangentialStresses_trial;

                    tangentialStresses_DU( 1 ) = friction_coeff
                                                 * normalStress_DU * tangentialStresses_trial[1] / normTangentialStresses_trial;
            } else
            {
              tangentialStresses_DU( 0 ) = 0;
              tangentialStresses_DU( 1 ) = 0;
            }
                    rLeftHandSideMatrix( pileNN*dimension+ prim*dimension+ i, pileNN*dimension+ sec*dimension+ j )     +=
                        ( tangentialStresses_DU[0] * XiPrim[0] +
                          tangentialStresses_DU[1] * XiPrim[1] )
                        * SoilIntegrationWeight * influence_area;
                }
            }
        }
    }
   }
   else
   return;

    //BEGIN OF CONTRIBUTION: SLAVE -SLAVE
//             //*****************************************************************
//             //END OF CONTRIBUTION DUE TO LINEARIZATION OF FRICTIONAL STRESS
//              // ONLY IF STICK IS FALSE
//             //*****************************************************************
    KRATOS_CATCH( "" )
} // CalculateAll

/**
 * This function calculates the system contributions to the global damp matrix due to the contact problem
 * with regard to the current master and slave partners.
 * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node
 */
/**
 * This function calculates the system contributions to the global damp matrix due to the contact problem
 * with regard to the current Pile and Soil partners.
 * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node
 */
void Pile_Kinematic_Linear::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY
    double penalty =1.0e+10;
    double penalty_T =1.0e+10;
    double friction_coeff = 0.5;

    if ( !( friction_coeff > 0.0 ) )
    {
        return;
    }

    //**********************************************
    //setting up the dimensions of the contributions
    //**********************************************
    unsigned int pileNN = mpPileElement->GetGeometry().size();

    unsigned int soilNN = mpSoilElement->GetGeometry().size();

    unsigned int dimension = 3;

    //resizing as needed the LHS
    int MatSize = ( pileNN  + soilNN ) * dimension;

    rDampMatrix.resize( MatSize, MatSize, false );

    noalias( rDampMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS

    //calculating shape function values for current Soil element
    Vector SoilShapeFunctionValues( mpSoilElement->GetGeometry().size() );

    noalias( SoilShapeFunctionValues ) = ZeroVector( mpSoilElement->GetGeometry().size() );

    for ( IndexType PointNumber = 0;
            PointNumber < mpSoilElement->GetGeometry().size(); PointNumber++ )
    {
        SoilShapeFunctionValues[PointNumber]
        = mpSoilElement->GetGeometry().ShapeFunctionValue( PointNumber,
                mSoilLocalPoint );
    }

    //calculating shape function values for current Pile element
    Vector PileShapeFunctionValues( mpPileElement->GetGeometry().size() );

    noalias( PileShapeFunctionValues ) = ZeroVector( mpPileElement->GetGeometry().size() );

    for ( IndexType PointNumber = 0;
            PointNumber < mpPileElement->GetGeometry().size(); PointNumber++ )
    {
        PileShapeFunctionValues[PointNumber]
        = mpPileElement->GetGeometry().ShapeFunctionValue( PointNumber,
                mPileLocalPoint );
    }

    //getting first order derivatives of Pile element shape functions
    Matrix PileDN = ZeroMatrix( pileNN, dimension - 2 );

    noalias( PileDN ) = mpPileElement->GetGeometry().ShapeFunctionsLocalGradients( PileDN,
                        mPileLocalPoint );

    Matrix SoilDN = ZeroMatrix( soilNN, dimension - 1 );

    noalias( SoilDN ) = mpSoilElement->GetGeometry().ShapeFunctionsLocalGradients( SoilDN,
                        mSoilLocalPoint );

    //********************************
    //Calculating contact area
    //********************************
    Matrix& inertia = mpPileElement->GetProperties()[INERTIA];
    double mArea = mpPileElement->GetProperties()[CROSS_AREA];
    double mInertia_x = inertia(0,0);
    double mInertia_y = inertia(1,1);
    double mInertia_Polar = inertia(2,2);

    double circumference_pile = sqrt ((12/mArea)*(mInertia_x*mInertia_y)+2*mArea) ;

    double DetJPile = mpPileElement->GetGeometry().DeterminantOfJacobian( mPileIntegrationPointIndex );

    double influence_area = mpPileElement->GetGeometry().IntegrationPoints()[mPileIntegrationPointIndex].Weight()
                            * DetJPile * circumference_pile;

    Vector vPile( 3 );

    noalias( vPile ) = ZeroVector( 3 );

    //getting tangential vectors
    Matrix TPileN( 2, 3 );

    noalias( TPileN ) = ZeroMatrix( 2, 3 );


    //calculating tangential vectors


    for ( unsigned int n = 0; n < mpPileElement->GetGeometry().size(); n++ )
    {


        TPileN( 0, 0 ) += ( mpPileElement->GetGeometry()[n].X())* PileDN( n, 0 );
        TPileN( 0, 1 ) += ( mpPileElement->GetGeometry()[n].Y())* PileDN( n, 0 );
        TPileN( 0, 2 ) += ( mpPileElement->GetGeometry()[n].Z())* PileDN( n, 0 );
        TPileN( 1, 0 ) +=0;
        TPileN( 1, 1 ) +=0;
        TPileN( 1, 2 ) +=0;


    }
    //calculating curent possition of the pile target point
    Vector vCurentPile( 3 );

    noalias( vCurentPile ) = ZeroVector( 3 );

    for ( unsigned int n = 0; n < mpPileElement->GetGeometry().size(); n++ )
    {

        {
            vCurentPile[0] += ( mpPileElement->GetGeometry()[n].X())*  PileShapeFunctionValues[n];
            vCurentPile[1] += ( mpPileElement->GetGeometry()[n].Y() )*  PileShapeFunctionValues[n];
            vCurentPile[2] += ( mpPileElement->GetGeometry()[n].Z())*  PileShapeFunctionValues[n];

        }
    }

    //calculating original (undeformd) possition of the pile target point

    Vector vOrigintPile( 3 );
    noalias( vOrigintPile ) = ZeroVector( 3 );

    for ( unsigned int n = 0; n < mpPileElement->GetGeometry().size(); n++ )
    {

        {
            vOrigintPile( 0 ) +=  (mpPileElement->GetGeometry()[n].X() - mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_X, 0 ))*  PileShapeFunctionValues[n];
            vOrigintPile( 1 ) +=   (mpPileElement->GetGeometry()[n].Y() - mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_Y, 0 ))*  PileShapeFunctionValues[n];
            vOrigintPile( 2 ) += ( mpPileElement->GetGeometry()[n].Z() - mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_Z, 0 ))*  PileShapeFunctionValues[n];
        }
    }

    //calculating projection of he curent possition of the pile target point

    Vector vOriginPileProj( 3 );
    {
        double T0 = ( vOrigintPile( 0 ) * TPileN( 0, 0 ) +  vOrigintPile( 1 ) * TPileN( 0, 1 ) +  vOrigintPile( 2 ) * TPileN( 0, 2 ) - (TPileN( 0, 0 )* vCurentPile[0] + TPileN( 0, 1 )* vCurentPile[1] + TPileN( 0, 2 )* vCurentPile[2])   ) / ( TPileN( 0, 0 ) * TPileN( 0, 0 ) + TPileN( 0, 1 ) * TPileN( 0, 1 ) + TPileN( 0, 2 ) * TPileN( 0, 2 ) ) ;

        noalias( vOriginPileProj ) = ZeroVector( 3 );

        vOriginPileProj( 0 ) += ( vOrigintPile( 0 ) - TPileN( 0, 0 ) * T0 );
        vOriginPileProj( 1 ) += ( vOrigintPile( 1 ) - TPileN( 0, 1 ) * T0 );
        vOriginPileProj( 2 ) += ( vOrigintPile( 2 ) - TPileN( 0, 2 ) * T0 );

    }


    noalias( mPileGlobalPoint ) = ZeroVector( 3 );
    noalias( mSoilGlobalPoint ) = ZeroVector( 3 );

    noalias( mPileGlobalPoint ) = GetGlobalCoordinates( mpPileElement, mPileGlobalPoint, mPileLocalPoint );

    noalias( mSoilGlobalPoint ) = GetGlobalCoordinates( mpSoilElement, mSoilGlobalPoint, mSoilLocalPoint );
    Vector relDispS( 3 );

    noalias( relDispS ) = ( mPileGlobalPoint - mSoilGlobalPoint );
    double NormrelDispS = MathUtils<double>::Norm3( relDispS );
    //KRATOS_WATCH( relDispS );

    //calculating normal vector on pile

    vPile = vCurentPile - vOriginPileProj;

    SD_MathUtils<double>::Normalize( vPile );
    double norm_vPile;
    norm_vPile = sqrt( vPile( 0 ) * vPile( 0 ) + vPile( 1 ) * vPile( 1 ) + vPile(2 ) * vPile( 2));

    if (norm_vPile!= 0)
    {
        TPileN( 1, 0 ) = TPileN( 0, 1 ) * vPile[2] - TPileN( 0, 2 ) * vPile[1];
        TPileN( 1, 1 ) = TPileN( 0, 2 ) * vPile[0] - TPileN( 0, 0 ) * vPile[2];
        TPileN( 1, 2 ) = TPileN( 0, 0 ) * vPile[1] - TPileN( 0, 1 ) * vPile[0];
    }
    else
    {
        TPileN( 1, 0 ) ==0;
        TPileN( 1, 1 ) ==0;
        TPileN( 1, 2 ) ==0;
    }

    noalias( mPileGlobalPoint ) = ZeroVector( 3 );
    noalias( mSoilGlobalPoint ) = ZeroVector( 3 );

    noalias( mPileGlobalPoint ) = GetGlobalCoordinates( mpPileElement, mPileGlobalPoint, mPileLocalPoint );

    noalias( mSoilGlobalPoint ) = GetGlobalCoordinates( mpSoilElement, mSoilGlobalPoint, mSoilLocalPoint );
    Vector relDisp( 3 );

    noalias( relDisp ) = ( mPileGlobalPoint - mSoilGlobalPoint );
    //KRATOS_WATCH( relDisp );

    double Gap;

    if (norm_vPile !=0)
    {
        Gap = inner_prod(
                  vPile, ( mPileGlobalPoint - mSoilGlobalPoint ) );
    }
    else
    {
        Gap ==0;
    }


    //calculating normal contact stress
    double normalStress = penalty * Gap;


    normalStress = sqrt( normalStress * normalStress );



    Matrix T( 2, 3 );
//       noalias(T) = TangentialVectors( mpPileElement, mPileLocalPoint );
    noalias( T ) = ZeroMatrix( 2, 3 );
    //calculating tangential vectors



    for ( unsigned int n = 0; n < mpPileElement->GetGeometry().size(); n++ )
    {


        T( 0, 0 ) += ( mpPileElement->GetGeometry()[n].X() )* PileDN( n, 0 );
        T( 0, 1 ) += ( mpPileElement->GetGeometry()[n].Y() )* PileDN( n, 0 );
        T( 0, 2 ) += ( mpPileElement->GetGeometry()[n].Z())* PileDN( n, 0 );

        T( 1, 0 ) +=  0;
        T( 1, 1 ) += 0;
        T( 1, 2 ) +=  0;
    }

    Vector vCurentPileB( 3 );

    noalias( vCurentPileB ) = ZeroVector( 3 );

    for ( unsigned int n = 0; n < pileNN; n++ )
    {
        vCurentPileB( 0 ) += ( mpPileElement->GetGeometry()[n].X())*  PileShapeFunctionValues[n];
        vCurentPileB( 1 ) += ( mpPileElement->GetGeometry()[n].Y())*  PileShapeFunctionValues[n];
        vCurentPileB( 2 ) += ( mpPileElement->GetGeometry()[n].Z())*  PileShapeFunctionValues[n];
    }

    Vector vOrigintPileB( 3 );

    noalias( vOrigintPileB ) = ZeroVector( 3 );

    for ( unsigned int n = 0; n < pileNN; n++ )
    {
        vOrigintPileB( 0 ) += (( mpPileElement->GetGeometry()[n].X() - mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_X, 0 ) )
                               *  PileShapeFunctionValues[n]);
        vOrigintPileB( 1 ) += (( mpPileElement->GetGeometry()[n].Y() - mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_Y, 0 ))
                               *  PileShapeFunctionValues[n]);
        vOrigintPileB( 2 ) += (( mpPileElement->GetGeometry()[n].Z() -mpPileElement->GetGeometry()[n].GetSolutionStepValue( DISPLACEMENT_Z, 0 ) )
                               *  PileShapeFunctionValues[n]);
    }

    Vector vOriginPileProjB( 3 );
    {
        double T0B = ( vOrigintPileB( 0 ) * T( 0, 0 ) +  vOrigintPileB( 1 ) * T( 0, 1 ) +  vOrigintPileB( 2 ) * T( 0, 2 ) - (T( 0, 0 )* vCurentPileB[0] + T( 0, 1 )* vCurentPileB[1] + T( 0, 2 )* vCurentPileB[2])   ) / ( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) ) ;

        double T0BX = ( vOrigintPileB( 0 ) * T( 0, 0 ) +  vOrigintPileB( 1 ) * T( 0, 1 ) +  vOrigintPileB( 2 ) * T( 0, 2 )   ) / ( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) ) ;


        noalias( vOriginPileProjB ) = ZeroVector( 3 );

        vOriginPileProjB[0] = ( vOrigintPileB[0] - T( 0, 0 ) * T0B );
        vOriginPileProjB[1] = ( vOrigintPileB[1] - T( 0, 1 ) * T0B  );
        vOriginPileProjB[2] = ( vOrigintPileB[2] - T( 0, 2 ) * T0B );

    }


    noalias( mPileGlobalPoint ) = ZeroVector( 3 );
    noalias( mSoilGlobalPoint ) = ZeroVector( 3 );

    noalias( mPileGlobalPoint ) = GetGlobalCoordinates( mpPileElement, mPileGlobalPoint, mPileLocalPoint );

    noalias( mSoilGlobalPoint ) = GetGlobalCoordinates( mpSoilElement, mSoilGlobalPoint, mSoilLocalPoint );


    Vector vPileNonNormalized = ZeroVector( 3 );
    {
        vPileNonNormalized[0] = vCurentPileB( 0 ) - vOriginPileProjB( 0 ) ;
        vPileNonNormalized[1] = vCurentPileB( 1 ) - vOriginPileProjB( 1 ) ;
        vPileNonNormalized[2] = vCurentPileB( 2 ) - vOriginPileProjB( 2 ) ;
    }

    //calculating tangent cetor as a produc of noraml and calculaed angent (nxt1)
    double NormvPile = MathUtils<double>::Norm3( vPileNonNormalized );
    if (NormvPile !=0)
    {
        T( 1, 0 ) = T( 0, 1 ) * vPileNonNormalized[2] - T( 0, 2 ) * vPileNonNormalized[1];
        T( 1, 1 ) = T( 0, 2 ) * vPileNonNormalized[0] - T( 0, 0 ) * vPileNonNormalized[2];
        T( 1, 2 ) = T( 0, 0 ) * vPileNonNormalized[1] - T( 0, 1 ) * vPileNonNormalized[0];
    }
    else
    {
        T( 1, 0 ) ==0;
        T( 1, 1 ) ==0;
        T( 1, 2 ) ==0;
    }
    Vector norm_T( 2 );

    norm_T( 0 ) = sqrt( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) );

    norm_T( 1 ) = sqrt( T( 1, 0 ) * T( 1, 0 ) + T( 1, 1 ) * T( 1, 1 ) + T( 1, 2 ) * T( 1, 2 ) );


    //calculating compatibiliy vector m (due to orogonaliy condition it is usnity vector )
    Matrix m( 2, 2 );
    noalias( m ) = ZeroMatrix( 2, 2 );//m[alpha][beta]=T[alpha]*T[beta]


      for ( int i = 0; i < 2; i++ )
      {
        for ( int j = 0; j < 2; j++ )
        {
            m( i, j ) ==1;
        }
      }


    Vector tangentialStresses( 2 );

    noalias( tangentialStresses ) = ZeroVector( 2 );
    Vector tangentialStresses_trial( 2 );
    noalias( tangentialStresses_trial ) = ZeroVector( 2 );
    Vector tangentialVelocity( 2 );
    noalias( tangentialVelocity ) = ZeroVector( 2 );
    double normTangentialStresses_trial;

    bool Stick = true;

    //TEST Calculation of relative velocity between Pile and Soil surface
// noalias(tangentialVelocity)= GetRelativTangentialVelocity(T);
//  Vector tangentialVelocity(2);
    if (NormrelDispS != 0)
    {
        noalias( tangentialVelocity ) = ZeroVector( 2 );
        {
            Vector soil_velo( 3 );
            Vector pile_velo( 3 );
            noalias( soil_velo ) = ZeroVector( 3 );
            noalias( pile_velo ) = ZeroVector( 3 );

            for ( IndexType i = 0 ; i < mpSoilElement->GetGeometry().size() ; i++ )
            {
                double shape_func = mpSoilElement->GetGeometry().ShapeFunctionValue( i, mSoilLocalPoint );
                soil_velo += (( mpSoilElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) ) *
                             shape_func;
            }

            for ( IndexType i = 0 ; i < mpPileElement->GetGeometry().size() ; i++ )
            {
                double shape_func = mpPileElement->GetGeometry().ShapeFunctionValue( i, mPileLocalPoint );
                pile_velo += (( mpPileElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) ) *
                             shape_func;
            }

            Vector norm_T( 2 );

            norm_T( 0 ) = sqrt( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) );

            norm_T( 1 ) = sqrt( T( 1, 0 ) * T( 1, 0 ) + T( 1, 1 ) * T( 1, 1 ) + T( 1, 2 ) * T( 1, 2 ) );

            tangentialVelocity( 0 ) = (( soil_velo( 0 ) - pile_velo( 0 ) ) * T( 0, 0 ) + ( soil_velo( 1 ) - pile_velo( 1 ) ) * T( 0, 1 )
                                       + ( soil_velo( 2 ) - pile_velo( 2 ) ) * T( 0, 2 ) ) / norm_T( 0 );

            tangentialVelocity( 1 ) = (( soil_velo( 0 ) - pile_velo( 0 ) ) * T( 1, 0 ) + ( soil_velo( 1 ) - pile_velo( 1 ) ) * T( 1, 1 )
                                       + ( soil_velo( 2 ) - pile_velo( 2 ) ) * T( 1, 2 ) ) / norm_T( 1 );
        }

        tangentialStresses_trial[0] += tangentialVelocity( 0 ) * penalty_T;

        tangentialStresses_trial[1] += tangentialVelocity( 1 ) * penalty_T;

        normTangentialStresses_trial =  sqrt( tangentialStresses_trial[0] * m( 0, 0 ) * tangentialStresses_trial[0]
                                              + tangentialStresses_trial[0] * m( 0, 1 ) * tangentialStresses_trial[1]
                                              + tangentialStresses_trial[1] * m( 1, 0 ) * tangentialStresses_trial[0]
                                              + tangentialStresses_trial[1] * m( 1, 1 ) * tangentialStresses_trial[1] );
    }
    else
    {
        tangentialStresses_trial[0] ==0;

        tangentialStresses_trial[1]== 0;

        normTangentialStresses_trial == 0;
    }

    if (NormrelDispS == 0)
    {
        tangentialStresses[0] == 0;
        tangentialStresses[1] == 0;
    }
    else if ( normTangentialStresses_trial > friction_coeff*normalStress )//Slip
    {
        tangentialStresses[0] = friction_coeff
                                * normalStress * tangentialStresses_trial[0] / normTangentialStresses_trial;//
        tangentialStresses[1] = friction_coeff
                                * normalStress * tangentialStresses_trial[1] / normTangentialStresses_trial;//

        Stick = false;
    }
    else//Stick
    {
        tangentialStresses[0] = tangentialStresses_trial[0];
        tangentialStresses[1] = tangentialStresses_trial[1];
        Stick = true;
    }
    //KRATOS_WATCH (tangentialStresses_trial);
    //KRATOS_WATCH (tangentialStresses);

    //Calculating Soil element's current integration weight
    double SoilIntegrationWeight = mpSoilElement->GetGeometry().IntegrationPoints()[mPileIntegrationPointIndex].Weight();


    //BEGIN OF CONTRIBUTION: MASTER-MASTER
    // std::cout << "MASTER-MASTER" << std::endl;
    if (Gap !=0)
    {
        for ( unsigned int prim = 0; prim < pileNN; prim++ )
        {
            for ( unsigned int i = 0; i < dimension; i++ )
            {
                Vector Xi = ZeroVector( 2 );
        if ((norm_T (0)) != 0)
        {
          Xi[0] = -PileShapeFunctionValues[prim] * T( 0, i ) / norm_T( 0 );
        }
        else
        {
        Xi[0]=0;
        }

        if ((norm_T (1)) != 0)
        {
          Xi[1] = -PileShapeFunctionValues[prim] * T( 1, i ) / norm_T( 1 );
        }
        else
        {
          Xi[1]= 0;
        }

                for ( unsigned int sec = 0; sec < pileNN; sec++ )
                {
                    for ( unsigned int j = 0; j < dimension; j++ )
                    {
                        //Derivative RelativeTangentialVelocity/MasterVelocity
                        Vector tangentialStresses_DV( 2 );

                        Vector tangentialStresses_trial_DV( 2 );
            if ((norm_T (0)) != 0)
            {
                        tangentialStresses_trial_DV( 0 ) =
                            penalty_T * T( 0, j ) / norm_T( 0 ) * ( -1 ) * PileShapeFunctionValues[sec];
            }else
            {tangentialStresses_trial_DV( 0 ) = 0;
            }
            if ((norm_T (1)) != 0)
            {
                        tangentialStresses_trial_DV( 1 ) =
                            penalty_T * T( 1, j ) / norm_T( 1 ) * ( -1 ) * PileShapeFunctionValues[sec];
            }else
            {tangentialStresses_trial_DV( 1 ) = 0;
            }

                        if ( Stick )
                        {
                            tangentialStresses_DV( 0 ) =
                                tangentialStresses_trial_DV( 0 );

                            tangentialStresses_DV( 1 ) =
                                tangentialStresses_trial_DV( 1 );
                        }
                        else
                        {
              if (normTangentialStresses_trial !=0 )
              {
                            tangentialStresses_DV( 0 ) =
                                ( penalty_T * T( 0, j ) / norm_T( 0 ) * ( friction_coeff
                                        * normalStress ) / normTangentialStresses_trial
                                  + penalty_T * ( tangentialStresses_trial[0] * friction_coeff
                                                  * normalStress ) / ( 2 * normTangentialStresses_trial )
                                  * 2 * ( T( 0, j ) / norm_T( 0 ) * m( 0, 0 ) * tangentialStresses_trial[0] + T( 0, j ) / norm_T( 0 ) * m( 0, 1 ) * tangentialStresses_trial[1] + T( 1, j ) / norm_T( 1 ) * m( 1, 0 ) * tangentialStresses_trial[0] + T( 1, j ) / norm_T( 1 ) * m( 1, 1 ) * tangentialStresses_trial[1] ) ) * PileShapeFunctionValues[sec];
                            tangentialStresses_DV( 1 ) =
                                ( penalty_T * T( 1, j ) / norm_T( 1 ) * ( friction_coeff
                                        * normalStress ) / normTangentialStresses_trial
                                  + penalty_T * ( tangentialStresses_trial[1] * friction_coeff
                                                  * normalStress ) / ( 2 * normTangentialStresses_trial )
                                  * 2 * ( T( 0, j ) / norm_T( 0 ) * m( 0, 0 ) * tangentialStresses_trial[0] + T( 0, j ) / norm_T( 0 ) * m( 0, 1 ) * tangentialStresses_trial[1] + T( 1, j ) / norm_T( 1 ) * m( 1, 0 ) * tangentialStresses_trial[0] + T( 1, j ) / norm_T( 1 ) * m( 1, 1 ) * tangentialStresses_trial[1] ) ) * PileShapeFunctionValues[sec];
              }
              else
              {
                            tangentialStresses_DV( 0 ) = 0;
                             tangentialStresses_DV( 1 ) = 0;
              }
                        }

                        rDampMatrix( prim*dimension+ i, sec*dimension+ j )

                        += ( tangentialStresses_DV[0] * Xi[0] + tangentialStresses_DV[1] * Xi[1] )
                           * SoilIntegrationWeight * influence_area;
                    }
                }
            }
        }


        //END OF CONTRIBUTION: MASTER-MASTER


        //BEGIN OF CONTRIBUTION MASTER-SLAVE
        for ( unsigned int prim = 0; prim < pileNN; prim++ )
        {
            for ( unsigned int i = 0; i < dimension; i++ )
            {
              Vector Xi = ZeroVector( 2 );
        if ((norm_T (0)) != 0)
        {
          Xi[0] = -PileShapeFunctionValues[prim] * T( 0, i ) / norm_T( 0 );
        }
        else
        {
        Xi[0]=0;
        }

        if ((norm_T (1)) != 0)
        {
          Xi[1] = -PileShapeFunctionValues[prim] * T( 1, i ) / norm_T( 1 );
        }
        else
        {
          Xi[1]= 0;
        }
                for ( unsigned int sec = 0; sec < soilNN; sec++ )
                {

                    for ( unsigned int j = 0; j < dimension; j++ )
                    {
                        Vector tangentialStresses_DV( 2 );

                        Vector tangentialStresses_trial_DV( 2 );
            if ((norm_T (0)) != 0)
            {

                        tangentialStresses_trial_DV( 0 ) =
                            penalty_T * T( 0, j ) / norm_T( 0 ) * SoilShapeFunctionValues[sec];
            }
            else
            {tangentialStresses_trial_DV( 0 ) = 0;
            }
            if ((norm_T (1)) != 0)
            {

                        tangentialStresses_trial_DV( 1 ) =
                            penalty_T * T( 1, j ) / norm_T( 1 ) * SoilShapeFunctionValues[sec];
            }else
            {tangentialStresses_trial_DV( 1 ) = 0;
            }

                        if ( Stick )
                        {
                            tangentialStresses_DV( 0 ) =
                                tangentialStresses_trial_DV( 0 );

                            tangentialStresses_DV( 1 ) =
                                tangentialStresses_trial_DV( 1 );
                        }
                        else
                        {
              if (normTangentialStresses_trial !=0 )
              {                            tangentialStresses_DV( 0 ) =
                                -( penalty_T * T( 0, j ) / norm_T( 0 ) * ( friction_coeff
                                        * normalStress ) / normTangentialStresses_trial
                                   + penalty_T * ( tangentialStresses_trial[0] * friction_coeff
                                                   * normalStress ) / ( 2 * normTangentialStresses_trial )
                                   * 2 * ( T( 0, j ) / norm_T( 0 ) * m( 0, 0 ) * tangentialStresses_trial[0] + T( 0, j ) / norm_T( 0 ) * m( 0, 1 ) * tangentialStresses_trial[1] + T( 1, j ) / norm_T( 1 ) * m( 1, 0 ) * tangentialStresses_trial[0] + T( 1, j ) / norm_T( 1 ) * m( 1, 1 ) * tangentialStresses_trial[1] ) ) * PileShapeFunctionValues[sec];
                            tangentialStresses_DV( 1 ) =
                                -( penalty_T * T( 1, j ) / norm_T( 1 ) * ( friction_coeff
                                        * normalStress ) / normTangentialStresses_trial
                                   + penalty_T * ( tangentialStresses_trial[1] * friction_coeff
                                                   * normalStress ) / ( 2 * normTangentialStresses_trial )
                                   * 2 * ( T( 0, j ) / norm_T( 0 ) * m( 0, 0 ) * tangentialStresses_trial[0] + T( 0, j ) / norm_T( 0 ) * m( 0, 1 ) * tangentialStresses_trial[1] + T( 1, j ) / norm_T( 1 ) * m( 1, 0 ) * tangentialStresses_trial[0] + T( 1, j ) / norm_T( 1 ) * m( 1, 1 ) * tangentialStresses_trial[1] ) ) * PileShapeFunctionValues[sec];
              }
              else
              {
                            tangentialStresses_DV( 0 ) = 0;
                             tangentialStresses_DV( 1 ) = 0;
              }
             }

                        rDampMatrix( prim*dimension+ i, pileNN*dimension+ sec*dimension+ j )

                        += ( tangentialStresses_DV[0] * Xi[0] + tangentialStresses_DV[1] * Xi[1] )
                           * SoilIntegrationWeight * influence_area;
                    }
                }
            }
        }

        //END OF CONTRIBUTION: MASTER-SLAVE
        //BEGIN OF CONTRIBUTION: SLAVE -MASTER
        for ( unsigned int prim = 0; prim < soilNN; prim++ )
        {
            for ( unsigned int i = 0; i < dimension; i++ )
            {

                Vector Xi = ZeroVector( 2 );

        if ((norm_T (0)) != 0)
        {

          Xi[0] = SoilShapeFunctionValues[prim] * T( 0, i ) / norm_T( 0 );
        }
        else
        {
        Xi[0]=0;
        }

        if ((norm_T (1)) != 0)
        {

        Xi[1] = SoilShapeFunctionValues[prim] * T( 1, i ) / norm_T( 1 );
        }
        else
        {
          Xi[1]= 0;
        }

                for ( unsigned int sec = 0; sec < pileNN;  sec++ )
                {

                    for ( unsigned int j = 0; j < dimension; j++ )
                    {
                        Vector tangentialStresses_DV( 2 );

                        Vector tangentialStresses_trial_DV( 2 );
            if ((norm_T (0)) != 0)
            {
                        tangentialStresses_trial_DV( 0 ) =
                            penalty_T * T( 0, j ) / norm_T( 0 ) * ( -1 ) * PileShapeFunctionValues[sec];
            }else
            {tangentialStresses_trial_DV( 0 ) = 0;
            }
            if ((norm_T (1)) != 0)
            {

                        tangentialStresses_trial_DV( 1 ) =
                            penalty_T * T( 1, j ) / norm_T( 1 ) * ( -1 ) * PileShapeFunctionValues[sec];
            }else

            {
              tangentialStresses_trial_DV( 1 ) = 0;
            }
                        if ( Stick )
                        {
                            tangentialStresses_DV( 0 ) =
                                tangentialStresses_trial_DV( 0 );

                            tangentialStresses_DV( 1 ) =
                                tangentialStresses_trial_DV( 1 );
                        }
                        else
                        {
              if (normTangentialStresses_trial !=0 )
              {
                            tangentialStresses_DV( 0 ) =
                                ( penalty_T * T( 0, j ) / norm_T( 0 ) * ( friction_coeff
                                        * normalStress ) / normTangentialStresses_trial
                                  + penalty_T * ( tangentialStresses_trial[0] * friction_coeff
                                                  * normalStress ) / ( 2 * normTangentialStresses_trial )
                                  * 2 * ( T( 0, j ) / norm_T( 0 ) * m( 0, 0 ) * tangentialStresses_trial[0] + T( 0, j ) / norm_T( 0 ) * m( 0, 1 ) * tangentialStresses_trial[1] + T( 1, j ) / norm_T( 1 ) * m( 1, 0 ) * tangentialStresses_trial[0] + T( 1, j ) / norm_T( 1 ) * m( 1, 1 ) * tangentialStresses_trial[1] ) ) * PileShapeFunctionValues[sec];
                            tangentialStresses_DV( 1 ) =
                                ( penalty_T * T( 1, j ) / norm_T( 1 ) * ( friction_coeff
                                        * normalStress ) / normTangentialStresses_trial
                                  + penalty_T * ( tangentialStresses_trial[1] * friction_coeff
                                                  * normalStress ) / ( 2 * normTangentialStresses_trial )
                                  * 2 * ( T( 0, j ) / norm_T( 0 ) * m( 0, 0 ) * tangentialStresses_trial[0] + T( 0, j ) / norm_T( 0 ) * m( 0, 1 ) * tangentialStresses_trial[1] + T( 1, j ) / norm_T( 1 ) * m( 1, 0 ) * tangentialStresses_trial[0] + T( 1, j ) / norm_T( 1 ) * m( 1, 1 ) * tangentialStresses_trial[1] ) ) * PileShapeFunctionValues[sec];
              }
              else
              {
                            tangentialStresses_DV( 0 ) = 0;
                             tangentialStresses_DV( 1 ) = 0;
              }
             }

                        rDampMatrix( pileNN*dimension+ prim*dimension+ i, sec*dimension+ j )

                        += ( tangentialStresses_DV[0] * Xi[0] + tangentialStresses_DV[1] * Xi[1] )
                           * SoilIntegrationWeight * influence_area;
                    }
                }
            }
        }

        //END OF CONTRIBUTION: SLAVE-MASTER
        //BEGIN OF CONTRIBUTION: SLAVE -SLAVE
        for ( unsigned int prim = 0;  prim < soilNN; prim++ )
        {
            for ( unsigned int i = 0; i < dimension; i++ )
            {
                Vector Xi = ZeroVector( 2 );

        if ((norm_T (0)) != 0)
        {

          Xi[0] = SoilShapeFunctionValues[prim] * T( 0, i ) / norm_T( 0 );
        }
        else
        {
        Xi[0]=0;
        }

        if ((norm_T (1)) != 0)
        {

        Xi[1] = SoilShapeFunctionValues[prim] * T( 1, i ) / norm_T( 1 );
        }
        else
        {
          Xi[1]= 0;
        }

                for ( unsigned int sec = 0; sec < soilNN; sec++ )
                {
                    for ( unsigned int j = 0; j < dimension; j++ )
                    {
                        Vector tangentialStresses_DV( 2 );

                        Vector tangentialStresses_trial_DV( 2 );
            if ((norm_T (0)) != 0)
            {
                        tangentialStresses_trial_DV( 0 ) =
                            penalty_T * T( 0, j ) / norm_T( 0 ) * SoilShapeFunctionValues[sec];
            }else
            {tangentialStresses_trial_DV( 0 ) = 0;
            }
            if ((norm_T (1)) != 0)
            {
                        tangentialStresses_trial_DV( 1 ) =
                            penalty_T * T( 1, j ) / norm_T( 1 ) * SoilShapeFunctionValues[sec];
            }else
            {
              tangentialStresses_trial_DV( 1 ) = 0;
            }


                        if ( Stick )
                        {
                            tangentialStresses_DV( 0 ) =
                                tangentialStresses_trial_DV( 0 );

                            tangentialStresses_DV( 1 ) =
                                tangentialStresses_trial_DV( 1 );
                        }
                        else
                        {
              if (normTangentialStresses_trial !=0 )
              {
                            tangentialStresses_DV( 0 ) =
                                -( penalty_T * T( 0, j ) / norm_T( 0 ) * ( friction_coeff
                                        * normalStress ) / normTangentialStresses_trial
                                   + penalty_T * ( tangentialStresses_trial[0] * friction_coeff
                                                   * normalStress ) / ( 2 * normTangentialStresses_trial )
                                   * 2 * ( T( 0, j ) / norm_T( 0 ) * m( 0, 0 ) * tangentialStresses_trial[0] + T( 0, j ) / norm_T( 0 ) * m( 0, 1 ) * tangentialStresses_trial[1] + T( 1, j ) / norm_T( 1 ) * m( 1, 0 ) * tangentialStresses_trial[0] + T( 1, j ) / norm_T( 1 ) * m( 1, 1 ) * tangentialStresses_trial[1] ) ) * PileShapeFunctionValues[sec];
                            tangentialStresses_DV( 1 ) =
                                -( penalty_T * T( 1, j ) / norm_T( 1 ) * ( friction_coeff
                                        * normalStress ) / normTangentialStresses_trial
                                   + penalty_T * ( tangentialStresses_trial[1] * friction_coeff
                                                   * normalStress ) / ( 2 * normTangentialStresses_trial )
                                   * 2 * ( T( 0, j ) / norm_T( 0 ) * m( 0, 0 ) * tangentialStresses_trial[0] + T( 0, j ) / norm_T( 0 ) * m( 0, 1 ) * tangentialStresses_trial[1] + T( 1, j ) / norm_T( 1 ) * m( 1, 0 ) * tangentialStresses_trial[0] + T( 1, j ) / norm_T( 1 ) * m( 1, 1 ) * tangentialStresses_trial[1] ) ) * PileShapeFunctionValues[sec];
              }
              else
              {
                            tangentialStresses_DV( 0 ) = 0;
                             tangentialStresses_DV( 1 ) = 0;
              }
            }

                        rDampMatrix( pileNN*dimension+ prim*dimension+ i, pileNN*dimension+ sec*dimension+ j )

                        += ( tangentialStresses_DV[0] * Xi[0] + tangentialStresses_DV[1] * Xi[1] )
                           * SoilIntegrationWeight * influence_area;
                    }
                }
            }
        }
    }
    else
        return;

    KRATOS_CATCH( "" )
}

//***********************************************************************
//***********************************************************************
/**
 * Calculates the residual contribution due to contact stresses on
 * both the Pile and Soil conditions. The contributions are stored
 * into the residual vector.
 * @param residualvector the right hand side vector of the current
 * linking condition
 * @param NPile the shape function values of the Pile element in
 * the current contact point
 * @param NSoil the shape function values of the Soil element in
 * the current contact point
 * @param vPile the normal vector to the Pile surface in
 * current contact point
 * @param normalStress the value of contact stress in current contact
 * point
 * @param SoilIntegrationWeight the integration weight in current Soil
 * integration point
 * @param dASoil the differential area element of the current Soil
 * integration point
 */
void Pile_Kinematic_Linear::CalculateAndAdd_RHS( Vector& residualvector,
        const Vector& NPile,
        const Vector& NSoil,
        const Vector& vPile,
        const Matrix& T,
        const Vector& tangentialStresses,
        double Gap,
        double normalStress,
        double SoilIntegrationWeight,
        double dASoil,
        double friction_coeff
                                               )
{
    //**********************************************
    //BEGIN OF Pile: CalculateAndAdd_PressureForce
    //**********************************************
//    //KRATOS_WATCH ("line://BEGIN OF Pile: CalculateAndAdd_PressureForce ")

    unsigned int pileNN = mpPileElement->GetGeometry().size();
    unsigned int dimension = 3;

    for ( unsigned int i = 0; i < pileNN; i++ )
    {
        int index = dimension * i;
        double coeff = -normalStress * NPile[i] * SoilIntegrationWeight * dASoil;
        residualvector[index]   = coeff * vPile[0];
        residualvector[index+1] = coeff * vPile[1];
        residualvector[index+2] = coeff * vPile[2];
//        //KRATOS_WATCH (normalStress);
//        //KRATOS_WATCH (residualvector);
    }


    //********************************************
    //END OF Pile: CalculateAndAdd_PressureForce
    //********************************************

    //*********************************************
    //BEGIN OF Soil: CalculateAndAdd_PressureForce
    //*********************************************

    unsigned int soilNN = mpSoilElement->GetGeometry().size();

    for ( unsigned int i = 0; i < soilNN; i++ )
    {
        int index =  pileNN * dimension + dimension * i;
        double coeff = normalStress * NSoil[i] * SoilIntegrationWeight * dASoil;
        residualvector[index]   = coeff * vPile[0];
        residualvector[index+1] = coeff * vPile[1];
        residualvector[index+2] = coeff * vPile[2];
    }

    //*******************************************
    //END OF Soil: CalculateAndAdd_PressureForce
    //*******************************************
//         //KRATOS_WATCH( residualvector );


    //*******************************************
    //CalculateAndAdd RHS-Vector due to frictional contact
    //*******************************************


    //**********************************************
    //BEGIN OF Pile: CalculateAndAdd_PressureForce
    //**********************************************

    if ( friction_coeff > 0.0 )
    {

        Vector norm_T( 2 );

        norm_T( 0 ) = sqrt( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) );
        norm_T( 1 ) = sqrt( T( 1, 0 ) * T( 1, 0 ) + T( 1, 1 ) * T( 1, 1 ) + T( 1, 2 ) * T( 1, 2 ) );

        for ( unsigned int i = 0; i < pileNN; i++ )
        {
            for ( unsigned int j = 0 ; j < dimension ; j++ )
            {
                Vector Xi = ZeroVector( 2 );

                Xi[0] = -NPile[i] * T( 0, j ) / norm_T( 0 );

                Xi[1] = -NPile[i] * T( 1, j ) / norm_T( 1 );

                residualvector[ i*dimension+j] -= ( tangentialStresses[0] * Xi[0] +
                                                      tangentialStresses[1] * Xi[1] ) * SoilIntegrationWeight * dASoil;
            }
        }

        //**********************************************
        //END OF Pile: CalculateAndAdd_PressureForce
        //**********************************************

        //**********************************************
        //BEGIN OF Soil: CalculateAndAdd_PressureForce
        //**********************************************


        for ( unsigned int i = 0; i < soilNN; i++ )
        {
            for ( unsigned int j = 0;j < dimension;j++ )
            {
                Vector Xi = ZeroVector( 2 );

                Xi[0] = NSoil[i] * T( 0, j ) / norm_T( 0 );

                Xi[1] = NSoil[i] * T( 1, j ) / norm_T( 1 );

                residualvector[pileNN*dimension +i*dimension+j] -= ( tangentialStresses[0] * Xi[0] +
                        tangentialStresses[1] * Xi[1] ) * SoilIntegrationWeight * dASoil;
            }
        }


        //**********************************************
        //END OF Soil: CalculateAndAdd_PressureForce
        //**********************************************


    }

}

/**
 * This function calculates updates the local and global coordinates
 * of the Pile contact partner in order to follow the movement of
 * the Soil surface along the Pile surface
 */
void Pile_Kinematic_Linear::UpdatePileLocalPoint( )
{
    double Xi1 = mPileLocalPoint[0];
    double Xi2 = mPileLocalPoint[1];
    double deltaXi1 = 0.0;
    double deltaXi2 = 0.0;

    for ( int k = 0; k < 1000; k++ )
    {
        // //KRATOS_WATCH( mPileGlobalPoint );
        //setting up tangential vectors
        Vector t1 = ZeroVector( 3 );//first tangential vector
        Vector t2 = ZeroVector( 3 );//second tangential vector
        //derivatives of tangential vectors
        Vector dt11 = ZeroVector( 3 );
        Vector dt12 = ZeroVector( 3 );
        Vector dt21 = ZeroVector( 3 );
        Vector dt22 = ZeroVector( 3 );

        //retrieving first order derivatives in current solution point
        Matrix DN = ZeroMatrix( mpPileElement->GetGeometry().size(), 2 );
        mpPileElement->GetGeometry().ShapeFunctionsLocalGradients( DN, mPileLocalPoint );
//        //KRATOS_WATCH( DN );
        //retrieving second order derivatives in current solution point
        GeometryType::ShapeFunctionsSecondDerivativesType D2N;

        mpPileElement->GetGeometry().ShapeFunctionsSecondDerivatives( D2N, mPileLocalPoint );


        for ( unsigned  int n = 0; n < mpPileElement->GetGeometry().size(); n++ )
        {
            //contribution to tangential vectors
            t1[0] += ( mpPileElement->GetGeometry().GetPoint( n ).X0()
                       + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                     * DN( n, 0 );
            t1[1] += ( mpPileElement->GetGeometry().GetPoint( n ).Y0()
                       + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                     * DN( n, 0 );
            t1[2] += ( mpPileElement->GetGeometry().GetPoint( n ).Z0()
                       + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                     * DN( n, 0 );
            t2[0] += ( mpPileElement->GetGeometry().GetPoint( n ).X0()
                       + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                     * DN( n, 1 );
            t2[1] += ( mpPileElement->GetGeometry().GetPoint( n ).Y0()
                       + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                     * DN( n, 1 );
            t2[2] += ( mpPileElement->GetGeometry().GetPoint( n ).Z0()
                       + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                     * DN( n, 1 );
            //contribution to derivatives of tangential vectors
            dt11[0] += ( mpPileElement->GetGeometry().GetPoint( n ).X0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                       * D2N[n]( 0, 0 );
            dt11[1] += ( mpPileElement->GetGeometry().GetPoint( n ).Y0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                       * D2N[n]( 0, 0 );
            dt11[2] += ( mpPileElement->GetGeometry().GetPoint( n ).Z0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                       * D2N[n]( 0, 0 );
            dt12[0] += ( mpPileElement->GetGeometry().GetPoint( n ).X0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                       * D2N[n]( 0, 1 );
            dt12[1] += ( mpPileElement->GetGeometry().GetPoint( n ).Y0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                       * D2N[n]( 0, 1 );
            dt12[2] += ( mpPileElement->GetGeometry().GetPoint( n ).Z0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                       * D2N[n]( 0, 1 );
            dt21[0] += ( mpPileElement->GetGeometry().GetPoint( n ).X0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                       * D2N[n]( 1, 0 );
            dt21[1] += ( mpPileElement->GetGeometry().GetPoint( n ).Y0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                       * D2N[n]( 1, 0 );
            dt21[2] += ( mpPileElement->GetGeometry().GetPoint( n ).Z0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                       * D2N[n]( 1, 0 );
            dt22[0] += ( mpPileElement->GetGeometry().GetPoint( n ).X0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                       * D2N[n]( 1, 1 );
            dt22[1] += ( mpPileElement->GetGeometry().GetPoint( n ).Y0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                       * D2N[n]( 1, 1 );
            dt22[2] += ( mpPileElement->GetGeometry().GetPoint( n ).Z0()
                         + mpPileElement->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                       * D2N[n]( 1, 1 );
        }

        //defining auxiliary terms
        double A1 = (( mSoilGlobalPoint[0] - mPileGlobalPoint[0] ) * t1[0] )
                    + (( mSoilGlobalPoint[1] - mPileGlobalPoint[1] ) * t1[1] )
                    + (( mSoilGlobalPoint[2] - mPileGlobalPoint[2] ) * t1[2] );

        double A2 = (( mSoilGlobalPoint[0] - mPileGlobalPoint[0] ) * t2[0] )
                    + (( mSoilGlobalPoint[1] - mPileGlobalPoint[1] ) * t2[1] )
                    + (( mSoilGlobalPoint[2] - mPileGlobalPoint[2] ) * t2[2] );

        double B11 = ( -t1[0] * t1[0] - t1[1] * t1[1] - t1[2] * t1[2] )
                     + (( mSoilGlobalPoint[0] - mPileGlobalPoint[0] ) * dt11[0] )
                     + (( mSoilGlobalPoint[1] - mPileGlobalPoint[1] ) * dt11[1] )
                     + (( mSoilGlobalPoint[2] - mPileGlobalPoint[2] ) * dt11[2] );

        double B12 = ( -t2[0] * t1[0] - t2[1] * t1[1] - t2[2] * t1[2] )
                     + (( mSoilGlobalPoint[0] - mPileGlobalPoint[0] ) * dt12[0] )
                     + (( mSoilGlobalPoint[1] - mPileGlobalPoint[1] ) * dt12[1] )
                     + (( mSoilGlobalPoint[2] - mPileGlobalPoint[2] ) * dt12[2] );

        double B21 = ( -t1[0] * t2[0] - t1[1] * t2[1] - t1[2] * t2[2] )
                     + (( mSoilGlobalPoint[0] - mPileGlobalPoint[0] ) * dt21[0] )
                     + (( mSoilGlobalPoint[1] - mPileGlobalPoint[1] ) * dt21[1] )
                     + (( mSoilGlobalPoint[2] - mPileGlobalPoint[2] ) * dt21[2] );

        double B22 = ( -t2[0] * t2[0] - t2[1] * t2[1] - t2[2] * t2[2] )
                     + (( mSoilGlobalPoint[0] - mPileGlobalPoint[0] ) * dt22[0] )
                     + (( mSoilGlobalPoint[1] - mPileGlobalPoint[1] ) * dt22[1] )
                     + (( mSoilGlobalPoint[2] - mPileGlobalPoint[2] ) * dt22[2] );

        //calculating update for Xi
        deltaXi1 = -A1 * B22 / ( B11 * B22 - B12 * B21 ) + A2 * B12 / ( B11 * B22 - B12 * B21 );

        deltaXi2 =  A2 * B21 / ( B11 * B22 - B12 * B21 ) - A2 * B11 / ( B11 * B22 - B12 * B21 );

        //updating Xi
        Xi1 += deltaXi1;

        Xi2 += deltaXi2;

        //updating LocalPoint
        mPileLocalPoint[0] = Xi1;

        mPileLocalPoint[1] = Xi2;

        //updating rResult
        mPileGlobalPoint = ZeroVector( 3 );

        mPileGlobalPoint = GetGlobalCoordinates( mpPileElement, mPileGlobalPoint, mPileLocalPoint );

        if ( fabs( deltaXi1 ) < 1e-7 && fabs( deltaXi2 ) < 1e-7 )
        {
            return;
        }
    }

    std::cout << "******** ATTENTION: NO MAPPING TO Pile Line FOUND ************" << std::endl;

    return;
}

//************************************************************************************
//************************************************************************************
/**
 * Setting up the EquationIdVector for the current partners.
 * All conditions are assumed to be defined in 3D space with 3 DOFs per node.
 * All Equation IDs are given Pile first, Soil second
 */
void Pile_Kinematic_Linear::EquationIdVector( EquationIdVectorType& rResult,
        ProcessInfo& CurrentProcessInfo )
{
    //determining size of DOF list
    //dimension of space
    unsigned int dimension = 3;
    unsigned int pileNN = mpPileElement->GetGeometry().size();
    unsigned int soilNN = mpSoilElement->GetGeometry().size();
    unsigned int ndofs = dimension * ( soilNN + pileNN );
    unsigned int index;
    rResult.resize( ndofs, false );


    for ( unsigned int i = 0; i < pileNN; i++ )
    {
        index = i * dimension;
        rResult[index]   = ( mpPileElement->GetGeometry()[i].GetDof( DISPLACEMENT_X ) ).EquationId();
        rResult[index+1] = ( mpPileElement->GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId() );
        rResult[index+2] = ( mpPileElement->GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId() );

    }

    for ( unsigned  int i = 0; i < soilNN; i++ )
    {
        index =  pileNN * dimension + i * dimension;
        rResult[index]   = ( mpSoilElement->GetGeometry()[i].GetDof( DISPLACEMENT_X ) ).EquationId();
        rResult[index+1] = ( mpSoilElement->GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId() );
        rResult[index+2] = ( mpSoilElement->GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId() );

    }
}



//************************************************************************************
//************************************************************************************
/**
 * Setting up the DOF list for the current partners.
 * All conditions are assumed to be defined in 3D space with 3 DOFs per Node.
 * All DOF are given Pile first, Soil second
 */
void Pile_Kinematic_Linear::GetDofList( DofsVectorType& ConditionalDofList,
                                        ProcessInfo& CurrentProcessInfo )
{

    //determining size of DOF list
    //dimension of space
    unsigned int dimension = 3;
    unsigned int pileNN = mpPileElement->GetGeometry().size();
    unsigned int soilNN = mpSoilElement->GetGeometry().size();
    unsigned int ndofs = dimension * ( soilNN + pileNN );

    ConditionalDofList.resize( ndofs );
    unsigned int index;
    //setting up Pile DOFs


    for ( unsigned int i = 0; i < pileNN; i++ )
    {
        index = i * dimension;
        ConditionalDofList[index]   = ( mpPileElement->GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ConditionalDofList[index+1] = ( mpPileElement->GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ConditionalDofList[index+2] = ( mpPileElement->GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

    }

    //setting up Soil DOFs
    for ( unsigned int i = 0; i < soilNN; i++ )
    {
        index =  pileNN * dimension + i * dimension;
        ConditionalDofList[index]   = ( mpSoilElement->GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ConditionalDofList[index+1] = ( mpSoilElement->GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ConditionalDofList[index+2] = ( mpSoilElement->GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}


/**
   * returns the relative tangential velocity between the quadrature point on the Soil surface and its
   * closest point projection
   * @param T Matrix of the Tanegntial Vectors on the Pile Surface in the current configuration
   * @return tangential velocity
   */
//new functions includes

Vector Pile_Kinematic_Linear::GetRelativTangentialVelocity( Matrix& T )
{
    Vector result( 2 );

    Vector pile_velo( 3 );
    Vector soil_velo( 3 );

    noalias( pile_velo ) = ZeroVector( 3 );
    noalias( soil_velo ) = ZeroVector( 3 );

    for ( IndexType i = 0 ; i < mpPileElement->GetGeometry().size() ; i++ )
    {
        double shape_func = mpPileElement->GetGeometry().ShapeFunctionValue( i, mPileLocalPoint );
        pile_velo += (( mpPileElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT ) ) *
                     shape_func;
    }

    for ( IndexType i = 0 ; i < mpSoilElement->GetGeometry().size() ; i++ )
    {
        double shape_func = mpSoilElement->GetGeometry().ShapeFunctionValue( i, mSoilLocalPoint );
        soil_velo += (( mpSoilElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT ) ) *
                     shape_func;
    }

    Vector norm_T( 2 );

    norm_T( 0 ) = sqrt( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) );

    norm_T( 1 ) = sqrt( T( 1, 0 ) * T( 1, 0 ) + T( 1, 1 ) * T( 1, 1 ) + T( 1, 2 ) * T( 1, 2 ) );

    result( 0 ) = (( pile_velo( 0 ) - soil_velo( 0 ) ) * T( 0, 0 ) + ( pile_velo( 1 ) - soil_velo( 1 ) ) * T( 0, 1 )
                   + ( pile_velo( 2 ) - soil_velo( 2 ) ) * T( 0, 2 ) ) / norm_T( 0 );

    result( 1 ) = (( pile_velo( 0 ) - soil_velo( 0 ) ) * T( 1, 0 ) + ( pile_velo( 1 ) - soil_velo( 1 ) ) * T( 1, 1 )
                   + ( pile_velo( 2 ) - soil_velo( 2 ) ) * T( 1, 2 ) ) / norm_T( 1 );

    return result;
}

/**
 * returns the relative velocity between the quadrature point on the slave surface and its
 * closest point projection
 * @return relative velocity
 */
Vector Pile_Kinematic_Linear::GetRelativVelocity()
{
    Vector result( 3 );

    Vector pile_velo( 3 );
    Vector soil_velo( 3 );

    noalias( pile_velo ) = ZeroVector( 3 );
    noalias( soil_velo ) = ZeroVector( 3 );

    for ( IndexType i = 0 ; i < mpPileElement->GetGeometry().size() ; i++ )
    {
        double shape_func = mpPileElement->GetGeometry().ShapeFunctionValue( i, mPileLocalPoint );
        pile_velo += (( mpPileElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT ) ) *
                     shape_func;
    }

    for ( IndexType i = 0 ; i < mpSoilElement->GetGeometry().size() ; i++ )
    {
        double shape_func = mpSoilElement->GetGeometry().ShapeFunctionValue( i, mSoilLocalPoint );
        soil_velo += (( mpSoilElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT ) ) *
                     shape_func;
    }

    result( 0 ) = ( pile_velo( 0 ) - soil_velo( 0 ) );

    result( 1 ) = ( pile_velo( 1 ) - soil_velo( 1 ) );

    result( 2 ) = ( pile_velo( 2 ) - soil_velo( 2 ) );

    return result;
}

/**
* Calculates global coordinates of a local point in a reference element
* @param rElement surface
* @param rResult global coordinates
* @param LocalCoordinates local coordinates
* @return global coordinates
*/
Point<3>& Pile_Kinematic_Linear::GetGlobalCoordinates( Element::Pointer rElement, Point<3>& rResult, const Point<3>& LocalCoordinates )
{
    noalias( rResult ) = ZeroVector( 3 );

    for ( IndexType i = 0 ; i < rElement->GetGeometry().size() ; i++ )
    {
        double shape_func = rElement->GetGeometry().ShapeFunctionValue( i, LocalCoordinates );

        rResult( 0 ) += shape_func *
                        (( rElement->GetGeometry()[i] ).X0()
                         + ( rElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_X ) );

        rResult( 1 ) += shape_func *
                        (( rElement->GetGeometry()[i] ).Y0()
                         + ( rElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_Y ) );

        rResult( 2 ) += shape_func *
                        (( rElement->GetGeometry()[i] ).Z0()
                         + ( rElement->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_Z ) );
    }

    return rResult;
}

Matrix Pile_Kinematic_Linear::CalculateTangentVectors( const Kratos::Element::Pointer rElement, const Kratos::Matrix& DN )
{
    Matrix result;
    return result;
}

/// Print information about this object.
void Pile_Kinematic_Linear::PrintInfo( std::ostream& rOStream ) const
{
    rOStream << "Condition #" << Id();
}

/// Print object's data.
void Pile_Kinematic_Linear::PrintData( std::ostream& rOStream ) const
{
    rOStream << "Pile Condition" << std::endl;
}
}
