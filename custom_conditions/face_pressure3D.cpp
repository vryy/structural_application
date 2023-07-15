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
/* **************************************************************************************
*
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2008-10-23 12:26:25 $
*   Revision:            $Revision: 1.1 $
*
* ***************************************************************************************/


// System includes


// External includes


// Project includes
#include "utilities/math_utils.h"
#include "custom_conditions/face_pressure3D.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"

namespace Kratos
{

//***********************************************************************************
//***********************************************************************************
// -------- //
//  PUBLIC  //
// -------- //

// Constructor
FacePressure3D::FacePressure3D()
{
}

// Constructor
FacePressure3D::FacePressure3D( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();//default method
}

// Constructor
FacePressure3D::FacePressure3D( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )

{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();//default method
}

//***********************************************************************************
//***********************************************************************************

Condition::Pointer FacePressure3D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new FacePressure3D( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Condition::Pointer FacePressure3D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new FacePressure3D( NewId, pGeom, pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
// Destructor
FacePressure3D::~FacePressure3D()
{
}

//***********************************************************************************
//***********************************************************************************

void FacePressure3D::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = number_of_nodes * 3;

    if ( rResult.size() != dim )
        rResult.resize( dim );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * 3;
        rResult[index]   = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void FacePressure3D::GetDofList(
    DofsVectorType& ElementalDofList,
    const ProcessInfo& rCurrentProcessInfo ) const
{
    ElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

//***********************************************************************************
//***********************************************************************************


//***********************************************************************************
//***********************************************************************************

void FacePressure3D::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo )

{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//***********************************************************************************
//***********************************************************************************

void FacePressure3D::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo )

{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

}

//***********************************************************************************
//***********************************************************************************

void FacePressure3D::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo )

{
    KRATOS_TRY

    rMassMatrix.resize( 0, 0, false );

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void FacePressure3D::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo )

{
    KRATOS_TRY

    rDampingMatrix.resize( 0, 0, false );

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void FacePressure3D::GetValuesVector(
    Vector& values,
    int Step ) const
{
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;

    if ( values.size() != MatSize )
        values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT, Step );
        unsigned int index = i * 3;
        values[index]   = disp[0];
        values[index+1] = disp[1];
        values[index+2] = disp[2];
    }
}

//***********************************************************************************
//***********************************************************************************

void FacePressure3D::GetFirstDerivativesVector(
    Vector& values,
    int Step ) const
{
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;

    if ( values.size() != MatSize )
        values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3>& vel = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY, Step );
        unsigned int index = i * 3;
        values[index]   = vel[0];
        values[index+1] = vel[1];
        values[index+2] = vel[2];
    }
}

//***********************************************************************************
//***********************************************************************************

void FacePressure3D::GetSecondDerivativesVector(
    Vector& values,
    int Step ) const
{
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;

    if ( values.size() != MatSize )
        values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3>& acc = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION, Step );
        unsigned int index = i * 3;
        values[index]   = acc[0];
        values[index+1] = acc[1];
        values[index+2] = acc[2];
    }
}

//***********************************************************************************
//***********************************************************************************
// --------- //
//  PRIVATE  //
// --------- //


//***********************************************************************************
//***********************************************************************************

void FacePressure3D::CalculateAndSubKp(
    Matrix& K,
    array_1d<double, 3>& ge,
    array_1d<double, 3>& gn,
    const Matrix& DN_De,
    const Vector& N,
    double pressure,
    double weight )

{
    KRATOS_TRY

    boost::numeric::ublas::bounded_matrix<double, 3, 3> Kij;
    boost::numeric::ublas::bounded_matrix<double, 3, 3> Cross_ge;
    boost::numeric::ublas::bounded_matrix<double, 3, 3> Cross_gn;
    double coeff;
    unsigned int number_of_nodes = GetGeometry().size();

    MakeCrossMatrix( Cross_ge, ge );
    MakeCrossMatrix( Cross_gn, gn );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int RowIndex = i * 3;

        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            unsigned int ColIndex = j * 3;

            coeff = pressure * N[i] * DN_De( j, 1 ) * weight;
            noalias( Kij )  = coeff * Cross_ge;

            coeff = pressure * N[i] * DN_De( j, 0 ) * weight;

            noalias( Kij ) -= coeff * Cross_gn;

            //TAKE CARE: the load correction matrix should be SUBTRACTED not added
            SubtractMatrix( K, Kij, RowIndex, ColIndex );
        }
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void FacePressure3D::MakeCrossMatrix(
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& M,
    array_1d<double, 3>& U )

{
    M( 0, 0 ) =  0.00;
    M( 0, 1 ) = -U[2];
    M( 0, 2 ) =  U[1];
    M( 1, 0 ) =  U[2];
    M( 1, 1 ) =  0.00;
    M( 1, 2 ) = -U[0];
    M( 2, 0 ) = -U[1];
    M( 2, 1 ) =  U[0];
    M( 2, 2 ) =  0.00;
}

//***********************************************************************************
//***********************************************************************************

void FacePressure3D::CrossProduct(
    array_1d<double, 3>& cross,
    array_1d<double, 3>& a,
    array_1d<double, 3>& b )

{
    cross[0] = a[1] * b[2] - a[2] * b[1];
    cross[1] = a[2] * b[0] - a[0] * b[2];
    cross[2] = a[0] * b[1] - a[1] * b[0];
}

//***********************************************************************************
//***********************************************************************************

void  FacePressure3D::ExpandReducedMatrix(
    Matrix& Destination,
    Matrix& ReducedMatrix )

{
    KRATOS_TRY

    unsigned int size = ReducedMatrix.size2();

    for ( unsigned int i = 0; i < size; i++ )
    {
        int rowindex = i * 3;

        for ( unsigned int j = 0; j < size; j++ )
        {
            unsigned int colindex = j * 3;

            for ( unsigned int ii = 0; ii < 3; ii++ )
                Destination( rowindex + ii, colindex + ii ) += ReducedMatrix( i, j );
        }
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void  FacePressure3D::SubtractMatrix(
    MatrixType& Destination,
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& InputMatrix,
    int InitialRow,
    int InitialCol )

{
    KRATOS_TRY

    for ( unsigned int i = 0; i < 3; i++ )
        for ( unsigned int j = 0; j < 3; j++ )
            Destination( InitialRow + i, InitialCol + j ) -= InputMatrix( i, j );

    KRATOS_CATCH( "" )
}




//***********************************************************************************
//***********************************************************************************

void FacePressure3D::CalculateAndAdd_PressureForce(
    VectorType& residualvector,
    const Vector& N,
    const array_1d<double, 3>& v3,
    double pressure,
    double weight,
    const ProcessInfo& rCurrentProcessInfo )

{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = 3 * i;
        double coeff = pressure * N[i] * weight;
        residualvector[index]   += coeff * v3[0];
        residualvector[index+1] += coeff * v3[1];
        residualvector[index+2] += coeff * v3[2];
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void FacePressure3D::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag )

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;
    const unsigned int dim= 3;
    //resizing as needed the LHS

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

    Vector t1(3), t2(3), v3(3), Load(3);

    double P = GetValue( NEGATIVE_FACE_PRESSURE ) - GetValue(POSITIVE_FACE_PRESSURE);

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        noalias(t1) = ZeroVector( 3 );//first tangential vector
        noalias(t2) = ZeroVector( 3 );//second tangential vector
        for ( unsigned int n = 0; n < GetGeometry().size(); ++n )
        {
            t1[0] += GetGeometry().GetPoint( n ).X0() * DN_De[PointNumber]( n, 0 );
            t1[1] += GetGeometry().GetPoint( n ).Y0() * DN_De[PointNumber]( n, 0 );
            t1[2] += GetGeometry().GetPoint( n ).Z0() * DN_De[PointNumber]( n, 0 );
            t2[0] += GetGeometry().GetPoint( n ).X0() * DN_De[PointNumber]( n, 1 );
            t2[1] += GetGeometry().GetPoint( n ).Y0() * DN_De[PointNumber]( n, 1 );
            t2[2] += GetGeometry().GetPoint( n ).Z0() * DN_De[PointNumber]( n, 1 );
        }

        //calculating normal
        v3[0] = t1[1] * t2[2] - t1[2] * t2[1];

        v3[1] = t1[2] * t2[0] - t1[0] * t2[2];

        v3[2] = t1[0] * t2[1] - t1[1] * t2[0];

        //calculating load vector
        noalias(Load) = P * v3;

        double IntegrationWeight = integration_points[PointNumber].Weight();

        // RIGHT HAND SIDE VECTOR
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
                for ( unsigned int i = 0; i < dim; ++i )
                    rRightHandSideVector( prim*dim + i ) +=
                        Ncontainer( PointNumber, prim ) * Load( i ) * IntegrationWeight;
        }
    }

//    KRATOS_WATCH(rRightHandSideVector)

    KRATOS_CATCH( "" )
}

int FacePressure3D::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
