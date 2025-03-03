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
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 22 Feb 2025 $
*   Revision:            $Revision: 1.0 $
*
* ***************************************************************************************/


// System includes


// External includes


// Project includes
#include "custom_conditions/line_pressure_distributed.h"
#include "structural_application_variables.h"

namespace Kratos
{
//***********************************************************************************
//***********************************************************************************
// -------- //
//  PUBLIC  //
// -------- //

// Constructor
LinePressureDistributed::LinePressureDistributed()
{
}

// Constructor
LinePressureDistributed::LinePressureDistributed( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
}

// Constructor
LinePressureDistributed::LinePressureDistributed( IndexType NewId, GeometryType::Pointer pGeometry,
                                                  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LinePressureDistributed::Create( IndexType NewId,
                                                    NodesArrayType const& ThisNodes,
                                                    PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new LinePressureDistributed( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LinePressureDistributed::Create( IndexType NewId,
                                                    GeometryType::Pointer pGeom,
                                                    PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new LinePressureDistributed( NewId, pGeom, pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
// Destructor
LinePressureDistributed::~LinePressureDistributed()
{
}

//***********************************************************************************
//***********************************************************************************
void LinePressureDistributed::Initialize( const ProcessInfo& rCurrentProcessInfo )
{
    const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints( ThisIntegrationMethod );

    if (mIntegrationPointPressure.size() != integration_points.size())
        mIntegrationPointPressure.resize(integration_points.size());

    const double P = this->GetValue(PRESSURE);
    for (std::size_t i = 0; i < mIntegrationPointPressure.size(); ++i)
    {
        mIntegrationPointPressure[i] = P;
    }
}

//***********************************************************************************
//***********************************************************************************
LinePressureDistributed::IntegrationMethod LinePressureDistributed::GetIntegrationMethod() const
{
    if(this->Has( INTEGRATION_ORDER ))
    {
        if(this->GetValue(INTEGRATION_ORDER) == 1)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 2)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 3)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 4)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 5)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_ERROR << Info() << " does not support for integration order " << this->GetValue(INTEGRATION_ORDER);
    }
    else if(GetProperties().Has( INTEGRATION_ORDER ))
    {
        if(GetProperties()[INTEGRATION_ORDER] == 1)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 2)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 3)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 4)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 5)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_ERROR << Info() << " does not support for integration order " << GetProperties()[INTEGRATION_ORDER];
    }
    else
        return GetGeometry().GetDefaultIntegrationMethod(); // default method
}

//***********************************************************************************
//***********************************************************************************
void LinePressureDistributed::EquationIdVector( EquationIdVectorType& rResult,
                                                const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    DofsVectorType ConditionalDofList;
    GetDofList(ConditionalDofList, rCurrentProcessInfo);

    if (rResult.size() != ConditionalDofList.size())
        rResult.resize(ConditionalDofList.size(), false);

    for(unsigned int i = 0; i < ConditionalDofList.size(); ++i)
    {
        rResult[i] = ConditionalDofList[i]->EquationId();
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void LinePressureDistributed::GetDofList( DofsVectorType& ElementalDofList,
                                          const ProcessInfo& rCurrentProcessInfo ) const
{
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    ElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
    }
}

//***********************************************************************************
//***********************************************************************************
void LinePressureDistributed::SetValuesOnIntegrationPoints( const Variable<double>& rVariable,
                                                            const std::vector<double>& rValues,
                                                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == PRESSURE )
    {
        if (mIntegrationPointPressure.size() != rValues.size())
            mIntegrationPointPressure.resize(rValues.size());
            // KRATOS_ERROR << "The input size array is incompatible";
        for (std::size_t i = 0; i < mIntegrationPointPressure.size(); ++i)
        {
            mIntegrationPointPressure[i] = rValues[i];
        }
    }
}

//***********************************************************************************
//***********************************************************************************
void LinePressureDistributed::CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
                                                            std::vector<array_1d<double, 3> >& rValues,
                                                            const ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    GetGeometry().Initialize(ThisIntegrationMethod);
    #endif

    const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints( ThisIntegrationMethod );

    if (rValues.size() != integration_points.size())
        rValues.resize(integration_points.size());

    if( rVariable == DISPLACEMENT )
    {
        const MatrixType& Ncontainer = GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            noalias(rValues[point]) = ZeroVector(3);
            for(std::size_t i = 0; i < GetGeometry().size(); ++i)
            {
                const array_1d<double, 3>& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                noalias(rValues[point]) += Ncontainer(point, i) * displacement;
            }
        }
    }
    else if( rVariable == INTEGRATION_POINT_GLOBAL || rVariable == INTEGRATION_POINT_GLOBAL_IN_CURRENT_CONFIGURATION )
    {
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            rValues[point] = GetGeometry().GlobalCoordinates(rValues[point], integration_points[point]);
        }
    }
    else if( rVariable == INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION )
    {
        VectorType N( GetGeometry().size() );

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            GetGeometry().ShapeFunctionsValues( N, integration_points[point] );

            noalias( rValues[point] ) = ZeroVector(3);
            for(std::size_t i = 0 ; i < GetGeometry().size() ; ++i)
                noalias( rValues[point] ) += N[i] * GetGeometry()[i].GetInitialPosition();
        }
    }
    else if( rVariable == INTEGRATION_POINT_LOCAL )
    {
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            noalias(rValues[point]) = integration_points[point];
        }
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif
}

//***********************************************************************************
//***********************************************************************************
void LinePressureDistributed::CalculateRightHandSide( VectorType& rRightHandSideVector,
                                                      const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    //resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatSize )
        rRightHandSideVector.resize( MatSize, false );
    rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

    const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    GetGeometry().Initialize(ThisIntegrationMethod);
    #endif

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(ThisIntegrationMethod);

    // DN_DeContainer is the array of shape function gradients at each integration points
    const GeometryType::ShapeFunctionsGradientsType& DN_De =
        GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod);

    // Ncontainer is the array of shape function values at each integration points
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(ThisIntegrationMethod);

    // loop over integration points
    Vector Load( dim );
    Vector t( dim );
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        // compute integration weight
        double IntegrationWeight = integration_points[PointNumber].Weight();

        IntegrationWeight *= GetProperties()[THICKNESS];

        // compute tangetial vector (include length)
        noalias( t ) = ZeroVector( dim );//tangential vector
        for ( unsigned int n = 0; n < GetGeometry().size(); ++n )
        {
            t[0] += GetGeometry().GetPoint( n ).X0() * DN_De[PointNumber]( n, 0 );
            t[1] += GetGeometry().GetPoint( n ).Y0() * DN_De[PointNumber]( n, 0 );
        }
//        KRATOS_WATCH(t)

        //calculating load
        Load[0] = -mIntegrationPointPressure[PointNumber] * t[1];
        Load[1] = mIntegrationPointPressure[PointNumber] * t[0];
//        KRATOS_WATCH(Load)

        // contribute to RIGHT HAND SIDE VECTOR
        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
            for ( unsigned int i = 0; i < dim; ++i )
                rRightHandSideVector( prim * dim + i ) +=
                    Ncontainer( PointNumber, prim ) * Load( i ) * IntegrationWeight;
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void LinePressureDistributed::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    const ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    rLeftHandSideMatrix = ZeroMatrix( dim * number_of_nodes, dim * number_of_nodes );
    CalculateRightHandSide( rRightHandSideVector, rCurrentProcessInfo);
}

//***********************************************************************************
//***********************************************************************************
void LinePressureDistributed::CalculateMassMatrix( MatrixType& rMassMatrix,
                                                   const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rMassMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void LinePressureDistributed::CalculateDampingMatrix( MatrixType& rDampingMatrix,
                                                      const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rDampingMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int LinePressureDistributed::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    return 0;
}

} // Namespace Kratos.
