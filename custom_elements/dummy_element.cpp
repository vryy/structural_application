//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 28 Mar 2017$
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/dummy_element.h"
#include "structural_application_variables.h"


namespace Kratos
{

//************************************************************************************
//************************************************************************************
DummyElement::DummyElement()
{
}

DummyElement::DummyElement( IndexType NewId,
                              GeometryType::Pointer pGeometry)
: Element( NewId, pGeometry )
{
}

DummyElement::DummyElement( IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Element( NewId, pGeometry, pProperties )
{
}

/**
 * Destructor. Never to be called manually
 */
DummyElement::~DummyElement()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Element::Pointer DummyElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DummyElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer DummyElement::Create(IndexType NewId, GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DummyElement(NewId, pGeom, pProperties));
}

void DummyElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // integration rule
    if(this->Has( INTEGRATION_ORDER ))
    {
        if(this->GetValue(INTEGRATION_ORDER) == 1)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_1;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 2)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_2;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 3)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_3;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 4)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_4;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 5)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_ERROR << "DummyElement does not support for integration order " << this->GetValue(INTEGRATION_ORDER);
    }
    else if(GetProperties().Has( INTEGRATION_ORDER ))
    {
        if(GetProperties()[INTEGRATION_ORDER] == 1)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_1;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 2)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_2;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 3)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_3;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 4)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_4;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 5)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_ERROR << "DummyElement does not support for integration order " << GetProperties()[INTEGRATION_ORDER];
    }
    else
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void DummyElement::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
void DummyElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
void DummyElement::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  const ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    rLeftHandSideMatrix.resize(0, 0, false);
    rRightHandSideVector.resize(0, false);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void DummyElement::EquationIdVector( EquationIdVectorType& rResult,
                                     const ProcessInfo& CurrentProcessInfo) const
{
    rResult.resize(0);
}

//************************************************************************************
//************************************************************************************
void DummyElement::GetDofList( DofsVectorType& ElementalDofList,
                               const ProcessInfo& CurrentProcessInfo) const
{
    ElementalDofList.resize(0);
}

//************************************************************************************
//************************************************************************************
void DummyElement::CalculateOnIntegrationPoints( const Variable<MatrixType>& rVariable, std::vector<MatrixType>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    // TODO
}

//************************************************************************************
//************************************************************************************
void DummyElement::CalculateOnIntegrationPoints( const Variable<VectorType>& rVariable, std::vector<VectorType>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    // TODO
}

//************************************************************************************
//************************************************************************************
void DummyElement::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if (rValues.size() != integration_points.size())
        rValues.resize(integration_points.size());

    if( rVariable == INTEGRATION_POINT_GLOBAL || rVariable == INTEGRATION_POINT_GLOBAL_IN_CURRENT_CONFIGURATION )
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
    else
    {
        const MatrixType& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            noalias(rValues[point]) = ZeroVector(3);
            for(std::size_t i = 0; i < GetGeometry().size(); ++i)
            {
                if (GetGeometry()[i].Has(rVariable))
                {
                    const array_1d<double, 3>& value = GetGeometry()[i].GetSolutionStepValue(rVariable);
                    noalias(rValues[point]) += Ncontainer(point, i) * value;
                }
            }
        }
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif
}

//************************************************************************************
//************************************************************************************
void DummyElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if (rValues.size() != integration_points.size())
        rValues.resize(integration_points.size());

    if( rVariable == JACOBIAN_0 )
    {
        // initializing the Jacobian in the reference configuration
        GeometryType::JacobiansType J0;

        MatrixType DeltaPosition(GetGeometry().size(), 3);

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates()
                            - GetGeometry()[node].GetInitialPosition();

        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        // compute the Jacobian determinant
        for( unsigned int i = 0; i < rValues.size(); ++i )
        {
            if (J0[i].size1() == J0[i].size2())
                rValues[i] = MathUtils<double>::Det(J0[i]);
            else
                rValues[i] = sqrt(MathUtils<double>::Det(Matrix(prod(trans(J0[i]), J0[i]))));
        }
    }
    else if( rVariable == INTEGRATION_WEIGHT )
    {
        for( unsigned int i = 0; i < rValues.size(); ++i )
        {
            rValues[i] = integration_points[i].Weight();
        }
    }
    else
    {
        const MatrixType& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            rValues[point] = 0.0;
            for(std::size_t i = 0; i < GetGeometry().size(); ++i)
            {
                if (GetGeometry()[i].Has(rVariable))
                {
                    const double value = GetGeometry()[i].GetSolutionStepValue(rVariable);
                    rValues[point] += Ncontainer(point, i) * value;
                }
            }
        }
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif
}


//************************************************************************************
//************************************************************************************
int DummyElement::Check( const ProcessInfo& CurrentProcessInfo) const
{
    return 0;
}

} // Namespace Kratos

