/*
*/
/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 6 Feb 2023 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "includes/process_info_with_dofs.h"
#include "utilities/math_utils.h"
#include "geometries/line_2d_3.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"
#include "custom_conditions/mean_displacement_constraint.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

template<int IDimIndex>
MeanDisplacementConstraint<IDimIndex>::MeanDisplacementConstraint( IndexType NewId,
        GeometryType::Pointer pGeometry) :
    Condition( NewId, pGeometry )
{
}

template<int IDimIndex>
MeanDisplacementConstraint<IDimIndex>::MeanDisplacementConstraint( IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) :
    Condition( NewId, pGeometry, pProperties )
{
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************

template<int IDimIndex>
Condition::Pointer MeanDisplacementConstraint<IDimIndex>::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new MeanDisplacementConstraint(NewId, GetGeometry().Create(ThisNodes),
                                pProperties));
}

template<int IDimIndex>
Condition::Pointer MeanDisplacementConstraint<IDimIndex>::Create( IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
{
     return Condition::Pointer( new MeanDisplacementConstraint(NewId, pGeom,
                                pProperties));
}

/**
 * Destructor. Never to be called manually
 */
template<int IDimIndex>
MeanDisplacementConstraint<IDimIndex>::~MeanDisplacementConstraint()
{
}

//************************************************************************************
//************************************************************************************

template<int IDimIndex>
void MeanDisplacementConstraint<IDimIndex>::CalculateRightHandSide( VectorType& rRightHandSideVector,
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

template<int IDimIndex>
void MeanDisplacementConstraint<IDimIndex>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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

template<int IDimIndex>
void MeanDisplacementConstraint<IDimIndex>::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes + 1;

    // resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize || rLeftHandSideMatrix.size2() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    // resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints();

    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        GetGeometry().ShapeFunctionsLocalGradients();

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    // calculating actual jacobian
    MatrixType DeltaPosition(GetGeometry().size(), 3);
    for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates()
                                        - GetGeometry()[node].GetInitialPosition();

    IntegrationMethod ThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian( J0, ThisIntegrationMethod, DeltaPosition );

    // obtain the Lagrange multiplier value
    const int mindex = this->GetValue(LAGRANGE_MULTIPLIER_INDEX);

    const auto* pprocess_info = dynamic_cast<const ProcessInfoWithDofs<>*>(&rCurrentProcessInfo);
    if (pprocess_info == nullptr)
        KRATOS_ERROR << "ProcessInfoWithDofs is required to use MeanDisplacementConstraint condition";

    const double lambda = pprocess_info->GetDof(LAGRANGE_MULTIPLIER_CONSTRAINT, mindex).GetSolutionStepValue();

    //loop over integration points
    VectorType N(number_of_nodes);
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        noalias(N) = row(Ncontainer, PointNumber);

        double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();

        const double dA = std::sqrt(MathUtils<double>::Det(Matrix(prod(trans(J0[PointNumber]), J0[PointNumber]))));

        if ( CalculateResidualVectorFlag )
        {
            // r_u
            for ( unsigned int prim = 0; prim < number_of_nodes; ++prim )
                rRightHandSideVector( prim ) -= lambda * N[prim] * IntegrationWeight * dA;

            // r_lambda
            double mean_u = 0.0;
            for ( unsigned int prim = 0; prim < number_of_nodes; ++prim )
            {
                const double u = GetGeometry()[prim].GetSolutionStepValue(DISPLACEMENT)[IDimIndex];
                mean_u += N[prim] * u * IntegrationWeight * dA;
            }

            rRightHandSideVector(number_of_nodes) -= mean_u;
        }

        if ( CalculateStiffnessMatrixFlag )
        {
            // K_uu is zero

            // K_u_lambda
            for ( unsigned int prim = 0; prim < number_of_nodes; ++prim )
                rLeftHandSideMatrix( prim, number_of_nodes ) += N[prim] * IntegrationWeight * dA;

            // K_lambda_u
            for ( unsigned int prim = 0; prim < number_of_nodes; ++prim )
                rLeftHandSideMatrix( number_of_nodes, prim ) += N[prim] * IntegrationWeight * dA;

            // K_lambda_lambda is zero
        }
    } // end loop over integration points

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

template<int IDimIndex>
void MeanDisplacementConstraint<IDimIndex>::EquationIdVector( EquationIdVectorType& rResult,
        const ProcessInfo& CurrentProcessInfo ) const
{
    //determining size of DOF list
    rResult.resize(GetGeometry().size()+1, false);
    for( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        if (IDimIndex == 0)
            rResult[i] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        else if (IDimIndex == 1)
            rResult[i] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        else if (IDimIndex == 2)
            rResult[i] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }
    if (!this->Has(LAGRANGE_MULTIPLIER_INDEX))
        KRATOS_ERROR << "LAGRANGE_MULTIPLIER_INDEX is not assigned to MeanDisplacementConstraint";

    const int mindex = this->GetValue(LAGRANGE_MULTIPLIER_INDEX);
    const auto* pprocess_info = dynamic_cast<const ProcessInfoWithDofs<>*>(&CurrentProcessInfo);
    if (pprocess_info == nullptr)
        KRATOS_ERROR << "ProcessInfoWithDofs is required to use MeanDisplacementConstraint condition";
    rResult[GetGeometry().size()] = pprocess_info->GetDof(LAGRANGE_MULTIPLIER_CONSTRAINT, mindex).EquationId();
}

//************************************************************************************
//************************************************************************************

template<int IDimIndex>
void MeanDisplacementConstraint<IDimIndex>::GetDofList( DofsVectorType& ConditionalDofList,
                                        const ProcessInfo& CurrentProcessInfo) const
{
    //determining size of DOF list
    ConditionalDofList.resize(GetGeometry().size()+1);
    for( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        if (IDimIndex == 0)
            ConditionalDofList[i] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        else if (IDimIndex == 1)
            ConditionalDofList[i] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        else if (IDimIndex == 2)
            ConditionalDofList[i] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
    }

    if (!this->Has(LAGRANGE_MULTIPLIER_INDEX))
        KRATOS_ERROR << "LAGRANGE_MULTIPLIER_INDEX is not assigned to MeanDisplacementConstraint";

    const int mindex = this->GetValue(LAGRANGE_MULTIPLIER_INDEX);
    const auto* pprocess_info = dynamic_cast<const ProcessInfoWithDofs<>*>(&CurrentProcessInfo);
    if (pprocess_info == nullptr)
        KRATOS_ERROR << "ProcessInfoWithDofs is required to use MeanDisplacementConstraint condition";
    auto* pprocess_info_no_const = const_cast<ProcessInfoWithDofs<>*>(pprocess_info);
    ConditionalDofList[GetGeometry().size()] = pprocess_info_no_const->pGetDof(LAGRANGE_MULTIPLIER_CONSTRAINT, mindex);
    // ConditionalDofList[GetGeometry().size()] = pprocess_info->pGetDof(LAGRANGE_MULTIPLIER_CONSTRAINT, mindex);
}

//************************************************************************************
//************************************************************************************

template class MeanDisplacementConstraint<0>;
template class MeanDisplacementConstraint<1>;
template class MeanDisplacementConstraint<2>;

} // Namespace Kratos
