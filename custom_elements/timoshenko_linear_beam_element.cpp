//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Jun 2020$
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/timoshenko_linear_beam_element.h"
#include "utilities/math_utils.h"
#include "structural_application.h"

// #define DEBUG_BEAM

namespace Kratos
{

//************************************************************************************
//***** Constructor and Destructor ***************************************************
//************************************************************************************
TimoshenkoLinearBeamElement::TimoshenkoLinearBeamElement()
{
}

TimoshenkoLinearBeamElement::TimoshenkoLinearBeamElement( IndexType NewId,
                              GeometryType::Pointer pGeometry)
: Element( NewId, pGeometry )
{
    if( pGeometry->GetGeometryType() != GeometryData::Kratos_Line2D2
     && pGeometry->GetGeometryType() != GeometryData::Kratos_Line3D2 )
        KRATOS_THROW_ERROR(std::logic_error, "This element only works with 2-node line geometry", "")
}

TimoshenkoLinearBeamElement::TimoshenkoLinearBeamElement( IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Element( NewId, pGeometry, pProperties )
{
    if( pGeometry->GetGeometryType() != GeometryData::Kratos_Line2D2
     && pGeometry->GetGeometryType() != GeometryData::Kratos_Line3D2 )
        KRATOS_THROW_ERROR(std::logic_error, "This element only works with 2-node line geometry", "")
}

/**
 * Destructor. Never to be called manually
 */
TimoshenkoLinearBeamElement::~TimoshenkoLinearBeamElement()
{
}

//********************************************************
//**** Operations ****************************************
//********************************************************

Element::Pointer TimoshenkoLinearBeamElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TimoshenkoLinearBeamElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer TimoshenkoLinearBeamElement::Create(IndexType NewId, GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TimoshenkoLinearBeamElement(NewId, pGeom, pProperties));
}

void TimoshenkoLinearBeamElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rCurrentProcessInfo[RESET_CONFIGURATION] == 0)
    {
        mInitialDisp.resize(GetGeometry().size(), 3, false);
        noalias(mInitialDisp) = ZeroMatrix(GetGeometry().size(), 3);
        mInitialRot.resize(GetGeometry().size(), 3, false);
        noalias(mInitialRot) = ZeroMatrix(GetGeometry().size(), 3);
    }
    else if (rCurrentProcessInfo[RESET_CONFIGURATION] == 1)
    {
        if (mInitialDisp.size1() != GetGeometry().size() || mInitialDisp.size2() != 3)
            mInitialDisp.resize( GetGeometry().size(), 3, false );

        if (mInitialRot.size1() != GetGeometry().size() || mInitialRot.size2() != 3)
            mInitialRot.resize( GetGeometry().size(), 3, false );

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            for ( unsigned int i = 0; i < 3; ++i )
            {
                mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];
                mInitialRot( node, i ) = GetGeometry()[node].GetSolutionStepValue( ROTATION )[i];
            }
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
 * calculates only the RHS vector
 */
void TimoshenkoLinearBeamElement::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
void TimoshenkoLinearBeamElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
//************************************************************************************
void TimoshenkoLinearBeamElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
}

//************************************************************************************
//************************************************************************************
void TimoshenkoLinearBeamElement::CalculateDampingMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if(rCurrentProcessInfo[QUASI_STATIC_ANALYSIS])
    {
        CalculateMassMatrix(rDampMatrix, rCurrentProcessInfo);
    }
    else
    {
        // TODO
        KRATOS_THROW_ERROR(std::logic_error, "TRANSIENT ANALYSIS is not implemented at TimoshenkoLinearBeamElement::", __FUNCTION__)
    }
}

//************************************************************************************
//************************************************************************************
/**
 * This function calculates all system contributions
 */
void TimoshenkoLinearBeamElement::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo,
                                  const bool& CalculateStiffnessMatrixFlag,
                                  const bool& CalculateResidualVectorFlag)
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int mat_size;

    if (dim == 2)
        mat_size = 3*number_of_nodes;
    else if (dim == 3)
        mat_size = 6*number_of_nodes;

    if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);
    noalias( rLeftHandSideMatrix ) = ZeroMatrix(mat_size, mat_size);

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size, false);
    noalias( rRightHandSideVector ) = ZeroVector(mat_size);

    Matrix transformation_matrix(mat_size, mat_size);
    this->CalculateInitialLocalCS(transformation_matrix);

    Matrix material_matrix;
    this->CreateElementStiffnessMatrix_Material(material_matrix);

    if (Id() == 1)
    {
        KRATOS_WATCH(transformation_matrix)
        KRATOS_WATCH(material_matrix)
    }

    noalias(rLeftHandSideMatrix) = prod( Matrix( prod( transformation_matrix, material_matrix ) ), trans(transformation_matrix));

    Vector Disp(mat_size);
    this->GetCurrentDisplacement(Disp);

    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, Disp);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
* Setting up the EquationIdVector
*/
void TimoshenkoLinearBeamElement::EquationIdVector( EquationIdVectorType& rResult,
                                      ProcessInfo& CurrentProcessInfo)
{
    unsigned int dofs_per_node;
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if (dim == 2)
        dofs_per_node = 3;
    else if (dim == 3)
        dofs_per_node = 6;

    unsigned int mat_size = GetGeometry().size() * dofs_per_node;

    if ( rResult.size() != mat_size )
        rResult.resize( mat_size, false );

    if (dim == 2)
    {
        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            int index = i * dofs_per_node;
            rResult[index    ] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof( ROTATION_Z ).EquationId();
        }
    }
    else if (dim == 3)
    {
        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            int index = i * dofs_per_node;
            rResult[index    ] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( ROTATION_X ).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof( ROTATION_Y ).EquationId();
            rResult[index + 5] = GetGeometry()[i].GetDof( ROTATION_Z ).EquationId();
        }
    }
}

//************************************************************************************
//************************************************************************************
/**
 * Setting up the DOF list
 */
void TimoshenkoLinearBeamElement::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    unsigned int dofs_per_node;
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if (dim == 2)
        dofs_per_node = 3;
    else if (dim == 3)
        dofs_per_node = 6;

    unsigned int mat_size = GetGeometry().size() * dofs_per_node;

    if ( ElementalDofList.size() != mat_size )
        ElementalDofList.resize( mat_size );

    if (dim == 2)
    {
        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            int index = i * dofs_per_node;
            ElementalDofList[index    ] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
            ElementalDofList[index + 1] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
            ElementalDofList[index + 2] = GetGeometry()[i].pGetDof( ROTATION_Z );
        }
    }
    else if (dim == 3)
    {
        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            int index = i * dofs_per_node;
            ElementalDofList[index    ] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
            ElementalDofList[index + 1] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
            ElementalDofList[index + 2] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
            ElementalDofList[index + 3] = GetGeometry()[i].pGetDof( ROTATION_X );
            ElementalDofList[index + 4] = GetGeometry()[i].pGetDof( ROTATION_Y );
            ElementalDofList[index + 5] = GetGeometry()[i].pGetDof( ROTATION_Z );
        }
    }
}

//************************************************************************************
//************************************************************************************
void TimoshenkoLinearBeamElement::GetCurrentDisplacement( Vector& Disp ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int dofs_per_node;
    unsigned int number_of_nodes = GetGeometry().size();

    if (dim == 2)
    {
        dofs_per_node = 3;
        if (Disp.size() != number_of_nodes*dofs_per_node)
            Disp.resize(number_of_nodes*dofs_per_node, false);

        for ( unsigned int i = 0 ; i < number_of_nodes ; ++i )
        {
            int index = i * dofs_per_node;
            Disp[index    ] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X ) - mInitialDisp(i, 0);
            Disp[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y ) - mInitialDisp(i, 1);
            Disp[index + 2] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z ) - mInitialRot(i, 2);
        }
    }
    else if (dim == 3)
    {
        dofs_per_node = 6;
        if (Disp.size() != number_of_nodes*dofs_per_node)
            Disp.resize(number_of_nodes*dofs_per_node, false);

        for ( unsigned int i = 0 ; i < number_of_nodes ; ++i )
        {
            int index = i * dofs_per_node;
            Disp[index    ] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X ) - mInitialDisp(i, 0);
            Disp[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y ) - mInitialDisp(i, 1);
            Disp[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z ) - mInitialDisp(i, 2);
            Disp[index + 3] = GetGeometry()[i].GetSolutionStepValue( ROTATION_X ) - mInitialRot(i, 0);
            Disp[index + 4] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y ) - mInitialRot(i, 1);
            Disp[index + 5] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z ) - mInitialRot(i, 2);
        }
    }
}

//************************************************************************************
//************************************************************************************
void TimoshenkoLinearBeamElement::CalculateInitialLocalCS(Matrix& transformation_matrix) const
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int mat_size;

    if (dim == 2)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Not implemented", "")
    }
    else if (dim == 3)
    {
        mat_size = 6*GetGeometry().size();

        if (transformation_matrix.size1() != mat_size || transformation_matrix.size2() != mat_size)
            transformation_matrix.resize(mat_size, mat_size, false);
        noalias(transformation_matrix) = ZeroMatrix(mat_size, mat_size);

        // const double numerical_limit = std::numeric_limits<double>::epsilon();
        // array_1d<double, 3> direction_vector_x = ZeroVector(3);
        // array_1d<double, 3> direction_vector_y = ZeroVector(3);
        // array_1d<double, 3> direction_vector_z = ZeroVector(3);
        // array_1d<double, 6> reference_coordinates = ZeroVector(6);

        // reference_coordinates[0] = GetGeometry()[0].X0();
        // reference_coordinates[1] = GetGeometry()[0].Y0();
        // reference_coordinates[2] = GetGeometry()[0].Z0();
        // reference_coordinates[3] = GetGeometry()[1].X0();
        // reference_coordinates[4] = GetGeometry()[1].Y0();
        // reference_coordinates[5] = GetGeometry()[1].Z0();

        // for (unsigned int i = 0; i < 3; ++i)
        // {
        //     direction_vector_x[i] =
        //         (reference_coordinates[i + 3] - reference_coordinates[i]);
        // }
        // Matrix temp_matrix = ZeroMatrix(3);

        // typedef array_1d<double, 3> arraydim;
        // arraydim global_z = ZeroVector(3);
        // global_z[2] = 1.0;

        // arraydim v2 = ZeroVector(3);
        // arraydim v3 = ZeroVector(3);

        // double vector_norm;
        // vector_norm = MathUtils<double>::Norm(direction_vector_x);
        // if (vector_norm > numerical_limit)
        // {
        //     direction_vector_x /= vector_norm;
        // }

        // if (std::abs(direction_vector_x[2] - 1.00) < numerical_limit)
        // {
        //     v2[1] = 1.0;
        //     v3[0] = -1.0;
        // }
        // else if (std::abs(direction_vector_x[2] + 1.00) < numerical_limit)
        // {
        //     v2[1] = 1.0;
        //     v3[0] = 1.0;
        // }
        // else
        // {
        //     noalias(v2) = MathUtils<double>::UnitCrossProduct(global_z, direction_vector_x);
        //     noalias(v3) = MathUtils<double>::UnitCrossProduct(direction_vector_x, v2);
        // }

        // for (int i = 0; i < 3; ++i) {
        //     temp_matrix(i, 0) = direction_vector_x[i];
        //     temp_matrix(i, 1) = v2[i];
        //     temp_matrix(i, 2) = v3[i];
        // }

        // // Create big rotation Matrix
        // for (unsigned int kk = 0; kk < 12; kk += 3)
        // {
        //     for (int i = 0; i < 3; ++i)
        //     {
        //         for (int j = 0; j < 3; ++j)
        //         {
        //             if (std::abs(temp_matrix(i, j)) <= numerical_limit)
        //             {
        //                 transformation_matrix(i + kk, j + kk) = 0.00;
        //             } else {
        //                 transformation_matrix(i + kk, j + kk) = temp_matrix(i, j);
        //             }
        //         }
        //     }
        // }


        const double numerical_limit = std::numeric_limits<double>::epsilon();
        array_1d<double, 3> direction_vector_x = ZeroVector(3);
        array_1d<double, 3> direction_vector_y = ZeroVector(3);
        array_1d<double, 3> direction_vector_z = ZeroVector(3);
        array_1d<double, 6> reference_coordinates = ZeroVector(6);

        reference_coordinates[0] = GetGeometry()[0].X0();
        reference_coordinates[1] = GetGeometry()[0].Y0();
        reference_coordinates[2] = GetGeometry()[0].Z0();
        reference_coordinates[3] = GetGeometry()[1].X0();
        reference_coordinates[4] = GetGeometry()[1].Y0();
        reference_coordinates[5] = GetGeometry()[1].Z0();

        for (unsigned int i = 0; i < 3; ++i)
        {
            direction_vector_x[i] = reference_coordinates[i + 3] - reference_coordinates[i];
        }
        Matrix temp_matrix = ZeroMatrix(3);

        // use orientation class 1st constructor
        array_1d<double, 3> global_z = ZeroVector(3);
        global_z[2] = 1.0;

        array_1d<double, 3> v2 = ZeroVector(3);
        array_1d<double, 3> v3 = ZeroVector(3);

        double vector_norm;
        vector_norm = MathUtils<double>::Norm(direction_vector_x);
        if (vector_norm > numerical_limit) {
            direction_vector_x /= vector_norm;
        }

        if (std::abs(direction_vector_x[2] - 1.0) < numerical_limit) {
            v2[1] = 1.0;
            v3[0] = -1.0;
        }

        else if (std::abs(direction_vector_x[2] + 1.0) < numerical_limit) {
            v2[1] = 1.0;
            v3[0] = 1.0;
        }

        else {
            noalias(v2) = MathUtils<double>::UnitCrossProduct(global_z, direction_vector_x);
            noalias(v3) = MathUtils<double>::UnitCrossProduct(direction_vector_x, v2);
        }

        for (int i = 0; i < 3; ++i) {
            temp_matrix(i, 0) = direction_vector_x[i];
            temp_matrix(i, 1) = v2[i];
            temp_matrix(i, 2) = v3[i];
        }

        // assemble rotation matrix
        for (int kk = 0; kk < 12; kk += 3) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    if (std::abs(temp_matrix(i, j)) <= numerical_limit) {
                        transformation_matrix(i + kk, j + kk) = 0.0;
                    } else {
                        transformation_matrix(i + kk, j + kk) = temp_matrix(i, j);
                    }
                }
            }
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void TimoshenkoLinearBeamElement::CreateElementStiffnessMatrix_Material(Matrix& local_stiffness_matrix) const
{

    KRATOS_TRY;

    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int mat_size;

    if (dim == 2)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
    }
    else if (dim == 3)
    {
        mat_size = 6*GetGeometry().size();

        if (local_stiffness_matrix.size1() != mat_size || local_stiffness_matrix.size2() != mat_size)
            local_stiffness_matrix.resize(mat_size, mat_size, false);
        noalias(local_stiffness_matrix) = ZeroMatrix(mat_size, mat_size);

        const double E = GetProperties()[YOUNG_MODULUS];
        const double nu = GetProperties()[POISSON_RATIO];
        const double G = E / (2.0 * (1.0 + nu));
        const double A = GetProperties()[AREA];
        const double L = MathUtils<double>::Norm3(GetGeometry()[1].GetInitialPosition().Coordinates()
                - GetGeometry()[0].GetInitialPosition().Coordinates());

        const double Ixx = GetProperties()[INERTIA_X];
        const double Iyy = GetProperties()[INERTIA_Y];
        const double Izz = GetProperties()[INERTIA_Z];

        double Ay = 0.00;
        if (GetProperties().Has(AREA_Y))
        {
            Ay = GetProperties()[AREA_Y];
        }

        double Az = 0.00;
        if (GetProperties().Has(AREA_Z))
        {
            Az = GetProperties()[AREA_Z];
        }

        // axial
        local_stiffness_matrix(0,0) =  E * A / L;
        local_stiffness_matrix(0,6) = -E * A / L;

        local_stiffness_matrix(3,3) =  E * Ixx / L;
        local_stiffness_matrix(3,9) = -E * Ixx / L;

        local_stiffness_matrix(6,6) =  E * A / L;
        local_stiffness_matrix(6,0) = -E * A / L;

        local_stiffness_matrix(9,3) =  E * Ixx / L;
        local_stiffness_matrix(9,9) = -E * Ixx / L;

        // bending
        local_stiffness_matrix(1,1)  =  G * Ay / L;
        local_stiffness_matrix(1,5)  =  G * Ay / 2.0;
        local_stiffness_matrix(1,7)  = -G * Ay / L;
        local_stiffness_matrix(1,11) =  G * Ay / 2.0;

        local_stiffness_matrix(2,2)  =  G * Az / L;
        local_stiffness_matrix(2,4)  =  G * Az / 2.0;
        local_stiffness_matrix(2,8)  = -G * Az / L;
        local_stiffness_matrix(2,10) =  G * Az / 2.0;

        local_stiffness_matrix(4,2)  =  G * Az / 2.0;
        local_stiffness_matrix(4,4)  =  E * Iyy / L + G * Az * L / 4.0;
        local_stiffness_matrix(4,8)  = -G * Az / 2.0;
        local_stiffness_matrix(4,10) = -E * Iyy / L + G * Az * L / 4.0;

        local_stiffness_matrix(5,1)  =  G * Ay / 2.0;
        local_stiffness_matrix(5,5)  =  E * Izz / L + G * Ay * L / 4.0;
        local_stiffness_matrix(5,7)  = -G * Ay / 2.0;
        local_stiffness_matrix(5,11) = -E * Izz / L + G * Ay * L / 4.0;

        local_stiffness_matrix(7,1)  = -G * Ay / L;
        local_stiffness_matrix(7,5)  = -G * Ay / 2.0;
        local_stiffness_matrix(7,7)  =  G * Ay / L;
        local_stiffness_matrix(7,11) = -G * Ay / 2.0;

        local_stiffness_matrix(8,2)  = -G * Az / L;
        local_stiffness_matrix(8,4)  = -G * Az / 2.0;
        local_stiffness_matrix(8,8)  =  G * Az / L;
        local_stiffness_matrix(8,10) = -G * Az / 2.0;

        local_stiffness_matrix(10,2)  =  G * Az / 2.0;
        local_stiffness_matrix(10,4)  = -E * Iyy / L + G * Az * L / 4.0;
        local_stiffness_matrix(10,8)  = -G * Az / 2.0;
        local_stiffness_matrix(10,10) =  E * Iyy / L + G * Az * L / 4.0;

        local_stiffness_matrix(11,1)  =  G * Ay / 2.0;
        local_stiffness_matrix(11,5)  = -E * Izz / L + G * Ay * L / 4.0;
        local_stiffness_matrix(11,7)  = -G * Ay / 2.0;
        local_stiffness_matrix(11,11) =  E * Izz / L + G * Ay * L / 4.0;
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void TimoshenkoLinearBeamElement::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
}

//************************************************************************************
//************************************************************************************
void TimoshenkoLinearBeamElement::GetValueOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
        std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo )
{
}

//************************************************************************************
//************************************************************************************
void TimoshenkoLinearBeamElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    this->GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void TimoshenkoLinearBeamElement::CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
        std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    this->GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

int TimoshenkoLinearBeamElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH("");
}

} // Namespace Kratos

#undef DEBUG_BEAM
