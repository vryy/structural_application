/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).
1
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
//   Project Name:        Kratos
//   Original Author:     $Author: janosch $
//   Date:                $Date: 2009-01-14 17:14:42 $
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 16 Feb 2021 $
//   Revision:            $Revision: 1.3 $
//

// System includes


// External includes


// Project includes
#include "custom_elements/crisfield_truss_element.h"
#include "structural_application_variables.h"

namespace Kratos
{

/**
* Constructor.
* This deals without DOFs
* @param NewId element ID
* @param pGeometry geometry pointer
* @see CrisfieldTrussElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
*/
CrisfieldTrussElement::CrisfieldTrussElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

/**
* Constructor.
* This deals with DOFs
* @param NewId element ID
* @param pGeometry geometry pointer
* @param pProperties properties pointer
* @see CrisfieldTrussElement(IndexType NewId, GeometryType::Pointer pGeometry)
*/
CrisfieldTrussElement::CrisfieldTrussElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

/**
* Destructor.
*/
CrisfieldTrussElement::~CrisfieldTrussElement()
{
}

/**
* Create a new Crisfield truss element.
* @return crisfield truss elment
* @param NewId element ID
* @param ThisNodes array of nodes
* @param pProperties properties pointer
*/
Element::Pointer CrisfieldTrussElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new CrisfieldTrussElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

/**
* Initialization of the Crisfield truss element.
* This initializes the cross-section, length, position vector and matrix A for the element
*/
void CrisfieldTrussElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

/**
* Calculation of the local system.
* This calculates both the elemental stiffness matrix and the elemental residual vector
* @param rLeftHandSideMatrix elemental stiffness matrix
* @param rRightHandSideVector elemental residual vetor
* @param rCurrentProcessInfo process info
* @see CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
* @see CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
*/
void CrisfieldTrussElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag, 0);
}

/**
* Calculation of the left hand side.
* This calculates only the elemental stiffness matrix
* @param rLeftHandSideMatrix elemental stiffness matrix
* @param rCurrentProcessInfo process info
* @see CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
* @see CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
*/
void CrisfieldTrussElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = false;
    VectorType temp = Vector();

    CalculateAll(rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag, 0);
}

/**
* Calculation of the right hand side.
* This calculates only the elemental residual vector
* @param rRightHandSideVector elemental residual vetor
* @param rCurrentProcessInfo process info
* @see CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
* @see CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
*/
void CrisfieldTrussElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag, 0);
}

/**
* Get the equation ID vector of the element.
* @param rResult equation ID vector
* @param rCurrentProcessInfo process info
*/
void CrisfieldTrussElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dof = number_of_nodes*dimension;

    if(rResult.size() != dof)
        rResult.resize(dof);

    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        int index = i*dimension;
        rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        if(dimension == 3)
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }
}

/**
* Get the DOF list of the element.
* @param ElementalDofList elemental DOF vector
* @param rCurrentProcessInfo process info
*/
void CrisfieldTrussElement::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();

    ElementalDofList.resize(0);

    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        if(dimension == 3)
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }
}

/**
* Get the mass matrix of the element.
* @param rMassMatrix mass matrix
* @param rCurrentProcessInfo process info
*/
void CrisfieldTrussElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    int MatSize = number_of_nodes * dimension;

    if(rMassMatrix.size1() != MatSize || rMassMatrix.size2() != MatSize)
        rMassMatrix.resize(MatSize, MatSize, false);
    noalias(rMassMatrix) = ZeroMatrix(MatSize, MatSize); //resetting LHS

    double rho = GetProperties()[DENSITY];
    double A = GetProperties()[CROSS_AREA];
    double Length = CalculateLength();
    double c = rho*A*Length / 6;

    for (unsigned int i = 0; i < dimension; ++i)
    {
        rMassMatrix(i, i) = 2*c;
        rMassMatrix(i+dimension, i+dimension) = 2*c;
        rMassMatrix(i, i+dimension) = c;
        rMassMatrix(i+dimension, i) = c;
    }

    KRATOS_CATCH("")
}

/**
* Get the damping matrix of the element.
* @param rDampMatrix mass matrix
* @param rCurrentProcessInfo process info
*/
void CrisfieldTrussElement::CalculateDampingMatrix(MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    double alpha = GetProperties()[RAYLEIGH_DAMPING_ALPHA];
    double beta = GetProperties()[RAYLEIGH_DAMPING_BETA];

    this->CalculateMassMatrix(rDampMatrix, rCurrentProcessInfo);

    rDampMatrix *= alpha;

    // compute the stiffness of the current step, to avoid linearization of the stiffness
    Matrix K;
    Vector dummy;
    this->CalculateAll(K, dummy, rCurrentProcessInfo, true, false, 1);

    noalias(rDampMatrix) += beta*K;
}

/**
* Get the displacement vector of the element
* @param values displacement vector
* @param Step solution step
* @see GetFirstDerivativesVector(Vector& values, int Step)
* @see GetSecondDerivativesVector(Vector& values, int Step)
*/
void CrisfieldTrussElement::GetValuesVector(Vector& values, int Step) const
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * dimension;
    if(values.size() != MatSize)
        values.resize(MatSize);
    for ( unsigned int i=0; i<number_of_nodes; i++)
    {
        int index = i*dimension;
        values[index] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X, Step);
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y, Step);
        if(dimension == 3)
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z, Step);
    }
    //KRATOS_WATCH( values );
}

/**
* Get the velocity vector of the element
* @param values velocity vector
* @param Step solution step
* @see GetValuesVector(Vector& values, int Step)
* @see GetSecondDerivativesVector(Vector& values, int Step)
*/
void CrisfieldTrussElement::GetFirstDerivativesVector(Vector& values, int Step) const
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * dimension;
    if(values.size() != MatSize)
        values.resize(MatSize);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        int index = i*dimension;
        values[index] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_DT_X, Step);
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_DT_Y, Step);
        if(dimension == 3)
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_DT_Z, Step);
    }
}

/**
* Get the acceleration vector of the element
* @param values acceleration vector
* @param Step solution step
* @see GetValuesVector(Vector& values, int Step)
* @see GetFirstDerivativesVector(Vector& values, int Step)
*/
void CrisfieldTrussElement::GetSecondDerivativesVector(Vector& values, int Step) const
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * dimension;
    if(values.size() != MatSize)
        values.resize(MatSize);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        int index = i*dimension;
        values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X, Step);
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y, Step);
        if(dimension == 3)
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z, Step);
    }
}

/*
void CrisfieldTrussElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo)
{
    if( rVariable == TRUSS_STRAIN )
    {
        if(Output.size() != 1)
                    Output.resize(1,false);
        Output[0] = msStrain;
    }
    if( rVariable == TRUSS_STRESS )
    {
        if(Output.size() != 1)
                    Output.resize(1,false);
        Output[0] = msStress;
    }
}*/

/**
* Auxiliary function.
* This calculates the elemental stiffness matrix and the elemental residual vector, when required
* @param rLeftHandSideMatrix elemental stiffness matrix
* @param rRightHandSideVector elemental residual vector
* @param rCurrentProcessInfo process info
* @param CalculateStiffnessMatrixFlag flag for elemental stiffness matrix
* @param CalculateResidualVectorFlag flag for elemental residual vector
*/
void CrisfieldTrussElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag,
        const int& index)
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * dimension;

    //resizing the LHS tangent stiffness matrix if required
    if (CalculateStiffnessMatrixFlag == true)
    {
        if(rLeftHandSideMatrix.size1() != MatSize || rLeftHandSideMatrix.size2() != MatSize)
            rLeftHandSideMatrix.resize(MatSize, MatSize, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize); //resetting LHS
    }

    //resizing the RHS residual vector if required
    if (CalculateResidualVectorFlag == true)
    {
        if(rRightHandSideVector.size() != MatSize)
            rRightHandSideVector.resize(MatSize, false);
        rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
    }

    Matrix A(MatSize, MatSize);
    Vector X(MatSize);
    Vector U(MatSize);

    CalculateA(A);
    CalculateX(X);
    CalculateU(U, index);

    //calculation of GREEN-LAGRANGE strain
    double Length = CalculateLength();
    double Area = GetProperties()[CROSS_AREA];
    double weight_strain = 1.0 / pow(Length, 2);
    double Strain = CalculateStrain(A, X, U, weight_strain);
    // KRATOS_WATCH( Strain )

    //Calculation of 2nd-Piola-Kirchhoff stress

    //** Method 1: according to St. Venant Model
    double E = GetProperties()[YOUNG_MODULUS];
    double Stress = E * Strain;
    if (this->Has(PRESTRESS)) // enable prestress (i.e. cable element)
    {
        const Vector& prestress = this->GetValue(PRESTRESS);
        Stress += prestress[0];
    }
    // KRATOS_WATCH( Stress )

    //calculation of the residual force vector if required
    if (CalculateResidualVectorFlag == true)
    {
        CalculateAndAdd_ExtForce(rRightHandSideVector, rCurrentProcessInfo);

        double weight_IntForce = ( Area / Length) * Stress;
    //KRATOS_WATCH(Stress)
    //KRATOS_WATCH(Length)
        CalculateAndMinus_IntForce(rRightHandSideVector, rCurrentProcessInfo, X, U, weight_IntForce);

//          KRATOS_WATCH( rRightHandSideVector );
    }

    //calculation of the tangent stiffness matrix if required
    //KRATOS_WATCH(Area)
    if (CalculateStiffnessMatrixFlag == true)
    {
        double weight_Km = E * Area / pow(Length, 3);
        CalculateAndAddKm(rLeftHandSideMatrix, A, X, U, weight_Km);

        double weight_Kg = Area * Stress / Length ;
        CalculateAndAddKg(rLeftHandSideMatrix, A, weight_Kg);

//          KRATOS_WATCH( rLeftHandSideMatrix );
    }

    KRATOS_CATCH("")
}

/**
* Auxiliary function.
* This adds the external force vector to the elemental residual vector
* @param rRightHandSideVector elemental residual vector
* @param CurrentProcessInfo process info
*/
void CrisfieldTrussElement::CalculateAndAdd_ExtForce(VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo) const
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        int index = dimension*i;
        for ( int j = 0; j < dimension; ++j )
        {
            rRightHandSideVector[index+j] += GetProperties()[BODY_FORCE][j];
        }
    }

    KRATOS_CATCH("")
}

/**
* Auxiliary function.
* This substracts the internal force vector to the elemental residual vector
* @param rRightHandSideVector elemental residual vector
* @param CurrentProcessInfo process info
* @param X position vector
* @param U displacement vector
* @param weight weighting factor
*/
void CrisfieldTrussElement::CalculateAndMinus_IntForce(VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo, Vector& X, Vector& U, double weight) const
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    for ( int i = 0; i < dimension; ++i )
    {
        rRightHandSideVector[i] -= weight*( X[i]+U[i]-X[dimension+i]-U[dimension+i] );
        rRightHandSideVector[dimension+i] -= rRightHandSideVector[i];
    }

    KRATOS_CATCH("")
}

/**
* Auxiliary function.
* This add the material element stiffness matrix to the elemental stiffness matrix
* @param rLeftHandSideMatrix elemental stiffness matrix
* @param A matrix A
* @param X position vector
* @param U displacement vector
* @param weight weighting factor
*/
void CrisfieldTrussElement::CalculateAndAddKm(MatrixType& rLeftHandSideMatrix, const Matrix& A, Vector& X, Vector& U, double weight) const
{
    KRATOS_TRY

    Vector V(X.size());
    noalias(V) = prod(A,(X+U));
    noalias(rLeftHandSideMatrix) += weight * outer_prod(V, V);

    KRATOS_CATCH("")
}

/**
* Auxiliary function.
* This add the geometrial element stiffness matrix to the elemental stiffness matrix
* @param rLeftHandSideMatrix elemental stiffness matrix
* @param A matrix A
* @param weight weighting factor
*/
void CrisfieldTrussElement::CalculateAndAddKg(MatrixType& rLeftHandSideMatrix, const Matrix& A, double weight) const
{
    KRATOS_TRY

    noalias(rLeftHandSideMatrix) += weight * A;

    KRATOS_CATCH("")
}

/**
* Auxiliary function.
* This calculates the GREEN-LAGRANGE strain.
* @return GREEN-LAGRANGE strain
* @param A matrix A
* @param X position vector
* @param U displacement vector
* @param weight weighting factor
*/
double CrisfieldTrussElement::CalculateStrain(const Matrix& A, const Vector& X, const Vector& U, double weight) const
{
    KRATOS_TRY

    Vector V = ZeroVector(X.size());
    noalias(V) = prod(A,U);
    double strain = 0;
    for ( unsigned int i = 0; i < X.size(); ++i )
    {
        strain += (X[i]+0.5*U[i])*V[i];
    }

    return strain * weight ;

    KRATOS_CATCH("")
}

/**
 * Calculate the matrix A
 */
void CrisfieldTrussElement::CalculateA(Matrix& A) const
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dof = number_of_nodes * dimension;
    noalias(A) = ZeroMatrix(dof, dof);
    for (unsigned int i = 0 ; i < dof ; i++)
    {
        for(unsigned int j = 0 ; j < dof ; j++)
        {
            if(i==j)
                A(i,j) = 1;
            if(abs(static_cast<int>(i-j))==dimension)
                A(i,j) = -1;
        }
    }
}

/**
 * Calculate the vector X, i.e. elemental position vector
 */
void CrisfieldTrussElement::CalculateX(Vector& X) const
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    //position vector X (in ref. config)
    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        int index= i*dimension;
        X[index] = GetGeometry()[i].X0();
        X[index+1] = GetGeometry()[i].Y0();
        if(dimension==3)
            X[index+2] = GetGeometry()[i].Z0();
    }
}

/**
 * Calculate the vector U, i.e. elemental deformation vector
 */
void CrisfieldTrussElement::CalculateU(Vector& U, const int& solution_index) const
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    //setting displacement vector U
    for ( unsigned int i = 0; i < number_of_nodes; ++i)
    {
        int index = i*dimension;
        U[index] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X, solution_index);
        U[index+1] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y, solution_index);
        if(dimension == 3)
            U[index+2] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z, solution_index);
    }
}

/**
 * Calculate the length in current configuration
 */
double CrisfieldTrussElement::CalculateLength() const
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //length Length (in ref. config)
    double Length = 0.0;
    Length += pow(GetGeometry()[1].X0()-GetGeometry()[0].X0(), 2);
    Length += pow(GetGeometry()[1].Y0()-GetGeometry()[0].Y0(), 2);
    if(dimension == 3)
        Length += pow(GetGeometry()[1].Z0()-GetGeometry()[0].Z0(), 2);
    return sqrt( Length );
}

void CrisfieldTrussElement::ApplyPrescribedDofs(const MatrixType& LHS_Contribution, VectorType& RHS_Constribution, const ProcessInfo& CurrentProcessInfo) const
{
    // modify the right hand side to account for prescribed displacement
    // according to the book of Bazant & Jirasek, this scheme is more stable than the total displacement scheme for prescribing displacement.
    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int mat_size = dim * GetGeometry().size();
    for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
    {
        if(GetGeometry()[node].IsFixed(DISPLACEMENT_X))
        {
            double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
            if (temp != 0.0)
                for( unsigned int i = 0; i < mat_size; ++i )
                    RHS_Constribution[i] -= LHS_Contribution(i, node * dim) * temp;
        }

        if(GetGeometry()[node].IsFixed(DISPLACEMENT_Y))
        {
            double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
            if (temp != 0.0)
                for( unsigned int i = 0; i < mat_size; ++i )
                    RHS_Constribution[i] -= LHS_Contribution(i, node * dim + 1) * temp;
        }

        if (dim > 2)
        {
            if(GetGeometry()[node].IsFixed(DISPLACEMENT_Z))
            {
                double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                if (temp != 0.0)
                    for( unsigned int i = 0; i < mat_size; ++i )
                        RHS_Constribution[i] -= LHS_Contribution(i, node * dim + 2) * temp;
            }
        }
    }
}

void CrisfieldTrussElement::ComputePrescribedForces(const MatrixType& LHS_Contribution, VectorType& Force, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int mat_size = dim * GetGeometry().size();
    for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
    {
        if(GetGeometry()[node].IsFixed(DISPLACEMENT_X))
        {
            double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DISPLACEMENT_X);
            if (temp != 0.0)
                for( unsigned int i = 0; i < mat_size; ++i )
                    Force[i] -= LHS_Contribution(i, node * dim) * temp;
        }

        if(GetGeometry()[node].IsFixed(DISPLACEMENT_Y))
        {
            double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DISPLACEMENT_Y);
            if (temp != 0.0)
                for( unsigned int i = 0; i < mat_size; ++i )
                    Force[i] -= LHS_Contribution(i, node * dim + 1) * temp;
        }

        if (dim > 2)
        {
            if(GetGeometry()[node].IsFixed(DISPLACEMENT_Z))
            {
                double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DISPLACEMENT_Z);
                if (temp != 0.0)
                    for( unsigned int i = 0; i < mat_size; ++i )
                        Force[i] -= LHS_Contribution(i, node * dim + 2) * temp;
            }
        }
    }
}

int CrisfieldTrussElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    //verify that the variables are correctly initialized

    if ( VELOCITY.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" );

    if ( DISPLACEMENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" );

    if ( ACCELERATION.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" );

    if ( DENSITY.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "" );

    if ( BODY_FORCE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "BODY_FORCE has Key zero! (check if the application is correctly registered", "" );

    if ( AREA.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "AREA has Key zero! (check if the application is correctly registered", "" );

    if ( this->GetProperties().Has( BODY_FORCE ) == false )
        KRATOS_THROW_ERROR( std::logic_error, "BODY_FORCE not provided for property ", this->GetProperties().Id())

    if ( this->GetProperties().Has( CROSS_AREA ) == false )
        KRATOS_THROW_ERROR( std::logic_error, "CROSS_AREA not provided for property ", this->GetProperties().Id())

    return 0;

    KRATOS_CATCH(" ")

}

} // Namespace Kratos
