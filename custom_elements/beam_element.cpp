
/* *********************************************************
*
*   Last Modified by:    $Author: nelson, hbui $
*   Date:                $Date: 14/3/2021 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/

// System includes

// External includes

// Project includes
#include "geometries/line_3d_2.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_elements/beam_element.h"
#include "structural_application_variables.h"



namespace Kratos
{
//*****************************************************************************
//*****************************************************************************

typedef GeometryData::IntegrationMethod IntegrationMethod;

BeamElement::BeamElement(IndexType NewId,GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
    //THIS IS THE DEFAULT CONSTRUCTOR
}

BeamElement::BeamElement(IndexType NewId,GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

BeamElement::BeamElement(IndexType NewId, GeometryType::PointType::Pointer pNode1,
        GeometryType::PointType::Pointer pNode2, PropertiesType::Pointer pProperties)
: Element(NewId, GeometryType::Pointer(new Line3D2<GeometryType::PointType>(pNode1, pNode2)), pProperties)
{
}

Element::Pointer BeamElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new BeamElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer BeamElement::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new BeamElement(NewId, pGeom, pProperties));
}

BeamElement::~BeamElement()
{
}

//************************************************************************************
//THIS IS THE INITIALIZATION OF THE ELEMENT (CALLED AT THE BEGIN OF EACH CALCULATION)
//************************************************************************************


void BeamElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
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

    mPreForces.resize(12, false);
    noalias(mPreForces) = ZeroVector(12);

    mCurrentForces.resize(12, false);
    noalias(mCurrentForces) = ZeroVector(12);

    CalculateSectionProperties();

    KRATOS_CATCH("")
}

void BeamElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    noalias(mCurrentForces) = mPreForces;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void BeamElement::InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo)
{

}
//************************************************************************************
//************************************************************************************

void BeamElement::FinalizeSolutionStep(const ProcessInfo& CurrentProcessInfo)
{
}


//************************************************************************************
//************************************************************************************


void BeamElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo,
                               bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (dimension != 3)
    {
        std::cout<<"this element works only with a 2 node line and 3D dimension"<<std::endl;
        return;
    }

    if (CalculateStiffnessMatrixFlag)
    {
        CalculateLHS(rLeftHandSideMatrix);
    }

    if (CalculateResidualVectorFlag)
    {
        CalculateRHS(rRightHandSideVector);
    }

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void  BeamElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************


void BeamElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                       VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,CalculateResidualVectorFlag);
    // KRATOS_WATCH(Id())
    // KRATOS_WATCH(rLeftHandSideMatrix)
    // KRATOS_WATCH(rRightHandSideVector)
}

//************************************************************************************
//************************************************************************************

void BeamElement::EquationIdVector(EquationIdVectorType& rResult,
                                   const ProcessInfo& CurrentProcessInfo) const
{
    if(rResult.size() != 12)
        rResult.resize(12,false);

    rResult[0]  = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1]  = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[2]  = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
    rResult[3]  = GetGeometry()[0].GetDof(ROTATION_X).EquationId();
    rResult[4]  = GetGeometry()[0].GetDof(ROTATION_Y).EquationId();
    rResult[5]  = GetGeometry()[0].GetDof(ROTATION_Z).EquationId();
    rResult[6]  = GetGeometry()[1].GetDof(DISPLACEMENT_X).EquationId();
    rResult[7]  = GetGeometry()[1].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[8]  = GetGeometry()[1].GetDof(DISPLACEMENT_Z).EquationId();
    rResult[9]  = GetGeometry()[1].GetDof(ROTATION_X).EquationId();
    rResult[10] = GetGeometry()[1].GetDof(ROTATION_Y).EquationId();
    rResult[11] = GetGeometry()[1].GetDof(ROTATION_Z).EquationId();
}

//************************************************************************************
//************************************************************************************

void BeamElement::GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo&
                             CurrentProcessInfo) const
{
    ElementalDofList.resize(0);

    ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_X));
    ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Y));
    ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Z));
    ElementalDofList.push_back(GetGeometry()[0].pGetDof(ROTATION_X));
    ElementalDofList.push_back(GetGeometry()[0].pGetDof(ROTATION_Y));
    ElementalDofList.push_back(GetGeometry()[0].pGetDof(ROTATION_Z));
    ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_X));
    ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Y));
    ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Z));
    ElementalDofList.push_back(GetGeometry()[1].pGetDof(ROTATION_X));
    ElementalDofList.push_back(GetGeometry()[1].pGetDof(ROTATION_Y));
    ElementalDofList.push_back(GetGeometry()[1].pGetDof(ROTATION_Z));
}

//************************************************************************************
//************************************************************************************

void BeamElement::CalculateLHS(Matrix& rLeftHandSideMatrix)
{
    Matrix LocalMatrix;
    Matrix Rotation;
    Matrix aux_matrix;

    LocalMatrix.resize(12, 12, false);
    Rotation.resize(12, 12, false);
    aux_matrix.resize(12, 12, false);

    if (rLeftHandSideMatrix.size1() != 12 || rLeftHandSideMatrix.size2() != 12)
        rLeftHandSideMatrix.resize(12, 12, false);

    CalculateLocalMatrix(LocalMatrix);
    CalculateTransformationMatrix(Rotation);
    noalias(aux_matrix) = prod(Rotation, LocalMatrix);
    noalias(rLeftHandSideMatrix)= prod(aux_matrix, Matrix(trans(Rotation)));

    /// adding the contribution from rotational stiffness if needed
    if (GetProperties().Has(ROTATIONAL_STIFFNESS))
    {
        const double Kr = GetProperties()[ROTATIONAL_STIFFNESS];

        array_1d<double, 3> axialVector;
        axialVector[0] = GetGeometry()[1].X0() - GetGeometry()[0].X0();
        axialVector[1] = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
        axialVector[2] = GetGeometry()[1].Z0() - GetGeometry()[0].Z0();
        double length = norm_2(axialVector);
        double aux = 0.25*Kr*length;

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                rLeftHandSideMatrix(i+3, j+3) += aux*axialVector[i]*axialVector[j];
                rLeftHandSideMatrix(i+3, j+6) += aux*axialVector[i]*axialVector[j];
                rLeftHandSideMatrix(i+6, j+3) += aux*axialVector[i]*axialVector[j];
                rLeftHandSideMatrix(i+6, j+6) += aux*axialVector[i]*axialVector[j];
            }
        }
    }

    return;
}


//************************************************************************************
//************************************************************************************
void BeamElement::CalculateRHS(Vector& rRightHandSideVector)
{
    Matrix Rotation;
    Matrix GlobalMatrix;
    Vector LocalBody;

    array_1d<double, 12> CurrentDisplacement;

    Rotation.resize(12, 12, false);
    LocalBody = ZeroVector(12);

    if (rRightHandSideVector.size() != 12)
        rRightHandSideVector.resize(12, false);
    noalias(rRightHandSideVector) = ZeroVector(12);

    CalculateTransformationMatrix(Rotation);
    CalculateBodyForce(Rotation, LocalBody, rRightHandSideVector);

    CurrentDisplacement(0)      =   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X) - mInitialDisp(0, 0);
    CurrentDisplacement(1)      =   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y) - mInitialDisp(0, 1);
    CurrentDisplacement(2)      =   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z) - mInitialDisp(0, 2);
    CurrentDisplacement(3)      =   GetGeometry()[0].GetSolutionStepValue(ROTATION_X) - mInitialRot(0, 0);
    CurrentDisplacement(4)      =   GetGeometry()[0].GetSolutionStepValue(ROTATION_Y) - mInitialRot(0, 1);
    CurrentDisplacement(5)      =   GetGeometry()[0].GetSolutionStepValue(ROTATION_Z) - mInitialRot(0, 2);
    CurrentDisplacement(6)      =   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X) - mInitialDisp(1, 0);
    CurrentDisplacement(7)      =   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y) - mInitialDisp(1, 1);
    CurrentDisplacement(8)      =   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z) - mInitialDisp(1, 2);
    CurrentDisplacement(9)      =   GetGeometry()[1].GetSolutionStepValue(ROTATION_X) - mInitialRot(1, 0);
    CurrentDisplacement(10)     =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Y) - mInitialRot(1, 1);
    CurrentDisplacement(11)     =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Z) - mInitialRot(1, 2);

    CalculateLHS(GlobalMatrix);
    noalias(mCurrentForces) = prod(GlobalMatrix, CurrentDisplacement);
    // KRATOS_WATCH(GlobalMatrix)
    // KRATOS_WATCH(CurrentDisplacement)
    // KRATOS_WATCH(GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT))
    // KRATOS_WATCH(GetGeometry()[0].GetSolutionStepValue(ROTATION))
    // KRATOS_WATCH(GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT))
    // KRATOS_WATCH(GetGeometry()[1].GetSolutionStepValue(ROTATION))
    // KRATOS_WATCH(mInitialDisp)
    // KRATOS_WATCH(mInitialRot)
    noalias(rRightHandSideVector) -= mCurrentForces;
    noalias(rRightHandSideVector) += mPreForces;

    return;
}


//************************************************************************************
//************************************************************************************

void BeamElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag  = false;
    Vector temp = Vector();
    CalculateAll(rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}



void BeamElement::CalculateSectionProperties()
{
    KRATOS_TRY

    array_1d<double, 3> x_0;
    array_1d<double, 3> x_1; // Vector que contiene coordenadas de los nodos.
    array_1d<double, 3> length;   // Vector que contiene la direccion de la barra.

//        double minimo, maximo, B;
//        const double b        = GetProperties()[BASE];
//        const double h        = GetProperties()[HEIGHT];

    if( GetProperties().Has(CROSS_AREA) )
        mArea = GetProperties()[CROSS_AREA];
    else
        mArea = GetValue(AREA);

    Matrix* inertia;
    if( GetProperties().Has(LOCAL_INERTIA) )
    {
        inertia = &(GetProperties()[LOCAL_INERTIA]);
    }
    else if( GetProperties().Has(INERTIA) )
    {
        inertia = &(GetProperties()[INERTIA]);
    }
    else if( Has(LOCAL_INERTIA) )
    {
        inertia = &(GetValue(LOCAL_INERTIA));
    }
    else if( Has(INERTIA) )
    {
        inertia = &(GetValue(INERTIA));
    }
    else
        KRATOS_ERROR << "The Inertia is not fully defined for the element";
    mInertia_x = (*inertia)(0,0);
    mInertia_y = (*inertia)(1,1);
    mInertia_Polar = (*inertia)(0,1);

//        mInertia_x     = b * h * h * h / 12.0;
//        mInertia_y     = b * b * b * h / 12.0;
//        minimo         = std::min( b, h );
//        maximo        = std::max( b, h );
//        B        = ( 1.00 - 0.63 * ( minimo / maximo ) * ( 1 - ( pow( minimo, 4 ) / ( 12 * pow( maximo, 4 ) ) ) ) ) / 3; // constante torsional. Solo para secciones rectangulares.
//        mInertia_Polar = B * minimo * minimo * minimo * maximo;
//        mArea        = b * h;


    x_0( 0 ) = GetGeometry()[0].X0() + mInitialDisp(0, 0);
    x_0( 1 ) = GetGeometry()[0].Y0() + mInitialDisp(0, 1);
    x_0( 2 ) = GetGeometry()[0].Z0() + mInitialDisp(0, 2);
    x_1( 0 ) = GetGeometry()[1].X0() + mInitialDisp(1, 0);
    x_1( 1 ) = GetGeometry()[1].Y0() + mInitialDisp(1, 1);
    x_1( 2 ) = GetGeometry()[1].Z0() + mInitialDisp(1, 2);

    noalias( length ) = x_1 - x_0;
    mlength = std::sqrt( inner_prod( length, length ) );

    if (mlength == 0.00)
        KRATOS_ERROR << "Zero length found in elemnet #" << this->Id();

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void BeamElement::CalculateLocalMatrix(Matrix& LocalMatrix)
{
    KRATOS_TRY

    if(LocalMatrix.size1()!=12 || LocalMatrix.size2()!=12)   // Matriz local de rigidez de la Estructura.
        LocalMatrix.resize(12,12,false);
    noalias(LocalMatrix)   = zero_matrix<double>(12,12);

    //Inicializando matriz local de la estructura

    //const double mlength = GetGeometry().Length();
    const double Poisson = GetProperties()[POISSON_RATIO];
    const double Youngs  = GetProperties()[YOUNG_MODULUS];
    const double Elasticidad_Cortante   = Youngs /(2.0*(1.0 + Poisson));

    const double L  =    mlength;
    const double LL   =  mlength* mlength;
    const double LLL    =    mlength* mlength * mlength;

    double const EA   =  mArea          * Youngs;
    double const EIx  =  mInertia_x     * Youngs;
    double const EIy  =  mInertia_y     * Youngs;
    double const JG   =  mInertia_Polar * Elasticidad_Cortante;

    LocalMatrix(0,0)    =   (EA)/(L);
    LocalMatrix(6,0)    =   -(EA)/(L);

    LocalMatrix(1,1)    =   (12*EIx)/(LLL);
    LocalMatrix(5,1)    =   (6*EIx)/(LL);
    LocalMatrix(7,1)    =   -(12*EIx)/(LLL);
    LocalMatrix(11,1)   =   (6*EIx)/(LL);

    LocalMatrix(2,2)    =   (12*EIy)/(LLL);
    LocalMatrix(4,2)    =   -(6*EIy)/(LL);
    LocalMatrix(8,2)    =   -(12*EIy)/(LLL);
    LocalMatrix(10,2)   =   -(6*EIy)/(LL);

    LocalMatrix(3,3)    =   (JG)/L;
    LocalMatrix(9,3)    =   -(JG)/L;

    LocalMatrix(2,4)    =   -(6*EIy)/(LL);
    LocalMatrix(4,4)    =   (4*EIy)/L;
    LocalMatrix(8,4)    =    (6*EIy)/(LL);
    LocalMatrix(10,4)   =    (2*EIy)/L;

    LocalMatrix(1,5)    =   (6*EIx)/(LL);
    LocalMatrix(5,5)    =   (4*EIx)/L;
    LocalMatrix(7,5)    =   -(6*EIx)/(LL);
    LocalMatrix(11,5)   =   (2*EIx)/L;

    LocalMatrix(0,6)    =    -(EA)/( L);
    LocalMatrix(6,6)    =   (EA)/( L);

    LocalMatrix(1,7)    =   -(12*EIx)/(LLL);
    LocalMatrix(5,7)    =   -(6*EIx)/(LL);
    LocalMatrix(7,7)    =    (12*EIx)/(LLL);
    LocalMatrix(11,7)   =   -(6*EIx)/(LL);

    LocalMatrix(2,8)    =    -(12*EIy)/(LLL);
    LocalMatrix(4,8)    =    (6*EIy)/(LL);
    LocalMatrix(8,8)    =    (12*EIy)/(LLL);
    LocalMatrix(10,8)   =    (6*EIy)/(LL);

    LocalMatrix(3,9)    =   -(JG)/L;
    LocalMatrix(9,9)    =    (JG)/L;

    LocalMatrix(2,10)   =   -(6*EIy)/(LL);
    LocalMatrix(4,10)   =  (2*EIy)/L;
    LocalMatrix(8,10)   =  (6*EIy)/(LL);
    LocalMatrix(10,10)=  (4*EIy)/L;

    LocalMatrix(1,11)   =   (6*EIx)/(LL);
    LocalMatrix(5,11)   =   (2*EIx)/L;
    LocalMatrix(7,11)   =  -(6*EIx)/(LL);
    LocalMatrix(11,11)=  (4*EIx)/L;

    KRATOS_CATCH("")
}

//*****************************************************************************
//*****************************************************************************
void BeamElement::CalculateTransformationMatrix(Matrix& Rotation)
{
    KRATOS_TRY

    Vector Normal_zero(9); // vector que contiene los cosenos directores.
    Vector x_zero(6);
    Vector Vector_zero(3);
    noalias(Normal_zero) =  zero_vector<double>(9);
    noalias(x_zero)      =  zero_vector<double>(6);
    noalias(Vector_zero) =  zero_vector<double>(3);
    noalias(Rotation)    =  zero_matrix<double>(12, 12);

    double nx, ny, nz,teta/*, phi*/;

    x_zero(0) = GetGeometry()[0].X0();
    x_zero(1) = GetGeometry()[0].Y0();
    x_zero(2) = GetGeometry()[0].Z0();
    x_zero(3) = GetGeometry()[1].X0();
    x_zero(4) = GetGeometry()[1].Y0();
    x_zero(5) = GetGeometry()[1].Z0();

    for (unsigned int i=0; i<3; i++)
    {
        Vector_zero[i] = x_zero[i+3] - x_zero[i];
    }

    double length_inverse = ( 1.00 / mlength );
    for ( unsigned int i = 0; i < 3; i++ )
    {
        Normal_zero[i] = Vector_zero[i] * length_inverse;
    }

    nx = Normal_zero[0];
    ny = Normal_zero[1];
    nz = Normal_zero[2];

    if (nx ==0.0)
    {
        teta = SD_MathUtils<double>::Pi()/2;
        if (ny == 0.0)
        {
            teta = 0.0;
//            phi  = SD_MathUtils<double>::Pi()/2.0;
        }
//        else
//        {
//            phi = atan(nz/sqrt(nx*nx+ny*ny));
//        }
    }
    else
    {
        teta = atan(ny/nx);
//        phi  = atan(nz/sqrt(nx*nx+ny*ny));
    }

    if(nx < 0.0)
        teta   = teta + SD_MathUtils<double>::Pi();

    Normal_zero[3] = -sin(teta);
    Normal_zero[4] =  cos(teta);
    Normal_zero[5] =  0.0;
    Normal_zero[6] = -nz*cos(teta);
    Normal_zero[7] = -nz*sin(teta);
    Normal_zero[8] =  nx*cos(teta) + ny*sin(teta);

    // Creacion de la matriz de transformacion.
    for (unsigned int kk=0; kk < 12; kk += 3)
    {
        for (unsigned int i=0; i<3; i++)
        {
            for(unsigned int j=0; j<3; j++)
            {
                Rotation(i+kk,j+kk)=Normal_zero(3*j+i);
            }
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void BeamElement::CalculateBodyForce(const Matrix& Rotation, Vector& LocalBody, Vector& GlobalBody)
{
    KRATOS_TRY
    //Creacion de los vectores de cargas externas.
    // Fuerzas externas Uniformente distriduida. Por lo general es un dato suministrado por el usuario.
    // Cambiaro de posicion una vez terminado el programa.

    double alpha  =  0.00;
    double signo  =  1.00;
    //const double mlength = GetGeometry().Length();
    double  sino;
    double  cose;

    array_1d<double, 3> Weight;
    Weight[0]        =  GetProperties()[GRAVITY](0) * GetProperties()[DENSITY] + GetProperties()[BODY_FORCE](0);
    Weight[1]        =  GetProperties()[GRAVITY](1) * GetProperties()[DENSITY] + GetProperties()[BODY_FORCE](1);
    Weight[2]        =  GetProperties()[GRAVITY](2) * GetProperties()[DENSITY] + GetProperties()[BODY_FORCE](2);


    array_1d<double, 12 > Cargas_X = ZeroVector(12);
    array_1d<double, 12 > Cargas_Y = ZeroVector(12);
    array_1d<double, 12 > Cargas_Z = ZeroVector(12);

    array_1d<double, 2 > Load;
    array_1d<double, 6 > x_zero;

    Vector Normal_Loads;
    Vector Vector_zero;

    Normal_Loads.resize(3,false);
    Vector_zero.resize(3,false);

    x_zero(0)= GetGeometry()[0].X0();
    x_zero(1)= GetGeometry()[0].Y0();
    x_zero(2)= GetGeometry()[0].Z0();
    x_zero(3)= GetGeometry()[1].X0();
    x_zero(4)= GetGeometry()[1].Y0();
    x_zero(5)= GetGeometry()[1].Z0();

    for (unsigned int i=0; i<3; i++)
    {
        Vector_zero[i] = x_zero[i+3] - x_zero[i];
    }

    //Fuerza En X
    //***********************************
    if(Weight[0]!=0.00)
    {
        Normal_Loads[0]   = 0.00;
        Normal_Loads[1]   = Vector_zero[1] ;
        Normal_Loads[2]   = Vector_zero[2] ;

        if (Vector_zero[0]<0)
        {
            signo =-1.00;
        }
        if( norm_2(Normal_Loads)==0 || norm_2( Vector_zero)==0  )
        {
            alpha = signo*SD_MathUtils<double>::Pi()/2;
        }
        else
        {
            alpha = inner_prod(Normal_Loads,Vector_zero)/(norm_2(Vector_zero)*norm_2( Normal_Loads));
            alpha   = signo*acos(alpha);
        }

        sino   = sin(alpha);
        cose   = cos(alpha);;
        if(fabs(sino) < 1E-7) sino = 0.00;
        if(fabs(cose) < 1E-7) cose = 0.00;

        // las fuerzas consideradas son las de peso propio.
        Load[0]= mArea*Weight[0]*sino;         // Carga Axialmente Distribuida.
        Load[1]= mArea*Weight[0]*cose;         // Carga en la Direccion gravedad

        Cargas_X[0]=   Load[0]*mlength/2.00;     // Fuerza en X;
        Cargas_X[1]=   -(Load[1]*mlength)/2.00;  // Fuerza en Y; graveded
        Cargas_X[2]=   0.00;                                                                                                     // Fuerza en Z
        Cargas_X[3]=   0.00;                                                                                                     // Momento Tersor X;
        Cargas_X[4]=   0.00;                                                                                                     // Momento Y
        Cargas_X[5]=  -(Load[1])*mlength*mlength/12.00;;                                        // Momento Z
        Cargas_X[6]=   Load[0]*mlength/2.00;
        Cargas_X[7]=   -(Load[1])*mlength/2.00;
        Cargas_X[8]=    0.00;
        Cargas_X[9]=    0.00;
        Cargas_X[10]=   0.00;
        Cargas_X[11]=   (Load[1])*mlength*mlength/12.00;

        noalias(GlobalBody) = prod(Rotation,Cargas_X);      // Cargas externas en coordenadas globales.
        noalias(LocalBody)  = Cargas_X;
    }

    //Fuerza En Z
    //***********************************
    if(Weight[2]!=0.00)
    {
        Normal_Loads[0] = Vector_zero[0] ;
        Normal_Loads[1] = Vector_zero[1] ;
        Normal_Loads[2] = 0.00;

        if (Vector_zero[2]<0)
        {
            signo =-1.00;
        }
        if( norm_2(Normal_Loads)==0 || norm_2( Vector_zero)==0  )
        {
            alpha = signo*SD_MathUtils<double>::Pi()/2;
        }
        else
        {
            alpha = inner_prod(Normal_Loads,Vector_zero)/(norm_2(Vector_zero)*norm_2( Normal_Loads));
            alpha   = signo*acos(alpha);
        }

        sino = sin(alpha);
        cose = cos(alpha);

        if(fabs(sino) < 1E-7) sino = 0.00;
        if(fabs(cose) < 1E-7) cose = 0.00;

        // las fuerzas consideradas son las de peso propio.
        Load[0]= mArea*Weight[2]*sino;         // Carga Axialmente Distribuida.
        Load[1]= mArea*Weight[2]*cose;         // Carga en la Direccion gravedad

        Cargas_Z[0]=   -Load[0]*mlength/2.00;     // Fuerza en X;
        Cargas_Z[1]=   0.00;
        Cargas_Z[2]=   -(Load[1]*mlength)/2.00;  // Fuerza en Z; graveded                                                                                                    // Fuerza en Z
        Cargas_Z[3]=   0.00;                                                                                                     // Momento Tersor X;
        Cargas_Z[4]=   -(Load[1])*mlength*mlength/12.00;                                                                                                 // Momento Y
        Cargas_Z[5]=    0.00;
        Cargas_Z[6]=   -Load[0]*mlength/2.00;
        Cargas_Z[7]=   0.00;
        Cargas_Z[8]=    -(Load[1])*mlength/2.00;
        Cargas_Z[9]=   0.00;
        Cargas_Z[10]= (Load[1])*mlength*mlength/12.00;
        Cargas_Z[11]=  0.00;

        noalias(GlobalBody) = prod(Rotation,Cargas_Z);      // Cargas externas en coordenadas globales.
        noalias(LocalBody)  = Cargas_Z;
    }

    //Fuerza En Y
    //***********************************
    if(Weight[1]!=0.00)
    {
        Normal_Loads      = ZeroVector(3);
        Normal_Loads[0] = Vector_zero[0] ;
        Normal_Loads[1] = 0.00 ;
        Normal_Loads[2] = Vector_zero[2];

        if (Vector_zero[1]<0)
        {
            signo =-1.00;
        }
        if( norm_2(Normal_Loads)==0 || norm_2( Vector_zero)==0  )
        {
            alpha = signo*SD_MathUtils<double>::Pi()/2;
        }
        else
        {
            alpha = inner_prod(Normal_Loads,Vector_zero)/(norm_2(Vector_zero)*norm_2( Normal_Loads));
            alpha   = signo*acos(alpha);
        }

        sino = sin(alpha);
        cose = cos(alpha);

        if(fabs(sino) < 1E-7) sino = 0.00;
        if(fabs(cose) < 1E-7) cose = 0.00;

        // las fuerzas consideradas son las de peso propio.
        Load[0]= mArea*Weight[1]*sino;         // Carga Axialmente Distribuida.
        Load[1]= mArea*Weight[1]*cose;         // Carga en la Direccion gravedad


        Cargas_Y[0]=   -Load[0]*mlength/2.00;     // Fuerza en X;
        Cargas_Y[1]=   -(Load[1]*mlength)/2.00;  // Fuerza en Y; graveded
        Cargas_Y[2]=   0.00;                                                                                                     // Fuerza en Z
        Cargas_Y[3]=   0.00;                                                                                                     // Momento Tersor X;
        Cargas_Y[4]=   0.00;                                                                                                     // Momento Y
        Cargas_Y[5]=  -(Load[1])*mlength*mlength/12.00;;                                        // Momento Z
        Cargas_Y[6]=   -Load[0]*mlength/2.00;
        Cargas_Y[7]=   -(Load[1])*mlength/2.00;
        Cargas_Y[8]=    0.00;
        Cargas_Y[9]=    0.00;
        Cargas_Y[10]=   0.00;
        Cargas_Y[11]=   (Load[1])*mlength*mlength/12.00;

        noalias(GlobalBody) = prod(Rotation,Cargas_Y);      // Cargas externas en coordenadas globales.
        noalias(LocalBody)  = Cargas_Y;
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void BeamElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int NumberOfNodes = GetGeometry().size();
    unsigned int MatSize = dimension * NumberOfNodes;
    if(rMassMatrix.size1() != MatSize)
        rMassMatrix.resize(MatSize,MatSize,false);

    rMassMatrix = ZeroMatrix(MatSize,MatSize);
    //const double& mlength = GetGeometry().Length();

    double TotalMass = mArea*mlength*GetProperties()[DENSITY];

    Vector LumpFact;
    LumpFact = GetGeometry().LumpingFactors(LumpFact);

    for(unsigned int i=0; i<NumberOfNodes; i++)
    {
        double temp = LumpFact[i]*TotalMass;
        for(unsigned int j=0; j<dimension; j++)
        {
            unsigned int index = i*dimension + j;
            rMassMatrix(index,index) = temp;
            if (index==3 || index==4 || index==5)
                rMassMatrix(index,index) = 0.00;
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void BeamElement::CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int mat_size = number_of_nodes * dim;

    if ( rDampMatrix.size1() != mat_size )
        rDampMatrix.resize( mat_size, mat_size, false );

    noalias( rDampMatrix ) = ZeroMatrix( mat_size, mat_size );

    Matrix StiffnessMatrix = ZeroMatrix( mat_size, mat_size );

    Vector RHS_Vector = ZeroVector( mat_size );

    //rayleigh damping
    CalculateAll( StiffnessMatrix, RHS_Vector, rCurrentProcessInfo, true, false );

    double alpha = 0.0, beta = 0.0;

    if(GetProperties().Has(RAYLEIGH_DAMPING_ALPHA))
    {
        alpha = GetProperties()[RAYLEIGH_DAMPING_ALPHA];
    }

    if(GetProperties().Has(RAYLEIGH_DAMPING_BETA))
    {
        beta = GetProperties()[RAYLEIGH_DAMPING_BETA];
    }

    CalculateMassMatrix( rDampMatrix, rCurrentProcessInfo );

    noalias( rDampMatrix ) = alpha * rDampMatrix + beta * StiffnessMatrix;

    // KRATOS_WATCH(rDampMatrix)

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void BeamElement::GetValuesVector(Vector& values, int Step) const
{
    if(values.size() != GetGeometry().size()*6)
        values.resize(GetGeometry().size()*6, false);

    for (std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        values(6*i)   = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X, Step);
        values(6*i+1) = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y, Step);
        values(6*i+2) = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z, Step);
        values(6*i+3) = GetGeometry()[i].GetSolutionStepValue(ROTATION_X, Step);
        values(6*i+4) = GetGeometry()[i].GetSolutionStepValue(ROTATION_Y, Step);
        values(6*i+5) = GetGeometry()[i].GetSolutionStepValue(ROTATION_Z, Step);
    }
}

//************************************************************************************
//************************************************************************************
void BeamElement::GetFirstDerivativesVector(Vector& values, int Step) const
{
    if(values.size() != GetGeometry().size()*6)
        values.resize(GetGeometry().size()*6, false);

    for (std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        values(6*i)   = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_DT_X, Step);
        values(6*i+1) = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_DT_Y, Step);
        values(6*i+2) = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_DT_Z, Step);
        values(6*i+3) = GetGeometry()[i].GetSolutionStepValue(ROTATION_DT_X, Step);
        values(6*i+4) = GetGeometry()[i].GetSolutionStepValue(ROTATION_DT_Y, Step);
        values(6*i+5) = GetGeometry()[i].GetSolutionStepValue(ROTATION_DT_Z, Step);
    }
}

//************************************************************************************
//************************************************************************************
void BeamElement::GetSecondDerivativesVector(Vector& values, int Step) const
{
    if(values.size() != GetGeometry().size()*6)
        values.resize(GetGeometry().size()*6, false);

    for (std::size_t i = 0; i < GetGeometry().size(); ++i)
    {
        values(6*i)   = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X, Step);
        values(6*i+1) = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y, Step);
        values(6*i+2) = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z, Step);
        values(6*i+3) = GetGeometry()[i].GetSolutionStepValue(ANGULAR_ACCELERATION_X, Step);
        values(6*i+4) = GetGeometry()[i].GetSolutionStepValue(ANGULAR_ACCELERATION_Y, Step);
        values(6*i+5) = GetGeometry()[i].GetSolutionStepValue(ANGULAR_ACCELERATION_Z, Step);
    }
}

//************************************************************************************
//************************************************************************************
void BeamElement::SetValuesOnIntegrationPoints( const Variable<Vector>& rVariable,
        const std::vector<Vector>& Output,
        const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == PRESTRESS || rVariable == INSITU_STRESS)
    {
        // KRATOS_WATCH(Id())
        // KRATOS_WATCH(Output.size())
        // if (Output.size() > 0)
        //     KRATOS_WATCH(Output[0])
        // KRATOS_WATCH(mPreForces)
        noalias(mPreForces) = Output[0];
        // std::cout << Id() << ": Preforces is set to " << mPreForces << std::endl;
    }
}

//************************************************************************************
//************************************************************************************
void BeamElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& Output,
        const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( GetIntegrationMethod() );

    if(Output.size() != integration_points.size())
        Output.resize(integration_points.size());

    for (std::size_t i = 0; i < integration_points.size(); ++i)
    {
        Output[i] = 0.0;
    }
}

//************************************************************************************
//************************************************************************************
void BeamElement::CalculateOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
        std::vector< array_1d<double,3> >& Output,
        const ProcessInfo& rCurrentProcessInfo)
{
    Vector Stress;
    Vector Load1;
    Vector Load2;
    Vector Load3;
    CalculateLocalNodalStress(Stress);

    for(unsigned int i = 0; i<Stress.size(); i++)
    {
        if( std::fabs(Stress[i])< 1E-6) Stress[i] = 0.00;
    }

    double factor = 1.0;

    // int factor = 1;
    // double x_toler     = GetGeometry()[1].X0() - GetGeometry()[0].X0();
    // double y_toler     = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
    // const double toler = 1E-6;

    // if(fabs(x_toler)>toler)
    // {
    //     if(GetGeometry()[1].X0() > GetGeometry()[0].X0())
    //         factor = 1;
    // }

    // else if(fabs(y_toler)>toler)
    // {
    //     if(GetGeometry()[1].Y0() > GetGeometry()[0].Y0())
    //         factor = 1;
    // }
    // else
    //     factor = 1; //-1;

    CalculateDistributedBodyForce(1, Load1);
    CalculateDistributedBodyForce(2, Load2);
    CalculateDistributedBodyForce(3, Load3);

    const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( GetIntegrationMethod() );

    if(Output.size() != integration_points.size())
        Output.resize(integration_points.size());

    if(rVariable==MOMENT)
    {
        /***** old *********/
//         /// Punto Inical
//         Output[0][0] = 0.00;  //Stress[3];
//         Output[0][1] = 0.00;  //Stress[4];
//         Output[0][2] = factor * CalculateInternalMoment(Stress[5], Stress[11], Load2[1], 1.00/4.00); //Stress[5];
// //           Output[0][2] = factor * CalculateInternalMoment(Stress[5], Stress[1], Load2[1], mlength/4.00); //Stress[5];
//         // hbui: It is noted that, the location of the integration point is not the one from Gauss quadrature

//         Output[1][0] = 0.00;
//         Output[1][1] = 0.00;
//         Output[1][2] = factor * CalculateInternalMoment(Stress[5], Stress[11], Load2[1], 1./2);
// //           Output[1][2] = factor * CalculateInternalMoment(Stress[5], Stress[1], Load2[1], mlength/2);


//         Output[2][0] = 0.00;
//         Output[2][1] = 0.00;
//         Output[2][2] = factor * CalculateInternalMoment(Stress[5], Stress[11], Load2[1], 3.00/4.00);
// //           Output[2][2] = factor * CalculateInternalMoment(Stress[5], Stress[1], Load2[1], 3.00*mlength/4.00);
        /*******end of old *********/

        /***** new *********/
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            double xi = (integration_points[point].X()+1.0) / 2;

            Output[point][0] = 0.00;  //Stress[3];
            Output[point][1] = 0.00;  //Stress[4];
            Output[point][2] = factor * CalculateInternalMoment(Stress[5], Stress[11], Load2[1], xi); //Stress[5];
        }
        /*******end of new *********/
    }
    else if(rVariable==FORCE)
    {
        /***** old *********/
        // Output[0][0] = factor * CalculateInternalAxil(Stress[0], Load2[0], mlength/4.00);
        // Output[0][1] = factor * CalculateInternalShear(Stress[1], Load2[1], mlength/4.00);
        // Output[0][2] = 0.00;

        // Output[1][0] = factor * CalculateInternalAxil(Stress[0], Load2[0], mlength/2.00);
        // Output[1][1] = factor * CalculateInternalShear(Stress[1], Load2[1], mlength/2.00);
        // Output[1][2] = 0.00;

        // Output[2][0] = factor * CalculateInternalAxil(Stress[0], Load2[0],  3.00 * mlength/4.00);
        // Output[2][1] = factor * CalculateInternalShear(Stress[1], Load2[1], 3.00 * mlength/4.00);
        // Output[2][2] = 0.00;
        /*******end of old *********/

        /***** new *********/
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            double xi = (integration_points[point].X()+1.0) / 2;

            Output[point][0] = factor * CalculateInternalAxil(Stress[0], Load2[0], mlength*xi);
            Output[point][1] = factor * CalculateInternalShear(Stress[1], Load2[1], mlength*xi);
            Output[point][2] = 0.00;
        }
        /*******end of new *********/
    }
    else
    {
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            noalias(Output[point]) = ZeroVector(3);
        }
    }
}

//************************************************************************************
//************************************************************************************
void BeamElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& Output,
        const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == STRESSES)
    {
        if (Output.size() != 1)
            Output.resize(1);

        // Matrix GlobalMatrix;
        // CalculateLHS(GlobalMatrix);

        // Vector CurrentDisplacement(12);
        // CurrentDisplacement(0)      =   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X) - mInitialDisp(0, 0);
        // CurrentDisplacement(1)      =   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y) - mInitialDisp(0, 1);
        // CurrentDisplacement(2)      =   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z) - mInitialDisp(0, 2);
        // CurrentDisplacement(3)      =   GetGeometry()[0].GetSolutionStepValue(ROTATION_X) - mInitialRot(0, 0);
        // CurrentDisplacement(4)      =   GetGeometry()[0].GetSolutionStepValue(ROTATION_Y) - mInitialRot(0, 1);
        // CurrentDisplacement(5)      =   GetGeometry()[0].GetSolutionStepValue(ROTATION_Z) - mInitialRot(0, 2);
        // CurrentDisplacement(6)      =   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X) - mInitialDisp(1, 0);
        // CurrentDisplacement(7)      =   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y) - mInitialDisp(1, 1);
        // CurrentDisplacement(8)      =   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z) - mInitialDisp(1, 2);
        // CurrentDisplacement(9)      =   GetGeometry()[1].GetSolutionStepValue(ROTATION_X) - mInitialRot(1, 0);
        // CurrentDisplacement(10)     =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Y) - mInitialRot(1, 1);
        // CurrentDisplacement(11)     =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Z) - mInitialRot(1, 2);

        // Output[0].resize(12, false);
        // noalias(Output[0]) = prod(GlobalMatrix, CurrentDisplacement);

        Output[0].resize(12, false);
        noalias(Output[0]) = mCurrentForces;
    }
    else if (rVariable == INSITU_STRESS || rVariable == PRESTRESS )
    {
        if (Output.size() != 1)
            Output.resize(1);

        Output[0].resize(12, false);
        noalias(Output[0]) = mPreForces;
    }
}

double BeamElement::CalculateInternalMoment(const double& Mo, const double& Vo, const double& Load, const double& X)
{
//       return Mo - Vo*X + 0.5 * Load * X * X;
    return Mo *(1-X) - Vo*X;
}

double BeamElement::CalculateInternalShear(const double& Vo, const double& Load, const double& X)
{
    return  -Vo + Load * X;
}

double BeamElement::CalculateInternalAxil(const double& Ao, const double& Load, const double& X)
{
    return  -Ao +  Load * X;
}

IntegrationMethod  BeamElement::GetIntegrationMethod() const
{
    return IntegrationMethod::GI_GAUSS_3;
}

void BeamElement::CalculateLocalNodalStress(Vector& Stress)
{
    Matrix Rotation;
    Matrix LocalMatrix;
    array_1d<double, 12 > CurrentDisplacement;
    array_1d<double, 12 > LocalDisplacement;
    Vector LocalBody  = ZeroVector(12);
    Vector GlobalBody = ZeroVector(12);
    Rotation.resize(12,12, false);
    Stress.resize(12, false);

    CurrentDisplacement(0)      =   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X);
    CurrentDisplacement(1)      =   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y);
    CurrentDisplacement(2)      =   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z);
    CurrentDisplacement(3)      =   GetGeometry()[0].GetSolutionStepValue(ROTATION_X);
    CurrentDisplacement(4)      =   GetGeometry()[0].GetSolutionStepValue(ROTATION_Y);
    CurrentDisplacement(5)      =   GetGeometry()[0].GetSolutionStepValue(ROTATION_Z);
    CurrentDisplacement(6)      =   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X);
    CurrentDisplacement(7)      =   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y);
    CurrentDisplacement(8)      =   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z);
    CurrentDisplacement(9)      =   GetGeometry()[1].GetSolutionStepValue(ROTATION_X);
    CurrentDisplacement(10)         =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Y);
    CurrentDisplacement(11)         =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Z);

    CalculateTransformationMatrix(Rotation);
    CalculateLocalMatrix(LocalMatrix);
    noalias(LocalDisplacement) = prod(Matrix(trans(Rotation)), CurrentDisplacement);
    CalculateBodyForce(Rotation, LocalBody, GlobalBody);
    noalias(Stress) = -LocalBody + prod(LocalMatrix, LocalDisplacement);
//      noalias(Stress) = -LocalBody + prod(Matrix(prod(Rotation,LocalMatrix)), LocalDisplacement);

    return;
}

void BeamElement::CalculateDistributedBodyForce(const int Direction, Vector& Load)
{
    array_1d<double, 3> Weight;
    Load.resize(2, false);
    Weight[0]        =  GetProperties()[GRAVITY](0) * GetProperties()[DENSITY] + GetProperties()[BODY_FORCE](0);
    Weight[1]        =  GetProperties()[GRAVITY](1) * GetProperties()[DENSITY] + GetProperties()[BODY_FORCE](1);
    Weight[2]        =  GetProperties()[GRAVITY](2) * GetProperties()[DENSITY] + GetProperties()[BODY_FORCE](2);

    double alpha  =  0.00;
    double signo  =  1.00;
    double  sino;
    double  cose;

    array_1d<double, 6 > x_zero;

    Vector Normal_Loads;
    Vector Vector_zero;
    Normal_Loads.resize(3,false);
    Vector_zero.resize(3,false);

    x_zero(0)= GetGeometry()[0].X0();
    x_zero(1)= GetGeometry()[0].Y0();
    x_zero(2)= GetGeometry()[0].Z0();
    x_zero(3)= GetGeometry()[1].X0();
    x_zero(4)= GetGeometry()[1].Y0();
    x_zero(5)= GetGeometry()[1].Z0();

    for (unsigned int i=0; i<3; i++)
    {
        Vector_zero[i] = x_zero[i+3] - x_zero[i];
    }

    if(Direction==1)
        Normal_Loads[0]   = 0.00;
    Normal_Loads[1]   = Vector_zero[1] ;
    Normal_Loads[2]   = Vector_zero[2] ;

    if (Vector_zero[0]<0)
    {
        signo =-1.00;
    }
    if( norm_2(Normal_Loads)==0 || norm_2( Vector_zero)==0  )
    {
        alpha = signo*SD_MathUtils<double>::Pi()/2;
    }
    else
    {
        alpha = inner_prod(Normal_Loads,Vector_zero)/(norm_2(Vector_zero)*norm_2( Normal_Loads));
        alpha   = signo*acos(alpha);
    }

    sino = sin(alpha);
    cose = cos(alpha);

    if(fabs(sino) < 1E-7) sino = 0.00;
    if(fabs(cose) < 1E-7) cose = 0.00;

    // las fuerzas consideradas son las de peso propio.
    Load[0]= mArea*Weight[0]*sino;         // Carga Axialmente Distribuida.
    Load[1]= mArea*Weight[0]*cose;         // Carga en la Direccion gravedad


    if(Direction==2) // 1=x, 2=y, 3=z
    {
        Normal_Loads    = ZeroVector(3);
        Normal_Loads[0] = Vector_zero[0] ;
        Normal_Loads[1] = 0.00 ;
        Normal_Loads[2] = Vector_zero[2];

        if (Vector_zero[1]<0)
        {
            signo =-1.00;
        }
        if( norm_2(Normal_Loads)==0 || norm_2( Vector_zero)==0  )
        {
            alpha = signo*SD_MathUtils<double>::Pi()/2;
        }
        else
        {
            alpha = inner_prod(Normal_Loads,Vector_zero)/(norm_2(Vector_zero)*norm_2( Normal_Loads));
            alpha = signo*acos(alpha);
        }

        sino = sin(alpha);
        cose = cos(alpha);

        if(fabs(sino) < 1E-7) sino = 0.00;
        if(fabs(cose) < 1E-7) cose = 0.00;


        // las fuerzas consideradas son las de peso propio.
        Load[0]= mArea*Weight[1]*sino;         // Carga Axialmente Distribuida.
        Load[1]= mArea*Weight[1]*cose;         // Carga en la Direccion gravedad
    }

    if(Direction==3) // 1=x, 2=y, 3=z
    {
        Normal_Loads[0] = Vector_zero[0] ;
        Normal_Loads[1] = Vector_zero[1] ;
        Normal_Loads[2] = 0.00;

        if (Vector_zero[2]<0)
        {
            signo =-1.00;
        }
        if( norm_2(Normal_Loads)==0 || norm_2( Vector_zero)==0  )
        {
            alpha = signo*SD_MathUtils<double>::Pi()/2;
        }
        else
        {
            alpha = inner_prod(Normal_Loads,Vector_zero)/(norm_2(Vector_zero)*norm_2( Normal_Loads));
            alpha   = signo*acos(alpha);
        }

        sino = sin(alpha);
        cose = cos(alpha);

        if(fabs(sino) < 1E-7) sino = 0.00;
        if(fabs(cose) < 1E-7) cose = 0.00;

        // las fuerzas consideradas son las de peso propio.
        Load[0]= mArea*Weight[2]*sino;         // Carga Axialmente Distribuida.
        Load[1]= mArea*Weight[2]*cose;         // Carga en la Direccion gravedad
    }
}

/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  BeamElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    //verify that the variables are correctly initialized
    if(VELOCITY.Key() == 0)
        KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered";
    if(DISPLACEMENT.Key() == 0)
        KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered";
    if(ACCELERATION.Key() == 0)
        KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered";
    if(DENSITY.Key() == 0)
        KRATOS_ERROR << "DENSITY has Key zero! (check if the application is correctly registered";
    if(GRAVITY.Key() == 0)
        KRATOS_ERROR << "GRAVITY has Key zero! (check if the application is correctly registered";
    if(BODY_FORCE.Key() == 0)
        KRATOS_ERROR << "BODY_FORCE has Key zero! (check if the application is correctly registered";
    if(CROSS_AREA.Key() == 0)
        KRATOS_ERROR << "CROSS_AREA has Key zero! (check if the application is correctly registered";
    if(LOCAL_INERTIA.Key() == 0)
        KRATOS_ERROR << "LOCAL_INERTIA has Key zero! (check if the application is correctly registered";
    if(ROTATION.Key() == 0)
        KRATOS_ERROR << "ROTATION has Key zero! (check if the application is correctly registered";

    //verify that the dofs exist
    for(unsigned int i=0; i<this->GetGeometry().size(); i++)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
            KRATOS_ERROR << "missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id();
        if(this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false)
            KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id();
    }

    //verify that the area and inertia is given by properties
    if( GetProperties().Has(CROSS_AREA) || Has(AREA) )
    {}
    else
        KRATOS_THROW_ERROR(std::logic_error, "The Area is not fully defined for the beam element", Id())

    const Matrix* inertia;
    if( GetProperties().Has(LOCAL_INERTIA) )
    {
        inertia = &(GetProperties()[LOCAL_INERTIA]);
    }
    else if( GetProperties().Has(INERTIA) )
    {
        inertia = &(GetProperties()[INERTIA]);
    }
    else if( Has(LOCAL_INERTIA) )
    {
        inertia = &(GetValue(LOCAL_INERTIA));
    }
    else if( Has(INERTIA) )
    {
        inertia = &(GetValue(INERTIA));
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "The Inertia is not fully defined for the beam element", Id())

    if(inertia->size1() < 2 || inertia->size2() < 2)
    {
        std::stringstream ss;
        ss << "Error on BeamElement " << Id() << ", the inertia is not of expected size, true size = ("
           << inertia->size1() << ", " << inertia->size2() << ")" << std::endl;
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    //Verify that the body force is defined
    if (this->GetProperties().Has(BODY_FORCE)==false)
    {
        KRATOS_THROW_ERROR(std::logic_error,"BODY_FORCE not provided for property ",this->GetProperties().Id())
    }

    return 0;

    KRATOS_CATCH("");
}

} // Namespace Kratos
