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
giang.bui@rub.de
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
//   Modified by:         $Author: virginia $
//   Date:                $Date: 2009-01-23 14:39:59 $
//   Revision:            $Revision: 1.27 $
//
//


// System includes

// External includes


// Project includes
#include "utilities/math_utils.h"
#include "custom_elements/total_lagrangian.h"
#include "custom_utilities/geometry_utility.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/eigen_utility.h"
#include "structural_application_variables.h"

// #define CHECK_DEFORMATION_GRADIENT

namespace Kratos
{

    TotalLagrangian::TotalLagrangian( IndexType NewId, GeometryType::Pointer pGeometry )
            : Element( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
        mIsInitialized = false;
    }

//************************************************************************************
//************************************************************************************

    TotalLagrangian::TotalLagrangian( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : Element( NewId, pGeometry, pProperties )
    {
        //         const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

        mIsInitialized = false;
    }

    Element::Pointer TotalLagrangian::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new TotalLagrangian( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    Element::Pointer TotalLagrangian::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new TotalLagrangian( NewId, pGeom, pProperties ) );
    }

    TotalLagrangian::~TotalLagrangian()
    {
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        //dimension of the problem
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

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
                KRATOS_THROW_ERROR(std::logic_error, "TotalLagrangian element does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
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
                KRATOS_THROW_ERROR(std::logic_error, "TotalLagrangian element does not support for integration points", GetProperties()[INTEGRATION_ORDER])
        }
        else
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

        //number of integration points used, mThisIntegrationMethod refers to the
        //integration method defined in the constructor
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        //Initialization of the constitutive law vector and
        // declaration, definition and initialization of the material
        // laws at each integration point
        if ( mConstitutiveLawVector.size() != integration_points.size() )
        {
            mConstitutiveLawVector.resize( integration_points.size() );
            InitializeMaterial( rCurrentProcessInfo );
        }

        if ( mIsInitialized )
            return;

        //initializing the Jacobian in the reference configuration
        GeometryType::JacobiansType J0;
        Matrix DeltaPosition(GetGeometry().size(), 3);

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
        }

        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        //calculating the domain size
        mTotalDomainInitialSize = 0.00;
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            //getting informations for integration
            double IntegrationWeight = integration_points[PointNumber].Weight();
            //calculating the total domain size
            mTotalDomainInitialSize += MathUtils<double>::Det(J0[PointNumber]) * IntegrationWeight;
        }
//        KRATOS_WATCH(mTotalDomainInitialSize)
        if ( mTotalDomainInitialSize < 0.0 )
        {
            std::stringstream ss;
            ss << "error on element -> " << this->Id() << std::endl;
            ss << ". Domain size can not be less than 0, mTotalDomainInitialSize = " << mTotalDomainInitialSize;
            KRATOS_THROW_ERROR( std::logic_error, ss.str(), "" );
        }
        this->SetValue(GEOMETRICAL_DOMAIN_SIZE, mTotalDomainInitialSize);

        mIsInitialized = true;

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector,
                                        const ProcessInfo& rCurrentProcessInfo,
                                        bool CalculateStiffnessMatrixFlag,
                                        bool CalculateResidualVectorFlag )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = this->GetStrainSize(dim);
        const unsigned int f_size = this->GetFSize(dim);

        Matrix B_Operator( strain_size, number_of_nodes * dim );

        Matrix F( f_size, f_size );

        Matrix D( strain_size, strain_size );

        Matrix C( f_size, f_size );

        Vector StrainVector( strain_size );

        Vector StressVector( strain_size );

        Vector N( number_of_nodes );

        Matrix DN_DX( number_of_nodes, dim );

        Matrix InvJ0(dim, dim);

        double DetJ0;

        double DetJn;

        double DetF;

        //constitutive law
        ConstitutiveLaw::Parameters const_params;
        const_params.SetStrainVector(StrainVector);
        const_params.SetDeformationGradientF(F);
        const_params.SetStressVector(StressVector);
        const_params.SetConstitutiveMatrix(D);
        const_params.SetProcessInfo(rCurrentProcessInfo);
        const_params.SetMaterialProperties(GetProperties());
        const_params.SetElementGeometry(GetGeometry());
        ConstitutiveLaw::StressMeasure stress_measure = ConstitutiveLaw::StressMeasure_PK2;

        //resizing as needed the LHS
        unsigned int MatSize = number_of_nodes * dim;

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

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        Matrix DeltaPosition(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

        //Current displacements
        Matrix CurrentDisp(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        //calculating jacobians
        MyJacobians Jacobians;
        this->CalculateJacobians(Jacobians, DeltaPosition, CurrentDisp);
        Jacobians.Check();

        const GeometryType::JacobiansType& J0 = Jacobians.J0();
        const GeometryType::JacobiansType& Jn = Jacobians.Jn();
        const GeometryType::JacobiansType& J = Jacobians.J();

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
            MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );
            noalias( N ) = row( Ncontainer, PointNumber );

            //deformation gradient
            this->CalculateF( F, N, DN_DX, CurrentDisp );
            DetF = MathUtils<double>::Det(F);
            const_params.SetDeterminantF(DetF);

            #ifdef CHECK_DEFORMATION_GRADIENT
            if (DetF < 0.0)
            {
                KRATOS_WATCH(F)
                KRATOS_ERROR << "Deformation gradient is negative at integration point " << PointNumber
                             << " of element " << Id();
            }
            #endif

            //strain calculation
            noalias( C ) = prod( trans( F ), F );

            this->CalculateStrain( C, StrainVector );
//            Comprobate_State_Vector( StrainVector ); // I don't understand why we need this

            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse( const_params, stress_measure );

            //calculating operator B
            this->CalculateB( B_Operator, F, N, DN_DX );

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight = this->GetIntegrationWeight(integration_points[PointNumber].Weight(), N) * DetJ0;

            if ( dim == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the lhs matrix is required
            {
                //contributions to stiffness matrix calculated on the reference config
                noalias( rLeftHandSideMatrix ) += prod( trans( B_Operator ), ( IntToReferenceWeight ) * Matrix( prod( D, B_Operator ) ) ); //to be optimized to remove the temporary
                CalculateAndAddKg( rLeftHandSideMatrix, N, DN_DX, StressVector, IntToReferenceWeight );
            }

            if ( CalculateResidualVectorFlag == true ) //calculation of the rhs vector is required
            {
                //contribution to external forces
                const Vector& BodyForce = GetProperties()[BODY_FORCE];

                // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
                CalculateAndAdd_ExtForceContribution( N, rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight );

                //contribution of gravity (if there is)
                AddBodyForcesToRHS( rRightHandSideVector, N, IntToReferenceWeight );

                // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
                noalias( rRightHandSideVector ) -= IntToReferenceWeight * prod( trans( B_Operator ), StressVector );
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //clean the internal data of the geometry
        GetGeometry().Clean();
        #endif

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::ApplyPrescribedDofs(const MatrixType& LHS_Contribution, VectorType& RHS_Constribution, const ProcessInfo& CurrentProcessInfo) const
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

    void TotalLagrangian::ComputePrescribedForces(const MatrixType& LHS_Contribution, VectorType& Force, const ProcessInfo& CurrentProcessInfo) const
    {
        // compute the force induced by the prescribed displacement.
        // Basically this function is identical to ApplyPrescribedDofs, but use PRESCRIBED_DISPLACEMENT instead of PRESCRIBED_DELTA_DISPLACEMENT
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int mat_size = dim * GetGeometry().size();
        for (unsigned int node = 0; node < GetGeometry().size(); ++node)
        {
            if (GetGeometry()[node].IsFixed(DISPLACEMENT_X))
            {
                double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DISPLACEMENT_X);
                if (temp != 0.0)
                    for (unsigned int i = 0; i < mat_size; ++i)
                        Force[i] -= LHS_Contribution(i, node * dim) * temp;
            }

            if (GetGeometry()[node].IsFixed(DISPLACEMENT_Y))
            {
                double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DISPLACEMENT_Y);
                if (temp != 0.0)
                    for (unsigned int i = 0; i < mat_size; ++i)
                        Force[i] -= LHS_Contribution(i, node * dim + 1) * temp;
            }

            if (dim > 2)
            {
                if (GetGeometry()[node].IsFixed(DISPLACEMENT_Z))
                {
                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DISPLACEMENT_Z);
                    if (temp != 0.0)
                        for (unsigned int i = 0; i < mat_size; ++i)
                            Force[i] -= LHS_Contribution(i, node * dim + 2) * temp;
                }
            }
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateJacobians( TotalLagrangian::MyJacobians& Jacobians,
            const Matrix& DeltaPosition, const Matrix& CurrentDisp ) const
    {
        // calculate the Jacobian in reference configuration no matter what the MoveMeshFlag is turned on or not
        Jacobians.J0Values = GetGeometry().Jacobian( Jacobians.J0Values, mThisIntegrationMethod, DeltaPosition );

        // calculate the Jacobian in current configuration no matter what the MoveMeshFlag is turned on or not
        // The coordinates for Jacobian are GetInitialPosition + CurrentDisp
        Matrix ShiftedDisp = DeltaPosition - CurrentDisp;
        Jacobians.JValues = GetGeometry().Jacobian( Jacobians.JValues, mThisIntegrationMethod, ShiftedDisp );

        // assign the correct pointers
        Jacobians.J0p = &Jacobians.J0Values;
        Jacobians.Jnp = &Jacobians.J0Values; // the Jacobians in previous configutation are taken as the reference ones
        Jacobians.Jp = &Jacobians.JValues;
    }

    void TotalLagrangian::CalculateJacobians( MatrixType& J0, MatrixType& J, const CoordinatesArrayType& rCoordinates,
            const Matrix& DeltaPosition, const Matrix& CurrentDisp ) const
    {
        // calculate the Jacobian in reference configuration no matter what the MoveMeshFlag is turned on or not
        J0 = GetGeometry().Jacobian( J0, rCoordinates, DeltaPosition );

        // calculate the Jacobian in current configuration no matter what the MoveMeshFlag is turned on or not
        Matrix ShiftedDisp = DeltaPosition - CurrentDisp;
        J = GetGeometry().Jacobian( J, rCoordinates, ShiftedDisp );
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;

        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::InitializeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
        int need_shape_function = 0, tmp;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            tmp = mConstitutiveLawVector[Point]->GetValue(IS_SHAPE_FUNCTION_REQUIRED, tmp);
            need_shape_function += tmp;
        }

        if (need_shape_function)
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->InitializeSolutionStep( GetProperties(), GetGeometry(), row(Ncontainer, Point), CurrentProcessInfo );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            GetGeometry().Clean();
            #endif
        }
        else
        {
            Vector dummy;
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->InitializeSolutionStep( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
            }
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::InitializeNonLinearIteration(const ProcessInfo& CurrentProcessInfo)
    {
        //reset all resistant forces at node
        for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
        {
            GetGeometry()[i].GetSolutionStepValue( REACTION_X ) = 0.0;
            GetGeometry()[i].GetSolutionStepValue( REACTION_Y ) = 0.0;
            GetGeometry()[i].GetSolutionStepValue( REACTION_Z ) = 0.0;
        }

        int need_shape_function = 0, tmp;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            tmp = mConstitutiveLawVector[Point]->GetValue(IS_SHAPE_FUNCTION_REQUIRED, tmp);
            need_shape_function += tmp;
        }

        if (need_shape_function)
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->InitializeNonLinearIteration( GetProperties(), GetGeometry(), row(Ncontainer, Point), CurrentProcessInfo );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            GetGeometry().Clean();
            #endif
        }
        else
        {
            Vector dummy;
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->InitializeNonLinearIteration( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
            }
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::FinalizeNonLinearIteration(const ProcessInfo& CurrentProcessInfo)
    {
        int need_shape_function = 0, tmp;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            tmp = mConstitutiveLawVector[Point]->GetValue(IS_SHAPE_FUNCTION_REQUIRED, tmp);
            need_shape_function += tmp;
        }

        if (need_shape_function)
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->FinalizeNonLinearIteration( GetProperties(), GetGeometry(), row(Ncontainer, Point), CurrentProcessInfo );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            GetGeometry().Clean();
            #endif
        }
        else
        {
            Vector dummy;
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->FinalizeNonLinearIteration( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
            }
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
        int need_shape_function = 0, tmp;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            tmp = mConstitutiveLawVector[Point]->GetValue(IS_SHAPE_FUNCTION_REQUIRED, tmp);
            need_shape_function += tmp;
        }

        if (need_shape_function)
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->FinalizeSolutionStep( GetProperties(), GetGeometry(), row(Ncontainer, Point), CurrentProcessInfo );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            GetGeometry().Clean();
            #endif
        }
        else
        {
            Vector dummy;
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->FinalizeSolutionStep( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
            }
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::InitializeMaterial( const ProcessInfo& CurrentProcessInfo )
    {
        KRATOS_TRY

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();

        int need_shape_function = 0, tmp;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            mConstitutiveLawVector[Point]->SetValue( PARENT_ELEMENT_ID, this->Id(), CurrentProcessInfo );
            mConstitutiveLawVector[Point]->SetValue( INTEGRATION_POINT_INDEX, Point, CurrentProcessInfo );
            tmp = mConstitutiveLawVector[Point]->GetValue(IS_SHAPE_FUNCTION_REQUIRED, tmp);
            need_shape_function += tmp;
        }

        if (need_shape_function)
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );

                //check constitutive law
                mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), CurrentProcessInfo );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif
        }
        else
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), Vector(1) );
            }
        }

        KRATOS_CATCH( "" )
    }

    void TotalLagrangian::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    inline void TotalLagrangian::AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, const double& Weight ) const
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        array_1d<double, 3> gravity;
        noalias( gravity ) = GetProperties()[GRAVITY];

        double density = 0.0;
        if( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
        {
            density = GetValue(DENSITY);
        }
        else
        {
            density = GetProperties()[DENSITY];
        }

        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                R( prim*dim + i ) += N_DISP( prim ) * density * gravity( i ) * Weight;
            }
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        const Vector& BodyForce,
        VectorType& rRightHandSideVector,
        double Weight
    ) const
    {
        KRATOS_TRY

        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            int index = dimension * i;

            for ( unsigned int j = 0; j < dimension; j++ ) rRightHandSideVector[index + j] += Weight * N[i] * BodyForce[j];
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateAndAddKg(
        MatrixType& K,
        const Vector& N,
        const Matrix& DN_DX,
        const Vector& StressVector,
        double Weight ) const
    {
        KRATOS_TRY

        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        Matrix StressTensor = MathUtils<double>::StressVectorToTensor( StressVector );
        Matrix ReducedKg = prod( DN_DX, Weight * Matrix( prod( StressTensor, trans( DN_DX ) ) ) ); //to be optimized
        MathUtils<double>::ExpandAndAddReducedMatrix( K, ReducedKg, dimension );

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateStrain(
        const Matrix& C,
        Vector& StrainVector ) const
    {
        KRATOS_TRY

        const unsigned int dimension = C.size1();

        if ( dimension == 2 )
        {
            StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

            StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

            StrainVector[2] = C( 0, 1 );
        }
        else if ( dimension == 3 )
        {
            StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

            StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

            StrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

            StrainVector[3] = C( 0, 1 ); // xy

            StrainVector[4] = C( 1, 2 ); // yz

            StrainVector[5] = C( 0, 2 ); // xz
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateF( Matrix& F, const Vector& N, const Matrix& DN_DX, const Matrix& CurrentDisp ) const
    {
        const unsigned int dim = F.size1();
        const unsigned int number_of_nodes = CurrentDisp.size1();

        for (unsigned int i = 0; i < dim; ++i)
        {
            for (unsigned int j = 0; j < dim; ++j)
            {
                F(i, j) = 0.0;
                for (unsigned int n = 0; n < number_of_nodes; ++n)
                    F(i, j) += DN_DX(n, j) * CurrentDisp(n, i);
            }
            F(i, i) += 1.0;
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateB(
        Matrix& B,
        const Matrix& F,
        const Matrix& DN_DX ) const
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = dimension * i;

            if ( dimension == 2 )
            {
                B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
                B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
                B( 1, index + 0 ) = F( 0, 1 ) * DN_DX( i, 1 );
                B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
                B( 2, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
                B( 2, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
            }
            else
            {
                B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
                B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
                B( 0, index + 2 ) = F( 2, 0 ) * DN_DX( i, 0 );
                B( 1, index + 0 ) = F( 0, 1 ) * DN_DX( i, 1 );
                B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
                B( 1, index + 2 ) = F( 2, 1 ) * DN_DX( i, 1 );
                B( 2, index + 0 ) = F( 0, 2 ) * DN_DX( i, 2 );
                B( 2, index + 1 ) = F( 1, 2 ) * DN_DX( i, 2 );
                B( 2, index + 2 ) = F( 2, 2 ) * DN_DX( i, 2 );
                B( 3, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
                B( 3, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
                B( 3, index + 2 ) = F( 2, 0 ) * DN_DX( i, 1 ) + F( 2, 1 ) * DN_DX( i, 0 );
                B( 4, index + 0 ) = F( 0, 1 ) * DN_DX( i, 2 ) + F( 0, 2 ) * DN_DX( i, 1 );
                B( 4, index + 1 ) = F( 1, 1 ) * DN_DX( i, 2 ) + F( 1, 2 ) * DN_DX( i, 1 );
                B( 4, index + 2 ) = F( 2, 1 ) * DN_DX( i, 2 ) + F( 2, 2 ) * DN_DX( i, 1 );
                B( 5, index + 0 ) = F( 0, 2 ) * DN_DX( i, 0 ) + F( 0, 0 ) * DN_DX( i, 2 );
                B( 5, index + 1 ) = F( 1, 2 ) * DN_DX( i, 0 ) + F( 1, 0 ) * DN_DX( i, 2 );
                B( 5, index + 2 ) = F( 2, 2 ) * DN_DX( i, 0 ) + F( 2, 0 ) * DN_DX( i, 2 );
            }
        }

        KRATOS_CATCH( "" )
    }



//************************************************************************************
//************************************************************************************

    void TotalLagrangian::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo ) const
    {
        int number_of_nodes = GetGeometry().size();
        int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int dim2 = number_of_nodes * dim;

        if ( rResult.size() != dim2 )
            rResult.resize( dim2, false );

        for ( int i = 0; i < number_of_nodes; i++ )
        {
            int index = i * dim;
            rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

            if ( dim == 3 )
                rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
        }

    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo ) const
    {
        ElementalDofList.resize( 0 );
        int dim = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

            if ( dim == 3 )
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
            }
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int MatSize = dimension * NumberOfNodes;

        if ( rMassMatrix.size1() != MatSize )
            rMassMatrix.resize( MatSize, MatSize, false );

        noalias( rMassMatrix ) = ZeroMatrix( MatSize, MatSize );

        double density;
        if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
            density = GetValue(DENSITY);
        else
            density = GetProperties()[DENSITY];

        /// Lumped mass

        // double TotalMass = mTotalDomainInitialSize * density;
        // if ( dimension == 2 ) TotalMass *= GetProperties()[THICKNESS];

        // Vector LumpFact;

        // LumpFact = GetGeometry().LumpingFactors( LumpFact );

        // for ( unsigned int i = 0; i < NumberOfNodes; i++ )
        // {
        //     double temp = LumpFact[i] * TotalMass;

        //     for ( unsigned int j = 0; j < dimension; j++ )
        //     {
        //         unsigned int index = i * dimension + j;
        //         rMassMatrix( index, index ) = temp;
        //     }
        // }

        /// Consistent mass

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        Matrix DeltaPosition(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

        //Current displacements
        Matrix CurrentDisp(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        //Jacobians
        Vector N(GetGeometry().size());
        GeometryType::JacobiansType J0;
        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );
        double DetJ0;

        //Mass contribution
        for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            noalias( N ) = row( Ncontainer, PointNumber );
            DetJ0 = MathUtils<double>::Det(J0[PointNumber]);

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight = density * this->GetIntegrationWeight(integration_points[PointNumber].Weight(), N) * DetJ0;

            //modify integration weight in case of 2D
            if ( dimension == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

            for ( unsigned int i = 0; i < NumberOfNodes; ++i )
            {
                for ( unsigned int j = 0; j < NumberOfNodes; ++j )
                {
                    for ( unsigned int k = 0; k < dimension; ++k )
                        rMassMatrix(dimension*i + k, dimension*j + k) += N(i) * N(j) * IntToReferenceWeight;
                }
            }
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateDampingMatrix( MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int MatSize = number_of_nodes * dim;

        if ( rDampingMatrix.size1() != MatSize )
            rDampingMatrix.resize( MatSize, MatSize, false );

        // no damping
        noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );

        KRATOS_CATCH( "" )
    }

    TotalLagrangian::IntegrationMethod TotalLagrangian::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::SetValuesOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        //   std::cout << mConstitutiveLawVector[0] << std::endl;
        //  if (rVariable==INSITU_STRESS)
        //  {
        //                    for( unsigned int PointNumber = 0; PointNumber<GetGeometry().IntegrationPoints(mThisIntegrationMethod).size(); PointNumber++ )
        //   {
        //    mConstitutiveLawVector[PointNumber]->SetValue(INSITU_STRESS, rValues[PointNumber],
        //      rCurrentProcessInfo );
        //   }
        //  }
        //                if (rVariable==MATERIAL_PARAMETERS)
        //                {
        //                    for( unsigned int PointNumber = 0; PointNumber<GetGeometry().IntegrationPoints(mThisIntegrationMethod).size(); PointNumber++ )
        //                    {
        //                        mConstitutiveLawVector[PointNumber]->SetValue( MATERIAL_PARAMETERS,
        //                                rValues[PointNumber], rCurrentProcessInfo );
        //                    }
        //                }

        for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
                    rValues[PointNumber], rCurrentProcessInfo );
        }

    }


//************************************************************************************
//************************************************************************************

    void TotalLagrangian::SetValuesOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
                    rValues[PointNumber], rCurrentProcessInfo );
        }

    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateOnIntegrationPoints( const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo )
    {
        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        if ( rValues.size() != integration_points.size() )
            rValues.resize( integration_points.size(), false );

        if ( rVariable == CURRENT_DEFORMATION_GRADIENT_DETERMINANT )
        {
            const unsigned int dim = GetGeometry().WorkingSpaceDimension();

            Matrix F( dim, dim );

            Matrix InvJ0(dim, dim);

            double DetJ0;

            double DetF;

            //initializing the Jacobian in the reference configuration
            Matrix DeltaPosition(GetGeometry().size(), 3);
            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
                noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

            //Current displacements
            Matrix CurrentDisp(GetGeometry().size(), 3);
            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
                noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

            //calculating jacobians
            MyJacobians Jacobians;
            this->CalculateJacobians(Jacobians, DeltaPosition, CurrentDisp);
            Jacobians.Check();

            const GeometryType::JacobiansType& J0 = Jacobians.J0();
            const GeometryType::JacobiansType& Jn = Jacobians.Jn();
            const GeometryType::JacobiansType& J = Jacobians.J();

            for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
            {
                //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
                MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );

                //deformation gradient
                noalias( F ) = prod( J[PointNumber], InvJ0 );
                DetF = MathUtils<double>::Det(F);

                rValues[PointNumber] = DetF;
            }

            return;
        }
        else
        {
            for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
                rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
            std::vector<array_1d<double, 3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        if( rVariable == DISPLACEMENT )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                noalias(rValues[point]) = ZeroVector(3);
                for(std::size_t i = 0; i < GetGeometry().size(); ++i)
                {
                    const array_1d<double, 3>& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                    noalias(rValues[point]) += Ncontainer(point, i) * displacement;
                }
            }

            return;
        }

        if( rVariable == INTEGRATION_POINT_GLOBAL )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                rValues[point] = GetGeometry().GlobalCoordinates(rValues[point], integration_points[point]);
            }

            return;
        }

        if( rVariable == INTEGRATION_POINT_LOCAL )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                noalias(rValues[point]) = integration_points[point];
            }

            return;
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Clean();
        #endif
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int size = GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = this->GetStrainSize(dim);

        if ( rValues.size() != size )
            rValues.resize( size );

        if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
                    PointNumber++ )
            {
                if ( rValues[PointNumber].size() != strain_size )
                    rValues[PointNumber].resize( strain_size, false );

                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
            }
        }
        else if ( rVariable == ELASTIC_STRAIN_VECTOR || rVariable == PLASTIC_STRAIN_VECTOR )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != strain_size )
                    rValues[i].resize( strain_size );
                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
            }
        }
        else if ( rVariable == STRESSES )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != strain_size )
                    rValues[i].resize( strain_size );
                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( STRESSES, rValues[i] );
            }
        }
        else if ( rVariable == MATERIAL_PARAMETERS )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
            {
                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( MATERIAL_PARAMETERS, rValues[PointNumber] );
            }
        }
        else if ( rVariable == INTERNAL_VARIABLES )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
                    PointNumber++ )
            {
                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( INTERNAL_VARIABLES, rValues[PointNumber] );

            }
        }
        else
        {
            if ( rValues.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
                rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );

            for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
                rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = this->GetStrainSize(dim);
        const unsigned int f_size = this->GetFSize(dim);

        Matrix F( f_size, f_size );

        Matrix D( strain_size, strain_size );

        Matrix strain_tensor( f_size, f_size );

        Matrix stress_tensor( f_size, f_size );

        Vector StrainVector( strain_size );

        Vector StressVector( strain_size );

        Matrix DN_DX( number_of_nodes, dim );

        Vector N( number_of_nodes );

        Matrix InvJ0(dim, dim);

        double DetJ0;

        const Matrix eye = IdentityMatrix(f_size);

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        Matrix DeltaPosition(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

        //Current displacements
        Matrix CurrentDisp(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        //calculating jacobians
        MyJacobians Jacobians;
        this->CalculateJacobians(Jacobians, DeltaPosition, CurrentDisp);
        Jacobians.Check();

        const GeometryType::JacobiansType& J0 = Jacobians.J0();
        const GeometryType::JacobiansType& Jn = Jacobians.Jn();
        const GeometryType::JacobiansType& J = Jacobians.J();

        if ( rValues.size() != integration_points.size() )
            rValues.resize( integration_points.size() );

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            //deformation gradient
            MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );
            noalias( N ) = row( Ncontainer, PointNumber );

            //deformation gradient
            this->CalculateF( F, N, DN_DX, CurrentDisp );

            if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            {
                if ( rValues[PointNumber].size1() != f_size || rValues[PointNumber].size2() != f_size )
                    rValues[PointNumber].resize( f_size, f_size, false );

                //strain calculation
                noalias( strain_tensor ) = 0.5 * (prod( trans( F ), F ) - eye);

                noalias(rValues[PointNumber]) = strain_tensor;
            }
            else if ( rVariable == LEFT_STRETCH_TENSOR || rVariable == RIGHT_STRETCH_TENSOR )
            {
                if ( rValues[PointNumber].size1() != f_size || rValues[PointNumber].size2() != f_size )
                    rValues[PointNumber].resize( f_size, f_size, false );

                Matrix F3d(3, 3);
                if (f_size == 2)
                {
                    for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                            F3d(i, j) = F(i, j);
                    F3d(2, 2) = 1.0;
                }
                else if (f_size == 3)
                    noalias(F3d) = F;

                //strain calculation
                Matrix Aux( 3, 3 );
                if (rVariable == RIGHT_STRETCH_TENSOR)
                    noalias( Aux ) = prod( trans( F3d ), F3d );
                else if (rVariable == LEFT_STRETCH_TENSOR)
                    noalias( Aux ) = prod( F3d, trans( F3d ) );

                //polar decomposition
                std::vector<double> pri(3);
                std::vector<Matrix> eigprj(3);
                EigenUtility::calculate_principle_stresses( Aux(0, 0), Aux(1, 1), Aux(2, 2),
                        Aux(0, 1), Aux(1, 2), Aux(0, 2),
                        pri[0], pri[1], pri[2], eigprj[0], eigprj[1], eigprj[2]);

                Matrix stretch_tensor(3, 3);
                SD_MathUtils<double>::ComputeIsotropicTensorFunction(EigenUtility::sqrt, stretch_tensor, pri, eigprj);

                if (f_size == 2)
                {
                    Matrix& value = rValues[PointNumber];
                    for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                            value(i, j) = stretch_tensor(i, j);
                }
                else if (f_size == 3)
                    noalias(rValues[PointNumber]) = stretch_tensor;
            }
            else if ( rVariable == CURRENT_DEFORMATION_GRADIENT )
            {
                if ( rValues[PointNumber].size1() != f_size || rValues[PointNumber].size2() != f_size )
                    rValues[PointNumber].resize( f_size, f_size, false );

                noalias(rValues[PointNumber]) = F;
            }
            else if ( rVariable == PK2_STRESS_TENSOR )
            {
                if ( rValues[PointNumber].size1() != f_size || rValues[PointNumber].size2() != f_size )
                    rValues[PointNumber].resize( f_size, f_size, false );

                mConstitutiveLawVector[PointNumber]->GetValue( rVariable, stress_tensor );

                noalias(rValues[PointNumber]) = stress_tensor;
            }
            else if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
            {
                if ( rValues[PointNumber].size1() != f_size || rValues[PointNumber].size2() != f_size )
                    rValues[PointNumber].resize( f_size, f_size, false );

                mConstitutiveLawVector[PointNumber]->GetValue( rVariable, stress_tensor );

                noalias(rValues[PointNumber]) = stress_tensor;
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //clean the internal data of the geometry
        GetGeometry().Clean();
        #endif

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::GetValuesVector( Vector& values, int Step ) const
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
        }
    }


//************************************************************************************
//************************************************************************************

    void TotalLagrangian::GetFirstDerivativesVector( Vector& values, int Step ) const
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT_Z, Step );
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::GetSecondDerivativesVector( Vector& values, int Step ) const
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo )
    {

        double lamda = 1.00; // parametro que depende del tipo de problema y del elemento pag 308 libro dinamica de Barbat
        double c1 = 0.00; //sqrt(GetProperties()[YOUNG_MODULUS]/GetProperties()[DENSITY]); velocidad del sonido en el medio
        double c2 = 0.00; // norma de la velocidad actual dentro del elemento
        double c = 0.00;
        double wmax = 0.00;
        Vector Values( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );
        Vector Velocities;

        GetFirstDerivativesVector( Velocities, 0 );

        if ( rVariable == DELTA_TIME )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
                    PointNumber++ )
            {
                mConstitutiveLawVector[PointNumber]-> GetValue( DELTA_TIME, c1 );
                Values[PointNumber] = c1;
            }
        }

        c1 = ( *std::max_element( Values.begin(), Values.end() ) );

        c2 = norm_2( Velocities );

        c = ( c1 > c2 ) ? c1 : c2;


        double le = GetGeometry().Length();
        //KRATOS_WATCH(le)

        /// maxima frecuencia de un elemento
        wmax = ( lamda * c ) / le;
        Output = 2.0 / wmax;
        //KRATOS_WATCH(Output)

    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangian::Comprobate_State_Vector( Vector& Result ) const
    {
        for ( unsigned int i = 0.00; i < Result.size(); i++ )
        {
            if ( fabs( Result( i ) ) < 1E-9 )
            {
                Result( i ) = 0.00;
            }
        }
    }


//************************************************************************************
//************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int  TotalLagrangian::Check( const ProcessInfo& rCurrentProcessInfo ) const
    {
        KRATOS_TRY

        unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

        //verify that the variables are correctly initialized

//        if ( VELOCITY.Key() == 0 )
//            KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" );

//        if ( DISPLACEMENT.Key() == 0 )
//            KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" );

//        if ( ACCELERATION.Key() == 0 )
//            KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" );

//        if ( DENSITY.Key() == 0 )
//            KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "" );

//        if ( BODY_FORCE.Key() == 0 )
//            KRATOS_THROW_ERROR( std::invalid_argument, "BODY_FORCE has Key zero! (check if the application is correctly registered", "" );

//        if ( THICKNESS.Key() == 0 )
//            KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" );

//        //verify that the dofs exist
//        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
//        {
//            if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
//                KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() );

//            if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
//                KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() );
//        }

        //verify that the constitutive law exists
        if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() );
        }

        // verify the strain measure
        ConstitutiveLaw::Features features;
        this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(features);
        if ( std::find(features.GetStrainMeasures().begin(), features.GetStrainMeasures().end(), ConstitutiveLaw::StrainMeasure_Infinitesimal) == features.GetStrainMeasures().end()
          && std::find(features.GetStrainMeasures().begin(), features.GetStrainMeasures().end(), ConstitutiveLaw::StrainMeasure_GreenLagrange) == features.GetStrainMeasures().end() )
        {
            std::stringstream ss;
            ss << "The constitutive law strain measures are not supported by this element";
            KRATOS_THROW_ERROR( std::logic_error, ss.str(), "" )
        }

        //Verify that the body force is defined
        if ( this->GetProperties().Has( BODY_FORCE ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "BODY_FORCE not provided for property ", this->GetProperties().Id() )
        }



        //check constitutive law
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
        }

        //check if it is in the XY plane for 2D case


        return 0;

        KRATOS_CATCH( "" );
    }


    void TotalLagrangian::save( Serializer& rSerializer ) const
    {
//  std::cout << "Saving the TotalLagrangian #" << Id() << std::endl;
        rSerializer.save( "Name", "TotalLagrangian" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    }

    void TotalLagrangian::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
//  std::cout << "Loading the TotalLagrangian #" << Id() << std::endl;
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    }

} // Namespace Kratos

#undef CHECK_DEFORMATION_GRADIENT
