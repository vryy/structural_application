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
//   Modified by:         $Author: hbui $
//   Date:                $Date: 9 Aug 2021 $
//   Revision:            $Revision: 1.27 $
//
//


// System includes

// External includes


// Project includes
#include "utilities/math_utils.h"
#include "custom_elements/finite_strain.h"
#include "custom_utilities/geometry_utility.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"

// #define ENABLE_DEBUG_CONSTITUTIVE_LAW
#define CHECK_DEFORMATION_GRADIENT
#define USE_DETERMINANT_F0_FOR_FBAR

#ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
#include <iomanip>
#endif

namespace Kratos
{

    FiniteStrain::FiniteStrain( IndexType NewId, GeometryType::Pointer pGeometry )
            : Element( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
        mIsInitialized = false;
    }

//************************************************************************************
//************************************************************************************

    FiniteStrain::FiniteStrain( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : Element( NewId, pGeometry, pProperties )
    {
        //         const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

        mIsInitialized = false;
    }

    Element::Pointer FiniteStrain::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new FiniteStrain( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    Element::Pointer FiniteStrain::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new FiniteStrain( NewId, pGeom, pProperties ) );
    }

    FiniteStrain::~FiniteStrain()
    {
    }

//************************************************************************************
//************************************************************************************

    void FiniteStrain::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        //dimension of the problem
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();

        // integration rule
        if(this->Has( INTEGRATION_ORDER ))
        {
            if(this->GetValue(INTEGRATION_ORDER) == 1)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 2)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 3)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 4)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 5)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "FiniteStrain element does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
        }
        else if(GetProperties().Has( INTEGRATION_ORDER ))
        {
            if(GetProperties()[INTEGRATION_ORDER] == 1)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 2)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 3)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 4)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 5)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "FiniteStrain element does not support for integration points", GetProperties()[INTEGRATION_ORDER])
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
            InitializeMaterial();
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

    void FiniteStrain::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                     VectorType& rRightHandSideVector,
                                     const ProcessInfo& rCurrentProcessInfo,
                                     bool CalculateStiffnessMatrixFlag,
                                     bool CalculateResidualVectorFlag )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = this->GetStrainSize(dim);
        const unsigned int g_size = this->GetGSize(dim);
        const unsigned int f_size = this->GetFSize(dim);

        Matrix B;
        if (CalculateResidualVectorFlag)
        {
            B.resize( strain_size, number_of_nodes * dim, false );
        }

        Matrix Gx;
        Matrix GX( g_size, number_of_nodes*dim );
        if (CalculateStiffnessMatrixFlag)
        {
            Gx.resize( g_size, number_of_nodes * dim, false );
        }

        Matrix F( f_size, f_size );

        Matrix InvF( f_size, f_size );

        Matrix A( g_size, g_size );

        Matrix B( dim, dim );

        Vector StrainVector( strain_size );

        Vector StressVector( strain_size );

        Matrix DN_DX( number_of_nodes, dim );
        Vector N( number_of_nodes );
        Matrix DN_Dx( number_of_nodes, dim );

        Matrix CurrentDisp( number_of_nodes, 3 );

        Matrix InvJ0(dim, dim), InvJ(dim, dim);

        double DetJ0, DetJ;

        double DetF;

        #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
        if (Id() == DEBUG_ELEMENT_ID)
        {
            Vector disp_e;
            this->GetValuesVector(disp_e, 0);
            KRATOS_WATCH(disp_e)
        }
        std::cout << std::setprecision(16) << std::scientific;
        #endif

        //constitutive law
        ConstitutiveLaw::Parameters const_params;
        const_params.SetStrainVector(StrainVector);
        const_params.SetDeformationGradientF(F);
        const_params.SetStressVector(StressVector);
        const_params.SetConstitutiveMatrix(A);
        const_params.SetProcessInfo(rCurrentProcessInfo);
        const_params.SetMaterialProperties(GetProperties());
        const_params.SetElementGeometry(GetGeometry());
        ConstitutiveLaw::StressMeasure stress_measure = ConstitutiveLaw::StressMeasure_Cauchy;

        //resizing as needed the LHS
        unsigned int MatSize = number_of_nodes * dim;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != MatSize || rLeftHandSideMatrix.size2() != MatSize )
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
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        GeometryType::JacobiansType J0;
        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        //calculating actual jacobian
        GeometryType::JacobiansType J;
        Matrix ShiftedDisp = DeltaPosition - CurrentDisp;
        J = GetGeometry().Jacobian( J, mThisIntegrationMethod, ShiftedDisp );

        /////////////////////////////////////////////////////////////////////////
        //// Compute FO, GO for Fbar formulation
        //// Reference: Souze de Neto, Computational Plasticity, Box 15.1
        /////////////////////////////////////////////////////////////////////////
        const bool is_fbar = GetProperties().Has(IS_FBAR) ? GetProperties()[IS_FBAR] : false;
        Matrix FO, GOx, Q, A3d, eye, stress_tensor;
        SD_MathUtils<double>::Fourth_Order_Tensor Atensor, Qtensor, ItimesI;
        double DetFO;
        if (is_fbar)
        {
            CoordinatesArrayType origin;
            GeometryUtility::ComputeOrigin(GetGeometry().GetGeometryFamily(), origin);

            Matrix J0O, JO;
            J0O = GetGeometry().Jacobian( J0O, origin, DeltaPosition );
            JO = GetGeometry().Jacobian( JO, origin, ShiftedDisp );

            Vector NO( number_of_nodes );
            NO = GetGeometry().ShapeFunctionsValues(NO, origin);

            Matrix DN_De_O( number_of_nodes, dim );
            DN_De_O = GetGeometry().ShapeFunctionsLocalGradients(DN_De_O, origin);

            Matrix InvJ0O(dim, dim), InvJO(dim, dim);
            double DetJ0O, DetJO;
            Matrix DN_Dx_O( number_of_nodes, dim );

            MathUtils<double>::InvertMatrix( J0O, InvJ0O, DetJ0O );
            MathUtils<double>::InvertMatrix( JO, InvJO, DetJO );
            noalias( DN_Dx_O ) = prod( DN_De_O, InvJO );

            // FO
            Matrix DN_DX_O( number_of_nodes, dim );
            noalias( DN_DX_O ) = prod( DN_De_O, InvJ0O );

            Matrix GOX(g_size, number_of_nodes * dim);
            this->CalculateG( GOX, NO, DN_DX_O );

            FO.resize(f_size, f_size, false);
            this->CalculateF( FO, GOX, CurrentDisp );
            DetFO = MathUtils<double>::Det(FO);

            #ifdef CHECK_DEFORMATION_GRADIENT
            if (DetFO < 0.0)
            {
                KRATOS_WATCH(origin)
                KRATOS_WATCH(J0O)
                KRATOS_WATCH(InvJ0O)
                KRATOS_WATCH(JO)
                KRATOS_WATCH(FO)
                KRATOS_WATCH(DetFO)
                KRATOS_ERROR << "Deformation gradient is negative at origin of element " << Id();
            }
            #endif

            // GOx
            GOx.resize( g_size, number_of_nodes * dim, false );
            this->CalculateG( GOx, NO, DN_Dx_O, CurrentDisp );

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
            if (Id() == DEBUG_ELEMENT_ID)
            {
                KRATOS_WATCH(GOx)
            }
            #endif

            // resize Q as needed
            Q.resize(g_size, g_size, false);

            // resize A3d as needed
            A3d.resize(9, 9, false);

            // initialize the tensors
            SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Atensor);
            SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Qtensor);

            eye = IdentityMatrix(3);
            if (f_size == 2) eye(2, 2) = 0.0; // this modification is required because we don't want to involve the zz component
            SD_MathUtils<double>::CalculateFourthOrderZeroTensor(ItimesI);
            SD_MathUtils<double>::OuterProductFourthOrderTensor(1.0, eye, eye, ItimesI);

            stress_tensor.resize(3, 3);
        }

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
            MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            MathUtils<double>::InvertMatrix( J[PointNumber], InvJ, DetJ );
            noalias( N ) = row( Ncontainer, PointNumber );
            noalias( DN_Dx ) = prod( DN_De[PointNumber], InvJ );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );

            //deformation gradient
            this->CalculateG( GX, N, DN_DX );
            this->CalculateF( F, GX, CurrentDisp );


            DetF = MathUtils<double>::Det(F);

            #ifdef CHECK_DEFORMATION_GRADIENT
            if (DetF < 0.0)
            {
                KRATOS_WATCH(F)
                KRATOS_ERROR << "Deformation gradient is negative at integration point " << PointNumber
                             << " of element " << Id();
            }
            #endif

            if (is_fbar)
            {
                if (f_size == 2) // plane strain
                {
                    F *= std::sqrt(DetFO/DetF);
                    #ifdef USE_DETERMINANT_FO_FOR_FBAR
                    const_params.SetDeterminantF(DetFO); // for the reason why to do this, see iffba2.f
                    #else
                    const_params.SetDeterminantF(std::sqrt(DetFO/DetF)*DetF);
                    #endif
                }
                else if (f_size == 3) // plane stress, axisymmetric, 3D
                {
                    F *= std::cbrt(DetFO/DetF);
                    #ifdef USE_DETERMINANT_FO_FOR_FBAR
                    const_params.SetDeterminantF(DetFO); // for the reason why to do this, see iffba2.f
                    #else
                    const_params.SetDeterminantF(std::cbrt(DetFO/DetF)*DetF);
                    #endif
                }
                else
                {
                    KRATOS_ERROR << "Invalid size " << f_size << " of deformation gradient";
                }
            }
            else
                const_params.SetDeterminantF(DetF);

            //strain calculation
            noalias( C ) = prod( trans( F ), F );

            CalculateStrain( C, StrainVector );

            //integrate the material
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse( const_params, stress_measure );

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight = this->GetIntegrationWeight(integration_points[PointNumber].Weight(), N, CurrentDisp) * DetJ;

            if ( dim == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

            if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
            {
                //calculating operator B
                CalculateB( B_Operator, N, DN_Dx, CurrentDisp );

                //contribution to external forces
                const Vector& BodyForce = GetProperties()[BODY_FORCE];
                CalculateAndAdd_ExtForceContribution( N, rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight );

                //contribution of gravity (if there is)
                AddBodyForcesToRHS( rRightHandSideVector, N, IntToReferenceWeight );

                //contribution of internal forces
                noalias( rRightHandSideVector ) -= IntToReferenceWeight * prod( trans( B_Operator ), StressVector );
            }

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                //calculating operator G
                this->CalculateG( Gx, N, DN_Dx, CurrentDisp );

                //contributions to stiffness
                noalias( rLeftHandSideMatrix ) += prod( trans( Gx ), ( IntToReferenceWeight ) * Matrix( prod( A, Gx ) ) ); //to be optimized to remove the temporary

                if (is_fbar)
                {
                    //compute Q matrix
                    mConstitutiveLawVector[PointNumber]->GetValue(CAUCHY_STRESS_TENSOR, stress_tensor);
                    if ( dim == 2)
                    {
                        mConstitutiveLawVector[PointNumber]->GetValue(THREED_ALGORITHMIC_TANGENT, A3d);
                        SD_MathUtils<double>::UnsymmetricMatrixToTensor(A3d, Atensor);
                    }
                    else if (dim == 3)
                    {
                        SD_MathUtils<double>::UnsymmetricMatrixToTensor(A, Atensor);
                    }

                    SD_MathUtils<double>::ZeroFourthOrderTensor(Qtensor);
                    if (f_size == 2) // plane strain
                    {
                        SD_MathUtils<double>::ProductFourthOrderTensor(0.5, Atensor, ItimesI, Qtensor);
                        SD_MathUtils<double>::OuterProductFourthOrderTensor(-0.5, stress_tensor, eye, Qtensor);
                    }
                    else if (f_size == 3) // plane stress, axisymmetric, 3D
                    {
                        SD_MathUtils<double>::ProductFourthOrderTensor(1.0/3, Atensor, ItimesI, Qtensor);
                        SD_MathUtils<double>::OuterProductFourthOrderTensor(-2.0/3, stress_tensor, eye, Qtensor);
                    }
                    else
                        KRATOS_ERROR << "Invalid size " << f_size << " of deformation gradient";

                    SD_MathUtils<double>::TensorToUnsymmetricMatrix(Qtensor, Q);

                    //contributions to stiffness
                    noalias( rLeftHandSideMatrix ) += prod( trans( Gx ), ( IntToReferenceWeight ) * Matrix( prod( Q, GOx - Gx ) ) );
                }
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

    void FiniteStrain::ApplyPrescribedDofs(const MatrixType& LHS_Contribution, VectorType& RHS_Constribution, const ProcessInfo& CurrentProcessInfo) const
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

//************************************************************************************
//************************************************************************************

    void FiniteStrain::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

//************************************************************************************
//************************************************************************************

    void FiniteStrain::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;

        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

//************************************************************************************
//************************************************************************************

    double FiniteStrain::CalculateIntegrationWeight( const double& GaussPointWeight, const double& DetJ )
    {
        //to permorm the integration over the reference domain we need to include
        // the thickness in 2D
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        double weight = GaussPointWeight;

        weight *= DetJ;

        if ( dimension == 2 ) weight *= GetProperties()[THICKNESS];

        return weight;
    }

//************************************************************************************
//************************************************************************************

    void FiniteStrain::InitializeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
        // check if the constitutive law need current deformation gradient
        bool need_current_deformation_gradient = false;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            if (mConstitutiveLawVector[Point]->Has(CURRENT_DEFORMATION_GRADIENT))
            {
                need_current_deformation_gradient = true;
                break;
            }
        }

        if ( need_current_deformation_gradient )
        {
            std::vector<MatrixType> Values;
            this->CalculateOnIntegrationPoints( CURRENT_DEFORMATION_GRADIENT, Values, CurrentProcessInfo );
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->SetValue( CURRENT_DEFORMATION_GRADIENT, Values[Point], CurrentProcessInfo );
            }
        }

        // check if the constitutive law needs shape function
        bool need_shape_function = false;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            int tmp;
            tmp = mConstitutiveLawVector[Point]->GetValue(IS_SHAPE_FUNCTION_REQUIRED, tmp);
            if (tmp)
            {
                need_shape_function = true;
                break;
            }
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

    void FiniteStrain::InitializeNonLinearIteration(const ProcessInfo& CurrentProcessInfo)
    {
        //reset all resistant forces at node
        for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
        {
            GetGeometry()[i].GetSolutionStepValue( REACTION_X ) = 0.0;
            GetGeometry()[i].GetSolutionStepValue( REACTION_Y ) = 0.0;
            GetGeometry()[i].GetSolutionStepValue( REACTION_Z ) = 0.0;
        }

        // check if the constitutive law needs shape function
        bool need_shape_function = false;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            int tmp;
            tmp = mConstitutiveLawVector[Point]->GetValue(IS_SHAPE_FUNCTION_REQUIRED, tmp);
            if (tmp)
            {
                need_shape_function = true;
                break;
            }
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

    void FiniteStrain::FinalizeNonLinearIteration(const ProcessInfo& CurrentProcessInfo)
    {
        // check if the constitutive law need current deformation gradient
        bool need_current_deformation_gradient = false;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            if (mConstitutiveLawVector[Point]->Has(CURRENT_DEFORMATION_GRADIENT))
            {
                need_current_deformation_gradient = true;
                break;
            }
        }

        if ( need_current_deformation_gradient )
        {
            std::vector<MatrixType> Values;
            this->CalculateOnIntegrationPoints( CURRENT_DEFORMATION_GRADIENT, Values, CurrentProcessInfo );
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->SetValue( CURRENT_DEFORMATION_GRADIENT, Values[Point], CurrentProcessInfo );
            }

            const bool is_fbar = GetProperties().Has(IS_FBAR) ? GetProperties()[IS_FBAR] : false;
            if (is_fbar)
            {
                std::vector<double> Values2;
                this->CalculateOnIntegrationPoints( CURRENT_DEFORMATION_GRADIENT_DETERMINANT, Values2, CurrentProcessInfo );
                for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
                {
                    mConstitutiveLawVector[Point]->SetValue( CURRENT_DEFORMATION_GRADIENT_DETERMINANT, Values2[Point], CurrentProcessInfo );
                }
            }
        }

        // check if the constitutive law needs shape function
        bool need_shape_function = false;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            int tmp;
            tmp = mConstitutiveLawVector[Point]->GetValue(IS_SHAPE_FUNCTION_REQUIRED, tmp);
            if (tmp)
            {
                need_shape_function = true;
                break;
            }
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

    void FiniteStrain::FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
        // check if the constitutive law need current deformation gradient
        bool need_current_deformation_gradient = false;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            if (mConstitutiveLawVector[Point]->Has(CURRENT_DEFORMATION_GRADIENT))
            {
                need_current_deformation_gradient = true;
                break;
            }
        }

        if ( need_current_deformation_gradient )
        {
            std::vector<MatrixType> Values;
            this->CalculateOnIntegrationPoints( CURRENT_DEFORMATION_GRADIENT, Values, CurrentProcessInfo );
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->SetValue( CURRENT_DEFORMATION_GRADIENT, Values[Point], CurrentProcessInfo );
            }
        }

        // check if the constitutive law needs shape function
        bool need_shape_function = false;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            int tmp;
            tmp = mConstitutiveLawVector[Point]->GetValue(IS_SHAPE_FUNCTION_REQUIRED, tmp);
            if (tmp)
            {
                need_shape_function = true;
                break;
            }
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

    void FiniteStrain::InitializeMaterial()
    {
        KRATOS_TRY

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();

        int need_shape_function = 0, tmp;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            mConstitutiveLawVector[Point]->SetValue( PARENT_ELEMENT_ID, this->Id(), *(ProcessInfo*)0);
            mConstitutiveLawVector[Point]->SetValue( INTEGRATION_POINT_INDEX, Point, *(ProcessInfo*)0);
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
                mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), *(ProcessInfo*)0 );
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

    void FiniteStrain::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    inline void FiniteStrain::AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, const double& Weight )
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

    inline void FiniteStrain::CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        const Vector& BodyForce,
        VectorType& rRightHandSideVector,
        const double& weight
    )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            int index = dimension * i;

            for ( unsigned int j = 0; j < dimension; j++ ) rRightHandSideVector[index + j] += weight * N[i] * BodyForce[j];
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void FiniteStrain::CalculateF( Matrix& F, const Matrix& G_Operator, const Matrix& CurrentDisp ) const
    {
        const unsigned int dim = F.size1();
        const unsigned int number_of_nodes = CurrentDisp.size1();

        for (unsigned int i = 0; i < dim; ++i)
        {
            for (unsigned int j = 0; j < dim; ++j)
            {
                F(i, j) = 0.0;
                for (unsigned int n = 0; n < number_of_nodes; ++n)
                    F(i, j) += G_Operator(j*dim, n*dim) * CurrentDisp(n, i);
            }
            F(i, i) += 1.0;
        }
    }

//************************************************************************************
//************************************************************************************

    void FiniteStrain::CalculateStrain(
        const Matrix& C,
        Vector& StrainVector ) const
    {
        KRATOS_TRY

        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if ( dimension == 2 )
        {
            if ( StrainVector.size() != 3 ) StrainVector.resize( 3, false );

            StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

            StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

            StrainVector[2] = C( 0, 1 );
        }

        if ( dimension == 3 )
        {
            if ( StrainVector.size() != 6 ) StrainVector.resize( 6, false );

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

    void FiniteStrain::CalculateB( Matrix& B_Operator, const VectorType& N, const Matrix& DN_DX, const Matrix& CurrentDisp ) const
    {
        KRATOS_TRY

        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_nodes = GetGeometry().size();

        B_Operator.clear();

        if(dim == 2)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                B_Operator( 0, i*2     ) = DN_DX( i, 0 );
                B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 );
                B_Operator( 2, i*2     ) = DN_DX( i, 1 );
                B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
            }
        }
        else if(dim == 3)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                B_Operator( 0, i*3     ) = DN_DX( i, 0 );
                B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
                B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 );
                B_Operator( 3, i*3     ) = DN_DX( i, 1 );
                B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
                B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
                B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
                B_Operator( 5, i*3     ) = DN_DX( i, 2 );
                B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
            }
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void FiniteStrain::CalculateG( Matrix& G_Operator, const VectorType& N, const Matrix& DN_DX ) const
    {
        KRATOS_TRY

        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_nodes = GetGeometry().size();

        G_Operator.clear();

        if(dim == 2)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                G_Operator( 0, i*2     ) = DN_DX( i, 0 );
                G_Operator( 1, i*2 + 1 ) = DN_DX( i, 0 );
                G_Operator( 2, i*2     ) = DN_DX( i, 1 );
                G_Operator( 3, i*2 + 1 ) = DN_DX( i, 1 );
            }
        }
        else if(dim == 3)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                G_Operator( 0, i*3     ) = DN_DX( i, 0 );
                G_Operator( 1, i*3 + 1 ) = DN_DX( i, 0 );
                G_Operator( 2, i*3 + 2 ) = DN_DX( i, 0 );
                G_Operator( 3, i*3     ) = DN_DX( i, 1 );
                G_Operator( 4, i*3 + 1 ) = DN_DX( i, 1 );
                G_Operator( 5, i*3 + 2 ) = DN_DX( i, 1 );
                G_Operator( 6, i*3     ) = DN_DX( i, 2 );
                G_Operator( 7, i*3 + 1 ) = DN_DX( i, 2 );
                G_Operator( 8, i*3 + 2 ) = DN_DX( i, 2 );
            }
        }

        KRATOS_CATCH( "" )
    }

    void FiniteStrain::CalculateG( Matrix& G_Operator, const VectorType& N, const Matrix& DN_DX, const Matrix& CurrentDisp ) const
    {
        this->CalculateG( G_Operator, N, DN_DX );
    }

//************************************************************************************
//************************************************************************************

    void FiniteStrain::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo ) const
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

    void FiniteStrain::GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo ) const
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

    void FiniteStrain::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int MatSize = dimension * number_of_nodes;

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

        // for ( unsigned int i = 0; i < number_of_nodes; i++ )
        // {
        //     double temp = LumpFact[i] * TotalMass;

        //     for ( unsigned int j = 0; j < dimension; j++ )
        //     {
        //         unsigned int index = i * dimension + j;
        //         rMassMatrix( index, index ) = temp;
        //     }
        // }

        /// Consistent mass
        // TODO compute the consistent mass in current configuration

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints( mThisIntegrationMethod );
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        Matrix DeltaPosition(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

        GeometryType::JacobiansType J0;
        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );
        double DetJ0;

        //Current displacements
        Matrix CurrentDisp(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        VectorType N(number_of_nodes);

        for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            DetJ0 = MathUtils<double>::Det(J0[PointNumber]);
            noalias(N) = row( Ncontainer, PointNumber );

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight = density * this->GetIntegrationWeight(integration_points[PointNumber].Weight(), N, CurrentDisp);

            //modify integration weight in case of 2D
            if ( dimension == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                for ( unsigned int j = 0; j < number_of_nodes; ++j )
                {
                    for ( unsigned int k = 0; k < dimension; ++k )
                        rMassMatrix(dimension*i + k, dimension*j + k) += N(i) * N(j) * IntToReferenceWeight * DetJ0;
                }
            }
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void FiniteStrain::CalculateDampingMatrix( MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo )
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

    FiniteStrain::IntegrationMethod FiniteStrain::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

//************************************************************************************
//************************************************************************************

    void FiniteStrain::SetValuesOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
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

    void FiniteStrain::SetValuesOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
                    rValues[PointNumber], rCurrentProcessInfo );
        }

    }

//************************************************************************************
//************************************************************************************

    void FiniteStrain::CalculateOnIntegrationPoints( const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int g_size = this->GetGSize(dim);
        const unsigned int f_size = this->GetFSize(dim);

        Matrix F( f_size, f_size );

        Matrix DN_DX( number_of_nodes, dim );

        Matrix GX( g_size, number_of_nodes*dim );

        Vector N( number_of_nodes );

        Matrix CurrentDisp( number_of_nodes, 3 );

        Matrix InvJ0(dim, dim), InvJ(dim, dim);

        double DetJ0, DetJ;

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        Matrix DeltaPosition(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

        //Current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
        if (Id() == DEBUG_ELEMENT_ID)
        {
            std::cout << std::setprecision(16) << std::scientific;
            KRATOS_WATCH(CurrentDisp)
        }
        #endif

        GeometryType::JacobiansType J0;
        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        //calculating actual jacobian
        Matrix ShiftedDisp = DeltaPosition - CurrentDisp;

        #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
        // KRATOS_WATCH(DeltaPosition)
        // KRATOS_WATCH(ShiftedDisp)
        #endif

        // GeometryType::JacobiansType J;
        // J = GetGeometry().Jacobian( J, mThisIntegrationMethod, ShiftedDisp );

        if ( rValues.size() != integration_points.size() )
            rValues.resize( integration_points.size() );

        const bool is_fbar = GetProperties().Has(IS_FBAR) ? GetProperties()[IS_FBAR] : false;

        Matrix FO;
        double DetFO, DetF;
        if (is_fbar)
        {
            CoordinatesArrayType origin;
            GeometryUtility::ComputeOrigin(GetGeometry().GetGeometryFamily(), origin);

            Matrix J0O, JO;
            J0O = GetGeometry().Jacobian( J0O, origin, DeltaPosition );
            JO = GetGeometry().Jacobian( JO, origin, ShiftedDisp );

            Matrix InvJ0O(dim, dim), InvJO(dim, dim);
            double DetJ0O, DetJO;
            Matrix DN_Dx_O( number_of_nodes, dim );

            MathUtils<double>::InvertMatrix( J0O, InvJ0O, DetJ0O );

            // FO
            Vector NO( number_of_nodes );
            NO = GetGeometry().ShapeFunctionsValues(NO, origin);

            Matrix DN_De_O( number_of_nodes, dim );
            DN_De_O = GetGeometry().ShapeFunctionsLocalGradients(DN_De_O, origin);

            Matrix DN_DX_O( number_of_nodes, dim );
            noalias( DN_DX_O ) = prod( DN_De_O, InvJ0O );

            Matrix GOX(g_size, number_of_nodes * dim);
            this->CalculateG( GOX, NO, DN_DX_O );

            FO.resize(f_size, f_size, false);
            this->CalculateF( FO, GOX, CurrentDisp );
            DetFO = MathUtils<double>::Det(FO);
        }

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            //deformation gradient
            MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            noalias( N ) = row( Ncontainer, PointNumber );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );

            /// method 1
            // noalias( F ) = prod( J[PointNumber], InvJ0 );

            /// method 2
            this->CalculateG( GX, N, DN_DX );
            this->CalculateF( F, GX, CurrentDisp );

            /// method 3
            // for (unsigned int i = 0; i < dim; ++i)
            // {
            //     for (unsigned int j = 0; j < dim; ++j)
            //     {
            //         F(i, j) = 0.0;
            //         for (unsigned int n = 0; n < GetGeometry().size(); ++n)
            //             F(i, j) += DN_DX(n, j) * GetGeometry()[n].GetSolutionStepValue(DISPLACEMENT)[i];
            //     }
            //     F(i, i) += 1.0;
            // }

            DetF = MathUtils<double>::Det(F);

            if ( rVariable == CURRENT_DEFORMATION_GRADIENT_DETERMINANT )
            {
                if (is_fbar)
                {
                    #ifdef USE_DETERMINANT_FO_FOR_FBAR
                    rValues[PointNumber] = DetFO;
                    #else
                    if (f_size == 2) // plane strain
                    {
                        rValues[PointNumber] = std::sqrt(DetFO/DetF)*DetF;
                    }
                    else if (f_size == 3) // plane stress, axisymmetric, 3D
                    {
                        rValues[PointNumber] = std::cbrt(DetFO/DetF)*DetF;
                    }
                    #endif
                }
                else
                {
                    rValues[PointNumber] = DetF;
                }
            }
            else
            {
                mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
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

    void FiniteStrain::CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
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

    void FiniteStrain::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int& size = GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
        unsigned int StrainSize;

        if ( GetGeometry().WorkingSpaceDimension() == 2 )
            StrainSize = 3;
        else
            StrainSize = 6;

        if ( rValues.size() != size )
            rValues.resize( size );

        if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
                    PointNumber++ )
            {
                if ( rValues[PointNumber].size() != StrainSize )
                    rValues[PointNumber].resize( StrainSize, false );

                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
            }
        }
        else if ( rVariable == PLASTIC_STRAIN_VECTOR )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != StrainSize )
                    rValues[i].resize( StrainSize );
                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( PLASTIC_STRAIN_VECTOR, rValues[i] );
            }
        }
        else if ( rVariable == STRESSES )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != StrainSize )
                    rValues[i].resize( StrainSize );
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

    void FiniteStrain::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int g_size = this->GetGSize(dim);
        const unsigned int f_size = this->GetFSize(dim);

        Matrix F( f_size, f_size );

        Matrix DN_DX( number_of_nodes, dim );

        Matrix GX( g_size, number_of_nodes*dim );

        Vector N( number_of_nodes );

        Matrix CurrentDisp( number_of_nodes, 3 );

        Matrix InvJ0(dim, dim), InvJ(dim, dim);

        double DetJ0, DetJ;

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        Matrix DeltaPosition(GetGeometry().size(), 3);
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

        //Current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        GeometryType::JacobiansType J0;
        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        //calculating actual jacobian
        Matrix ShiftedDisp = DeltaPosition - CurrentDisp;

        GeometryType::JacobiansType J;
        J = GetGeometry().Jacobian( J, mThisIntegrationMethod, ShiftedDisp );

        if ( rValues.size() != integration_points.size() )
            rValues.resize( integration_points.size() );

        const bool is_fbar = GetProperties().Has(IS_FBAR) ? GetProperties()[IS_FBAR] : false;

        Matrix FO;
        double DetFO, DetF;
        if (is_fbar)
        {
            CoordinatesArrayType origin;
            GeometryUtility::ComputeOrigin(GetGeometry().GetGeometryFamily(), origin);

            Matrix J0O, JO;
            J0O = GetGeometry().Jacobian( J0O, origin, DeltaPosition );
            JO = GetGeometry().Jacobian( JO, origin, ShiftedDisp );

            Matrix InvJ0O(dim, dim), InvJO(dim, dim);
            double DetJ0O, DetJO;
            Matrix DN_Dx_O( number_of_nodes, dim );

            MathUtils<double>::InvertMatrix( J0O, InvJ0O, DetJ0O );

            // FO
            Vector NO( number_of_nodes );
            NO = GetGeometry().ShapeFunctionsValues(NO, origin);

            Matrix DN_De_O( number_of_nodes, dim );
            DN_De_O = GetGeometry().ShapeFunctionsLocalGradients(DN_De_O, origin);

            Matrix DN_DX_O( number_of_nodes, dim );
            noalias( DN_DX_O ) = prod( DN_De_O, InvJ0O );

            Matrix GOX(g_size, number_of_nodes * dim);
            this->CalculateG( GOX, NO, DN_DX_O );

            FO.resize(f_size, f_size, false);
            this->CalculateF( FO, GOX, CurrentDisp );
            DetFO = MathUtils<double>::Det(FO);
        }

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            //deformation gradient
            MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            noalias( N ) = row( Ncontainer, PointNumber );
            // MathUtils<double>::InvertMatrix( J[PointNumber], InvJ, DetJ );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );

            // method 1
            // noalias( F ) = prod( J[PointNumber], InvJ0 );

            // method 2
            this->CalculateG( GX, N, DN_DX );
            this->CalculateF( F, GX, CurrentDisp );

            // method 3
            // for (unsigned int i = 0; i < dim; ++i)
            // {
            //     for (unsigned int j = 0; j < dim; ++j)
            //     {
            //         F(i, j) = 0.0;
            //         for (unsigned int n = 0; n < GetGeometry().size(); ++n)
            //             F(i, j) += DN_DX(n, j) * GetGeometry()[n].GetSolutionStepValue(DISPLACEMENT)[j];
            //     }
            //     F(i, i) += 1.0;
            // }

            if (is_fbar)
            {
                if (f_size == 2) // plane strain
                {
                    DetF = MathUtils<double>::Det(F);
                    F *= std::sqrt(DetFO/DetF);
                }
                else if (f_size == 3) // plane stress, axisymmetric, 3D
                {
                    DetF = MathUtils<double>::Det(F);
                    F *= std::cbrt(DetFO/DetF);
                }
                else
                    KRATOS_ERROR << "Invalid size " << f_size << " of deformation gradient";
            }


            if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            {
                if ( rValues[PointNumber].size1() != f_size || rValues[PointNumber].size2() != f_size )
                    rValues[PointNumber].resize( f_size, f_size, false );

                //strain calculation
                Matrix C( f_size, f_size );
                noalias( C ) = prod( trans( F ), F );

                noalias(rValues[PointNumber]) = C;
            }
            else if ( rVariable == CURRENT_DEFORMATION_GRADIENT )
            {
                if ( rValues[PointNumber].size1() != f_size || rValues[PointNumber].size2() != f_size )
                    rValues[PointNumber].resize( f_size, f_size, false );

                noalias(rValues[PointNumber]) = F;
            }
            else
            {
                mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
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

    void FiniteStrain::GetValuesVector( Vector& values, int Step ) const
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

    void FiniteStrain::GetFirstDerivativesVector( Vector& values, int Step ) const
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

    void FiniteStrain::GetSecondDerivativesVector( Vector& values, int Step ) const
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
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int  FiniteStrain::Check( const ProcessInfo& rCurrentProcessInfo ) const
    {
        KRATOS_TRY

        unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();


        //verify that the constitutive law exists
        if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() );
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

        return 0;

        KRATOS_CATCH( "" );
    }


    void FiniteStrain::save( Serializer& rSerializer ) const
    {
//  std::cout << "Saving the FiniteStrain #" << Id() << std::endl;
        rSerializer.save( "Name", "FiniteStrain" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    }

    void FiniteStrain::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
//  std::cout << "Loading the FiniteStrain #" << Id() << std::endl;
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    }

} // Namespace Kratos

#undef ENABLE_DEBUG_CONSTITUTIVE_LAW
#undef CHECK_DEFORMATION_GRADIENT
#undef USE_DETERMINANT_F0_FOR_FBAR
