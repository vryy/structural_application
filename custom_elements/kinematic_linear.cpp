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
/* *********************************************************
 *
 *   Modified by:         $Author: hbui $
 *   Date:                $Date: 2013-02-22 16:16:48 $
 *   Revision:            $Revision: 1.2 $
 *
 * ***********************************************************/

// System includes


// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"
#include "custom_elements/kinematic_linear.h"
#include "custom_utilities/bathe_recover_stress_utility.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"

//#define ENABLE_DEBUG_CONSTITUTIVE_LAW

//TODO: there is a potential bug at the CalculateRightHandSide, which is used for calculating the reaction. In principal, it should not tell the material to update them-self, however, CalculateRightHandSide indirectly call CalculateMaterialResponse. THis should be fixed, by introducing another abstract layer to update the material in the input parameters for CalculateAll

namespace Kratos
{
    KinematicLinear::KinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry )
        : Element( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
        //THIS IS THE DEFAULT CONSTRUCTOR
    }

    /**
     * A simple kinematic linear 3D element for the solution
     * of the momentum balance in structural mechanics.
     * This element is used for students training at the Ruhr University Bochum.
     * Therefore it may includes comments that are obvious for the
     * experienced user.
     */
    KinematicLinear::KinematicLinear( IndexType NewId,
            GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : Element( NewId, pGeometry, pProperties )
    {
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();//default method
    }

    Element::Pointer KinematicLinear::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new KinematicLinear( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    Element::Pointer KinematicLinear::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new KinematicLinear( NewId, pGeom, pProperties ) );
    }

    KinematicLinear::~KinematicLinear()
    {
    }

    /**
     * Initialization of the element, called at the begin of each simulation.
     * Member variables and the Material law are initialized here
     */
    void KinematicLinear::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY //EXCEPTION HANDLING (see corresponding KRATOS_CATCH("") )

        //dimension of the problem
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        if (rCurrentProcessInfo[RESET_CONFIGURATION] == 0)
        {
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
                    KRATOS_THROW_ERROR(std::logic_error, "KinematicLinear element does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
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
                    KRATOS_THROW_ERROR(std::logic_error, "KinematicLinear element does not support for integration points", GetProperties()[INTEGRATION_ORDER])
            }
            else
                mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

            // // use the default integration rule in the case of finite cell geometry
            // This is not necessary if the INTEGRATION_ORDER is set for finite cell geometry
            // std::string geo_name = typeid(GetGeometry()).name();
            // if ( geo_name.find("FiniteCellGeometry") != std::string::npos )
            //     mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

            // number of integration points used, mThisIntegrationMethod refers to the
            // integration method defined in the constructor
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            // Initialization of the constitutive law vector and
            // declaration, definition and initialization of the material
            // laws at each integration point
            mConstitutiveLawVector.resize( integration_points.size() );
            InitializeMaterial(rCurrentProcessInfo);

            // initialize zero displacement
            mInitialDisp.resize( GetGeometry().size(), dim, false );
            noalias(mInitialDisp) = ZeroMatrix(GetGeometry().size(), dim);

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
    //            std::cout << "Element " << Id() << " mInitialDisp is reinitialized to " << mInitialDisp << std::endl;
            #endif

            // initializing the Jacobian in the reference configuration
            GeometryType::JacobiansType J0;
            Matrix DeltaPosition(GetGeometry().size(), 3);

            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
                noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

            J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

            // calculating the domain size
            double TotalDomainInitialSize = 0.00;
            for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
            {
                //getting informations for integration
                double IntegrationWeight = integration_points[PointNumber].Weight();
                //calculating the total domain size
                TotalDomainInitialSize += MathUtils<double>::Det(J0[PointNumber]) * IntegrationWeight;
            }
    //        KRATOS_WATCH(TotalDomainInitialSize)
            if ( TotalDomainInitialSize < 0.0 )
            {
                if ( TotalDomainInitialSize < -1.0e-10 )
                {
                    std::stringstream ss;
                    ss << "Error on element -> " << this->Id() << std::endl;
                    ss << "Properties " << GetProperties().Id() << ": " << GetProperties() << std::endl;
                    ss << "Domain size can not be less than 0, TotalDomainInitialSize = " << TotalDomainInitialSize;
                    KRATOS_THROW_ERROR( std::logic_error, ss.str(), "" );
                }
                else
                {
                    std::cout << "Warning on element -> " << this->Id();
                    std::cout << ". Domain size is small, TotalDomainInitialSize = " << TotalDomainInitialSize << " < -1e-10";
                    std::cout << ". This element will be deactivated." << std::endl;
                    this->SetValue(ACTIVATION_LEVEL, -1);
                    this->SetValue(IS_INACTIVE, true);
                    this->Set(ACTIVE, false);
                }
            }
            this->SetValue(GEOMETRICAL_DOMAIN_SIZE, TotalDomainInitialSize);
        }
        else if (rCurrentProcessInfo[RESET_CONFIGURATION] == 1)
        {
            //Set Up Initial displacement for StressFreeActivation of Elements
            if (mInitialDisp.size1() != GetGeometry().size() || mInitialDisp.size2() != dim)
                mInitialDisp.resize( GetGeometry().size(), dim, false );

            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
                for ( unsigned int i = 0; i < dim; ++i )
                    mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
    //        std::cout << "Element " << Id() << " mInitialDisp is initialized to " << mInitialDisp << std::endl;
            #endif
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Initialization of the Material law at each integration point
     */
    void KinematicLinear::InitializeMaterial(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();

        int need_shape_function = 0, tmp;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            mConstitutiveLawVector[Point]->SetValue( PARENT_ELEMENT_ID, this->Id(), rCurrentProcessInfo);
            mConstitutiveLawVector[Point]->SetValue( INTEGRATION_POINT_INDEX, Point, rCurrentProcessInfo);
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

                 //verify that the constitutive law has the correct dimension
//                if ( dimension == 2 )
//                {
//                    if ( this->GetProperties().Has( THICKNESS ) == false )
//                        KRATOS_THROW_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() );

//                    if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Getstrain_size() != 3 )
//                        KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) ", this->Id() );
//                }
//                else if(dimension == 3)
//                {
//                    if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Getstrain_size() != 6 )
//                        KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id() );
//                }

                //check constitutive law
                mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
//                if( mConstitutiveLawVector[i]->IsIncremental() )
//                    KRATOS_THROW_ERROR( std::logic_error, "This element does not provide incremental strains!", "" );
//                if( mConstitutiveLawVector[i]->GetStrainMeasure() != ConstitutiveLaw::StrainMeasure_Linear )
//                    KRATOS_THROW_ERROR( std::logic_error, "This element formulated in linear strain measure", "" );
//                if( mConstitutiveLawVector[i]->GetStressMeasure() != ConstitutiveLaw::StressMeasure_PK1 )
//                    KRATOS_THROW_ERROR( std::logic_error, "This element is formulated in PK1 stresses", "" );
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

    void KinematicLinear::ResetConstitutiveLaw()
    {
        KRATOS_TRY

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            ProcessInfo DummyProcessInfo;

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i]->SetValue( PARENT_ELEMENT_ID, this->Id(), DummyProcessInfo);
                mConstitutiveLawVector[i]->SetValue( INTEGRATION_POINT_INDEX, i, DummyProcessInfo);
                mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif

        KRATOS_CATCH( "" )
    }

    /**
     * THIS is the main method here the integration in space (loop over the integration points) is done,
     * the algorithmic tangent and the (inner and outer) load vector is computed
     * @param rLeftHandSideMatrix algorithmic tangent, size (number_of_nodes*dim)*(number_of_nodes*dim)
     * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
     * @param rCurrentProcessInfo
     * @param CalculateStiffnessMatrixFlag true: algorithmic tangent has to be computed
     * @param CalculateResidualVectorFlag true: load vector has to be computed
     */
    void KinematicLinear::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                       VectorType& rRightHandSideVector,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       bool CalculateStiffnessMatrixFlag,
                                       bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

//        std::cout << "start computing element " << Id() << std::endl;

        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);

        //this is the size of the elements stiffness matrix/force vector
        unsigned int mat_size = GetGeometry().size() * dim;

        //Initialize local variables
        Matrix B( strain_size, mat_size );
        Matrix TanC( strain_size, strain_size );
        Vector StrainVector( strain_size );
        Vector StressVector( strain_size );
        Matrix DN_DX( number_of_nodes, dim );
        Matrix CurrentDisp( number_of_nodes, dim );
        Matrix InvJ0(dim, dim);
        double DetJ0;
        double IntToReferenceWeight;

        //constitutive law
        ConstitutiveLaw::Parameters const_params;
        const_params.SetStrainVector(StrainVector);
        if (CalculateResidualVectorFlag)
            const_params.SetStressVector(StressVector);
        if (CalculateStiffnessMatrixFlag)
            const_params.SetConstitutiveMatrix(TanC);
        const_params.SetProcessInfo(rCurrentProcessInfo);
        const_params.SetMaterialProperties(GetProperties());
        const_params.SetElementGeometry(GetGeometry());
        ConstitutiveLaw::StressMeasure stress_measure = ConstitutiveLaw::StressMeasure_Cauchy;

        //resize the LHS=StiffnessMatrix if its size is not correct
        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != mat_size )
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            //resize the RHS=force vector if its size is not correct
            if ( rRightHandSideVector.size() != mat_size )
                rRightHandSideVector.resize( mat_size, false );

            noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
        }

        //int thread_id = omp_get_thread_num();

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        // KRATOS_WATCH(mThisIntegrationMethod)
        // KRATOS_WATCH(integration_points.size())
        // std::cout << "quadrature listing:" << std::endl;
        // for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        // {
        //    std::cout << "integration point " << PointNumber << ": " << integration_points[PointNumber] << std::endl;
        //    Vector shape_values;
        //    GetGeometry().ShapeFunctionsValues(shape_values, integration_points[PointNumber]);
        //    std::cout << "shape values: " << shape_values << std::endl;
        // }

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

//        for(std::size_t i = 0; i < GetGeometry().size(); ++i)
//            std::cout << " " << GetGeometry()[i].Id();
//        std::cout << std::endl;

//        std::vector<std::size_t> equation_ids;
//        this->EquationIdVector(equation_ids, rCurrentProcessInfo);
//        std::cout << "  equation_ids:";
//        for(std::size_t i = 0; i < equation_ids.size(); ++i)
//            std::cout << " " << equation_ids[i];
//        std::cout << std::endl;

        //auxiliary terms
        const Vector& BodyForce = GetProperties()[BODY_FORCE];

        /////////////////////////////////////////////////////////////////////////
        //// Compute the B-dilatational operator
        //// Reference: Thomas Hughes, The Finite Element Method
        /////////////////////////////////////////////////////////////////////////
        Matrix Bdil_bar;
        double TotalDomainInitialSize = this->GetValue(GEOMETRICAL_DOMAIN_SIZE);
        if(GetProperties().Has(IS_BBAR))
        {
            if(GetProperties()[IS_BBAR] == true)
            {
                Bdil_bar.resize(number_of_nodes, dim, false);
                noalias(Bdil_bar) = ZeroMatrix(number_of_nodes, dim);
                for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
                {
                    IntToReferenceWeight = this->GetIntegrationWeight(integration_points, PointNumber, Ncontainer);
                    MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
                    noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );
                    noalias(Bdil_bar) += DN_DX * IntToReferenceWeight * DetJ0;
                }
                Bdil_bar /= TotalDomainInitialSize;
            }
        }
    //    KRATOS_WATCH(Bdil_bar / dim)
    //    KRATOS_WATCH(TotalDomainInitialSize)

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );

            //Initializing B_Operator at the current integration point
            if(GetProperties().Has(IS_BBAR))
            {
                if(GetProperties()[IS_BBAR] == true)
                    CalculateBBaroperator( B, DN_DX, Bdil_bar );
                else
                    CalculateBoperator( B, row(Ncontainer, PointNumber), DN_DX );
            }
            else
                CalculateBoperator( B, row(Ncontainer, PointNumber), DN_DX );

            //calculate strain
            CalculateStrain( B, CurrentDisp, StrainVector );

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
            mConstitutiveLawVector[PointNumber]->SetValue(PARENT_ELEMENT_ID, this->Id(), rCurrentProcessInfo);
            mConstitutiveLawVector[PointNumber]->SetValue(INTEGRATION_POINT_INDEX, PointNumber, rCurrentProcessInfo);
            std::cout << "At element " << Id() << " integration point " << PointNumber << ":" << std::endl;
            std::cout << "B: " << B << std::endl;
            std::cout << "CurrentDisp: " << CurrentDisp << std::endl;
            std::cout << "mInitialDisp: " << mInitialDisp << std::endl;
            std::cout << "StrainVector: " << StrainVector << std::endl;
            #endif

            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(const_params, stress_measure);

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
            if (CalculateResidualVectorFlag)
                std::cout << "StressVector: " << StressVector << std::endl;
            #endif

            //calculating weights for integration on the reference configuration
            IntToReferenceWeight = this->GetIntegrationWeight(integration_points, PointNumber, Ncontainer);

            //modify integration weight in case of 2D
            if ( dim == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                //calculate stiffness matrix
                noalias( rLeftHandSideMatrix ) +=
                    prod( trans( B ), ( IntToReferenceWeight * DetJ0 ) * Matrix( prod( TanC, B ) ) );
            }

            if ( CalculateResidualVectorFlag == true )
            {
                //contribution of external forces
                CalculateAndAdd_ExtForceContribution( row( Ncontainer, PointNumber ), rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight, DetJ0);

                //contribution of gravity (if there is)
                AddBodyForcesToRHS( rRightHandSideVector, row( Ncontainer, PointNumber ), IntToReferenceWeight, DetJ0 );

                //contribution of internal forces
                AddInternalForcesToRHS( rRightHandSideVector, B, StressVector, IntToReferenceWeight, DetJ0 );
            }
        }//loop over integration points

        // modify the right hand side to account for prescribed displacement
        // according to the book of Bazant & Jirasek, this scheme is more stable than the current scheme for prescribing displacement.
        // // However, I have to temporarily disable it to keep the consistency.
        if(CalculateStiffnessMatrixFlag && CalculateResidualVectorFlag)
        {
            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            {
                if(GetGeometry()[node].IsFixed(DISPLACEMENT_X))
                {
                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
                    for( unsigned int i = 0; i < mat_size; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim) * temp;
                }
                if(GetGeometry()[node].IsFixed(DISPLACEMENT_Y))
                {
                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
                    for( unsigned int i = 0; i < mat_size; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 1) * temp;
                }
                if(GetGeometry()[node].IsFixed(DISPLACEMENT_Z) && dim == 3)
                {
                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                    for( unsigned int i = 0; i < mat_size; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 2) * temp;
                }
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //clean the internal data of the geometry
        GetGeometry().Clean();
        #endif

        // if (Id() == 1)
        // {
        //     KRATOS_WATCH(rRightHandSideVector)
        //     KRATOS_WATCH(rLeftHandSideMatrix)
        // }

        KRATOS_CATCH( "" )
    }

    /**
     * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
     * method with CalculateStiffnessMatrixFlag = true and CalculateResidualVectorFlag = false
     * @param rLeftHandSideMatrix (inner and outer) stiffness matrix, size (number_of_nodes*dim x number_of_nodes*dim)
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = false;
        VectorType temp = Vector();

        CalculateAll( rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
    }

    /**
     * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
     * method with CalculateStiffnessMatrixFlag = false and CalculateResidualVectorFlag = true
     * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );

//        //calculation flags
//        bool CalculateStiffnessMatrixFlag = true;
//        bool CalculateResidualVectorFlag = true;
//        MatrixType Ke = Matrix();

//        CalculateAll( Ke, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );

//        const unsigned int Dim = GetGeometry().WorkingSpaceDimension();
//        VectorType u(Ke.size1());
//        std::size_t cnt = 0;
//        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
//        {
//            u(cnt++) = GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT_X);
//            u(cnt++) = GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT_Y);
//            if (Dim == 3)
//                u(cnt++) = GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT_Z);
//        }

//        noalias(rRightHandSideVector) = -prod(Ke, u);
//////        KRATOS_WATCH(rRightHandSideVector)
    }

    /**
     * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
     * method with CalculateStiffnessMatrixFlag = true and CalculateResidualVectorFlag = true
     * @param rLeftHandSideMatrix algorithmic tangent, size (number_of_nodes*dim)*(number_of_nodes*dim)
     * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    /**
     * THIS method is called from the scheme at the start of each solution step
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::InitializeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
        if(mConstitutiveLawVector.size() > 0)
        {
            int need_current_strain_vector = 0, tmp;
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                if (mConstitutiveLawVector[Point]->Has(CURRENT_STRAIN_VECTOR))
                    ++need_current_strain_vector;
            }

            if ( need_current_strain_vector )
            {
                std::vector<Vector> Values;
                this->CalculateOnIntegrationPoints( CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo );
                for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
                {
                    mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
                }
            }
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

    void KinematicLinear::InitializeNonLinearIteration(const ProcessInfo& CurrentProcessInfo)
    {
        //reset all resistant forces at node
        for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
        {
            GetGeometry()[i].GetSolutionStepValue( REACTION_X ) = 0.0;
            GetGeometry()[i].GetSolutionStepValue( REACTION_Y ) = 0.0;
            GetGeometry()[i].GetSolutionStepValue( REACTION_Z ) = 0.0;
        }

//        if( mConstitutiveLawVector[0]->Has( CURRENT_STRAIN_VECTOR ) )
//        {
//            std::vector<Vector> Values;
//            this->CalculateOnIntegrationPoints( CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo );
//            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
//            {
//                mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
//            }
//        }

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

    void KinematicLinear::FinalizeNonLinearIteration(const ProcessInfo& CurrentProcessInfo)
    {
        if(mConstitutiveLawVector.size() > 0)
        {
            if(mConstitutiveLawVector[0]->Has(CURRENT_STRAIN_VECTOR))
            {
                std::vector<Vector> Values;
                this->CalculateOnIntegrationPoints(CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo);
                for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
                {
                    mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
                }
            }
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

    /**
     * THIS method is called from the scheme after each solution step, here the time step
     * start and end point variables can be transferred n --> n+1
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
        if(mConstitutiveLawVector.size() > 0)
        {
            int need_current_strain_vector = 0, tmp;
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                if (mConstitutiveLawVector[Point]->Has(CURRENT_STRAIN_VECTOR))
                    ++need_current_strain_vector;
            }

            if ( need_current_strain_vector )
            {
                std::vector<Vector> Values;
                this->CalculateOnIntegrationPoints( CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo );
                for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
                {
                    mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
                }
            }
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

    void KinematicLinear::MassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_THROW_ERROR(std::logic_error, "Deprecated method", __FUNCTION__)
    }

    void KinematicLinear::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int mat_size = dimension * NumberOfNodes;

        if ( rMassMatrix.size1() != mat_size )
            rMassMatrix.resize( mat_size, mat_size, false );

        noalias( rMassMatrix ) = ZeroMatrix( mat_size, mat_size );

        double density = 0.0;
        if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
            density = GetValue(DENSITY);
        else
            density = GetProperties()[DENSITY];

        /// Lumped mass

        // double TotalDomainInitialSize = this->GetValue(GEOMETRICAL_DOMAIN_SIZE);
        // double TotalMass = density * TotalDomainInitialSize;
        // if ( dimension == 2 ) TotalMass *= GetProperties()[THICKNESS];

        // Vector LumpFact;

        // LumpFact = GetGeometry().LumpingFactors( LumpFact );

        // for ( unsigned int i = 0; i < NumberOfNodes; ++i )
        // {
        //     double temp = LumpFact[i] * TotalMass;

        //     for ( unsigned int j = 0; j < dimension; ++j )
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

        GeometryType::JacobiansType J0;
        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );
        double DetJ0;
        double IntToReferenceWeight;

        for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            DetJ0 = MathUtils<double>::Det(J0[PointNumber]);

            //calculating weights for integration on the reference configuration
            IntToReferenceWeight = this->GetIntegrationWeight(integration_points, PointNumber, Ncontainer);
            IntToReferenceWeight *= density;

            //modify integration weight in case of 2D
            if ( dimension == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

            for ( unsigned int i = 0; i < NumberOfNodes; ++i )
            {
                for ( unsigned int j = 0; j < NumberOfNodes; ++j )
                {
                    for ( unsigned int k = 0; k < dimension; ++k )
                        rMassMatrix(dimension*i + k, dimension*j + k) += Ncontainer(PointNumber, i) * Ncontainer(PointNumber, j)
                            * IntToReferenceWeight * DetJ0;
                }
            }
        }

        // KRATOS_WATCH(Id())
        // KRATOS_WATCH(rMassMatrix)

        KRATOS_CATCH( "" )
    }

    void KinematicLinear::DampMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_THROW_ERROR(std::logic_error, "Deprecated method", __FUNCTION__)
    }

    void KinematicLinear::CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int mat_size = number_of_nodes * dim;

        if ( rDampMatrix.size1() != mat_size )
            rDampMatrix.resize( mat_size, mat_size, false );

        noalias( rDampMatrix ) = ZeroMatrix( mat_size, mat_size );

        //rayleigh damping

        double alpha = 0.0, beta = 0.0;

        if(GetProperties().Has(RAYLEIGH_DAMPING_ALPHA))
        {
            alpha = GetProperties()[RAYLEIGH_DAMPING_ALPHA];
        }

        if(GetProperties().Has(RAYLEIGH_DAMPING_BETA))
        {
            beta = GetProperties()[RAYLEIGH_DAMPING_BETA];
        }

        if (alpha > 0.0)
        {
            CalculateMassMatrix( rDampMatrix, rCurrentProcessInfo );

            rDampMatrix *= alpha;
        }

        if (beta > 0.0)
        {
            Matrix StiffnessMatrix = ZeroMatrix( mat_size, mat_size );

            Vector RHS_Vector = ZeroVector( mat_size );

            CalculateAll( StiffnessMatrix, RHS_Vector, rCurrentProcessInfo, true, false );

            noalias( rDampMatrix ) += beta * StiffnessMatrix;
        }

        // KRATOS_WATCH(Id())
        // KRATOS_WATCH(rDampMatrix)

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    void KinematicLinear::GetValuesVector( Vector& values, int Step ) const
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size )
            values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
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
    void KinematicLinear::GetFirstDerivativesVector( Vector& values, int Step ) const
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size )
            values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
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
    void KinematicLinear::GetSecondDerivativesVector( Vector& values, int Step ) const
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size )
            values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
        }
    }

    /**
     * returns the used integration method
     */
    KinematicLinear::IntegrationMethod KinematicLinear::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

    /**
     * Informations for the builder and solver to assemble the global vectors and matrices.
     * Here a Vector containing the EquationIds of the differnt Dofs is created
     * @param rResult Vector of the EquationIds
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::EquationIdVector( EquationIdVectorType& rResult,
            const ProcessInfo& CurrentProcessInfo ) const
    {
        unsigned int dim = ( GetGeometry().WorkingSpaceDimension() );
        unsigned int mat_size = GetGeometry().size() * dim;

        if ( rResult.size() != mat_size )
            rResult.resize( mat_size, false );

        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            int index = i * dim;
            rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
            if(dim == 3)
                rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
        }
    }

    /**
     * Informations for the builder and solver to assemble the global vectors and matrices.
     * Here a Container containing the pointers rto the DOFs of this element is created
     * @param ElementalDofList Container with of the DOFs associated with the nodes
     *                           of this element
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo&
            CurrentProcessInfo ) const
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            if(dim == 3)
                ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

    /**
     * Adds the Body Forces to the load vector
     * @param R RHS Vector
     * @param N_DISP shape function values at the current integration points
     * @param Weight current integration weight
     * @param detJ current Determinant of the Jacobian
     */
    inline void KinematicLinear::AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, double Weight, double detJ )
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
                R( prim*dim + i ) += N_DISP( prim ) * density * gravity( i ) * detJ * Weight;
            }
        }

        KRATOS_CATCH( "" )
    }

    inline void KinematicLinear::CalculateAndAdd_ExtForceContribution(
            const Vector& N,
            const ProcessInfo& CurrentProcessInfo,
            const Vector& BodyForce,
            VectorType& rRightHandSideVector,
            double weight,
            double detJ
            )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            int index = dimension * i;

            for ( unsigned int j = 0; j < dimension; ++j )
                rRightHandSideVector[index + j] += weight * detJ * N[i] * BodyForce[j];
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Adds the Internal Forces to the load vector
     * @param R RHS Vector
     * @param B_Operator B-Operator at the current integration point
     * @param StressVector current stress vector
     * @param Weight current integration weight
     * @param detJ current Determinant of the Jacobian
     */
//    void KinematicLinear::AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
//    {
//        KRATOS_TRY

//        unsigned int dim = GetGeometry().WorkingSpaceDimension();
//        unsigned int strain_size = dim * (dim + 1) / 2;
//
//        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
//        {
//            for ( unsigned int i = 0; i < dim; ++i )
//            {
//                for ( unsigned int gamma = 0; gamma < strain_size; ++gamma )
//                {
//                    R( prim * dim + i ) += ( -1 ) * ( B_Operator( gamma, prim * dim + i ) * StressVector( gamma ) * detJ * Weight );
//                }
//            }
//        }
//        //         noalias(R) -= detJ * Weight * prod(trans(B_Operator), StressVector);
//
//        KRATOS_CATCH( "" )
//    }

    void KinematicLinear::AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);
        Vector InternalForces(3);

        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                InternalForces(i) = 0.0;
                for ( unsigned int gamma = 0; gamma < strain_size; ++gamma )
                {
                    InternalForces(i) += B_Operator( gamma, prim * dim + i ) * StressVector( gamma ) * detJ * Weight;
                }

                R( prim * dim + i ) -= InternalForces(i);
            }
//            GetGeometry()[prim].GetSolutionStepValue( REACTION ) += InternalForces;
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Adds the Contribution of the current quadrature point to the load vector
     * @param K LHS Matrix
     * @param tan_C 6*6 algorithmic tangent of the materia law (derivation of stresses
     *               regarding strains
     * @param B_Operator B-Operator at the current integration point
     * @param Weight current integration weight
     * @param detJ current Determinant of the Jacobian
     */
    void KinematicLinear::CalculateStiffnesMatrix( Matrix& K, const Matrix& tan_C, const Matrix& B_Operator, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);

        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                for ( unsigned int sec = 0; sec < GetGeometry().size(); ++sec )
                {
                    for ( unsigned int j = 0; j < dim; ++j )
                    {
                        for ( unsigned int alpha = 0; alpha < strain_size; ++alpha )
                            for ( unsigned int beta = 0; beta < strain_size; ++beta )
                                K( prim*dim + i, sec*dim + j ) += B_Operator( alpha, dim * prim + i )
                                    * tan_C( alpha, beta ) * B_Operator( beta, dim * sec + j ) * detJ * Weight;
                    }
                }
            }
        }

//        noalias(K) -= prod( trans(B_Operator), (Weight * detJ) * Matrix(prod(C, B_Operator)) );

        KRATOS_CATCH( "" )
    }

    /**
     * Computes the stress vector and the algorithmic tangent at the current quadrature points
     * regarding the current strain vector
     * @param StressVector output of the method, Cauchy stress vector
     * @param tanC_U algorithmic tangen (derivative Cauchy stress vector/StrainVector) [strain_size*strain_size]
     * @param StrainVector output: current strain vector
     * @param B_Operator current B-operator
     * @param PointNumber number of the current integration point
     * @param CurrentProcessInfo
     */
    void KinematicLinear::CalculateStressAndTangentialStiffness( Vector& StressVector, Matrix& tanC_U,
            Vector& StrainVector, const Matrix& B_Operator, int PointNumber,
            const ProcessInfo& CurrentProcessInfo )
    {
        KRATOS_TRY

        KRATOS_CATCH( "" )
    }

    /**
     * Computes the strain vector
     */
    void KinematicLinear::CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector )
    {
        KRATOS_TRY
        unsigned int Dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(Dim);
        noalias( StrainVector ) = ZeroVector( strain_size );

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            for ( unsigned int item = 0; item < strain_size; ++item )
                for ( unsigned int dim = 0; dim < Dim; ++dim )
                    StrainVector[item] += B( item, Dim * node + dim ) * ( Displacements( node, dim ) - mInitialDisp( node, dim ) );
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Computes the B-Operator at the current quadrature point
     * @param B_Operator current B-operator
     * @param DN_DX shape function values at the current integration point
     */
    void KinematicLinear::CalculateBoperator( Matrix& B_Operator, const Vector& N, const Matrix& DN_DX )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);

        noalias( B_Operator ) = ZeroMatrix( strain_size, number_of_nodes * dim );

        if(dim == 2)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                B_Operator( 0, i*2 ) = DN_DX( i, 0 );
                B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 );
                B_Operator( 2, i*2 ) = DN_DX( i, 1 );
                B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
            }
        }
        else if(dim == 3)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                B_Operator( 0, i*3 ) = DN_DX( i, 0 );
                B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
                B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 );
                B_Operator( 3, i*3 ) = DN_DX( i, 1 );
                B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
                B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
                B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
                B_Operator( 5, i*3 ) = DN_DX( i, 2 );
                B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
            }
        }

        KRATOS_CATCH( "" )
    }

    void KinematicLinear::CalculateBBaroperator( Matrix& B_Operator, const Matrix& DN_DX, const Matrix& Bdil_bar )
    {
        KRATOS_TRY

        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        noalias( B_Operator ) = ZeroMatrix( dim*(dim+1)/2, number_of_nodes * dim );

        if(dim == 2)
        {
            // KRATOS_THROW_ERROR(std::logic_error, "Bbar formulation for 2D is not supported", "")
            // const double R1R3 = (1.0 / 3);
            // double tmp1;
            // double tmp2;

            // for ( unsigned int i = 0; i < number_of_nodes; ++i )
            // {
            //    tmp1 = R1R3 * (Bdil_bar( i, 0 ) - DN_DX( i, 0 ));
            //    tmp2 = R1R3 * (Bdil_bar( i, 1 ) - DN_DX( i, 1 ));

            //    B_Operator( 0, i*2 ) = DN_DX( i, 0 ) + tmp1;
            //    B_Operator( 1, i*2 ) = tmp1;

            //    B_Operator( 0, i*2 + 1) = tmp2;
            //    B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 ) + tmp2;

            //    B_Operator( 2, i*2 ) = DN_DX( i, 1 );
            //    B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
            // }
        }
        else if(dim == 3)
        {
            const double R1R3 = (1.0 / 3);
            double tmp1;
            double tmp2;
            double tmp3;
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                tmp1 = R1R3 * (Bdil_bar( i, 0 ) - DN_DX( i, 0 ));
                tmp2 = R1R3 * (Bdil_bar( i, 1 ) - DN_DX( i, 1 ));
                tmp3 = R1R3 * (Bdil_bar( i, 2 ) - DN_DX( i, 2 ));

                B_Operator( 0, i*3 ) = DN_DX( i, 0 ) + tmp1;
                B_Operator( 1, i*3 ) = tmp1;
                B_Operator( 2, i*3 ) = tmp1;

                B_Operator( 0, i*3 + 1)  = tmp2;
                B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 ) + tmp2;
                B_Operator( 2, i*3 + 1 ) = tmp2;

                B_Operator( 0, i*3 + 2)  = tmp3;
                B_Operator( 1, i*3 + 2 ) = tmp3;
                B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 ) + tmp3;

                B_Operator( 3, i*3 )     = DN_DX( i, 1 );
                B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
                B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
                B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
                B_Operator( 5, i*3 )     = DN_DX( i, 2 );
                B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
            }
        }

        KRATOS_CATCH( "" )
    }

    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        // TODO: This needs to be reviewed (BUI)

        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();;
        unsigned int strain_size = this->GetStrainSize(dim);

        //Initialize local variables
        Matrix B( strain_size, number_of_nodes*dim );
        Matrix TanC( strain_size, strain_size );
        Vector StrainVector( strain_size );
        Vector StressVector( strain_size );
        Matrix DN_DX( number_of_nodes, dim );
        Matrix CurrentDisp( number_of_nodes, dim );
        Matrix InvJ0(dim, dim);
        double DetJ0;

        //constitutive law
        ConstitutiveLaw::Parameters const_params;
        const_params.SetStrainVector(StrainVector);
        const_params.SetStressVector(StressVector);
        const_params.SetConstitutiveMatrix(TanC);
        const_params.SetProcessInfo(rCurrentProcessInfo);
        const_params.SetMaterialProperties(GetProperties());
        const_params.SetElementGeometry(GetGeometry());
        ConstitutiveLaw::StressMeasure stress_measure = ConstitutiveLaw::StressMeasure_Cauchy;

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints( mThisIntegrationMethod );
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        if ( rValues.size() != integration_points.size() )
            rValues.resize( integration_points.size() );

        const GeometryType::ShapeFunctionsGradientsType& DN_De =
            GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        GeometryType::JacobiansType J0;
        Matrix DeltaPosition(GetGeometry().size(), 3);

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
        }

        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        //Current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        //Declaration of the integration weight
        //    double Weight;

        //loop over all integration points
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );

            //Initializing B_Operator at the current integration point
            CalculateBoperator( B, row(Ncontainer, PointNumber), DN_DX );

            if ( rVariable == STRAIN_INTERPOLATION_OPERATOR )
            {
                rValues[PointNumber].resize( B.size1(), B.size2(), false );
                noalias(rValues[PointNumber]) = B;
                continue;
            }

            //calculate strain
            CalculateStrain( B, CurrentDisp, StrainVector );
            //assign the integration weight at the current integration point
            //        Weight = integration_points[PointNumber].Weight();

            //calculate material response
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(const_params, stress_measure);

            if ( rValues[PointNumber].size2() != StrainVector.size() )
                rValues[PointNumber].resize( 1, StrainVector.size(), false );

            if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            {
                for ( unsigned int ii = 0; ii < StrainVector.size(); ++ii )
                    rValues[PointNumber]( 0, ii ) = StrainVector[ii];
            }
            else if ( rVariable == PK2_STRESS_TENSOR )
            {
                for ( unsigned int ii = 0; ii < StrainVector.size(); ++ii )
                    rValues[PointNumber]( 0, ii ) = StressVector[ii];
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry
        GetGeometry().Clean();
        #endif

        KRATOS_CATCH( "" )
    }

    /**
     * Calculate Matrix Variables at each integration point, used for post-processing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the quadrature points, output of the method
     * @param rCurrentProcessInfo
     *
     * Calculate Vector Variables at each integration point, used for post-processing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the quadrature points, output of the method
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int mat_size = dim * number_of_nodes;

        if ( rVariable == RECOVERY_STRESSES )
        {
            /////////////////////////////////////////////////////////////////////////
            //// Calculate recover stresses
            /////////////////////////////////////////////////////////////////////////
            if( GetProperties().Has( STRESS_RECOVERY_TYPE ) == true )
            {
                int Type = GetProperties()[ STRESS_RECOVERY_TYPE ];

                if(Type == 0)
                {
                    // no recovery
                    this->CalculateOnIntegrationPoints( STRESSES, rValues, rCurrentProcessInfo );
                }
                else if(Type == 1)
                {
                    // new recovery method from Bathe

                    int ExpansionLevel = GetProperties()[ NEIGHBOUR_EXPANSION_LEVEL ];

                    BatheRecoverStressUtility StressUtils(ExpansionLevel);
                    StressUtils.CalculateImprovedStressOnIntegrationPoints( *this, rValues, rCurrentProcessInfo );
                }
                else
                    KRATOS_THROW_ERROR(std::logic_error, "The stress recovery type is not supported on element", Id());

            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "The stress recovery method is not defined for element", Id());

        }
        else if ( rVariable == STRAIN || rVariable == CURRENT_STRAIN_VECTOR )
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            // calculate shape function values and local gradients
            Matrix B(strain_size, mat_size);
            Vector StrainVector(strain_size);
            Matrix DN_DX(number_of_nodes, dim);
            Matrix CurrentDisp(number_of_nodes, dim);
            Matrix InvJ0(dim, dim);
            double DetJ0;

            const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

            //initializing the Jacobian in the reference configuration
            GeometryType::JacobiansType J0;
            Matrix DeltaPosition(GetGeometry().size(), 3);

            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            {
                noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
            }

            J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

            // extract current displacements
            for (unsigned int node = 0; node < GetGeometry().size(); ++node)
                noalias(row(CurrentDisp, node)) =
                    GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

            for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            {
                MathUtils<double>::InvertMatrix( J0[i], InvJ0, DetJ0 );
                noalias(DN_DX) = prod(DN_De[i], InvJ0);
                // compute B_Operator at the current integration point
                CalculateBoperator(B, row(Ncontainer, i), DN_DX);

                // compute the strain at integration point
                CalculateStrain(B, CurrentDisp, StrainVector);

                if (rValues[i].size() != strain_size)
                    rValues[i].resize(strain_size);
                noalias(rValues[i]) = StrainVector;
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif
        }
        else if ( rVariable == PRINCIPAL_STRESS )
        {
            std::vector<Vector> stresses;
            this->CalculateOnIntegrationPoints( STRESSES, stresses, rCurrentProcessInfo );

            Matrix sigma(3, 3);
            const double conv = 1.e-8;
            const double zero = 1.0e-12;
            for(std::size_t i = 0; i < stresses.size(); ++i)
            {
                if(rValues[i].size() != 3)
                    rValues[i].resize(3, false);

                SD_MathUtils<double>::StressVectorToTensor(stresses[i], sigma);
                Vector eigenvalues = SD_MathUtils<double>::EigenValuesSym3x3(sigma);

                noalias(rValues[i]) = SD_MathUtils<double>::OrganizeEigenvalues(eigenvalues);
            }
        }
        else if ( rVariable == PRINCIPAL_STRAIN )
        {
            std::vector<Vector> strain;
            this->CalculateOnIntegrationPoints( STRAIN, strain, rCurrentProcessInfo );

            Matrix epsilon(3, 3);
            const double conv = 1.e-8;
            const double zero = 1.0e-12;
            for(std::size_t i = 0; i < strain.size(); ++i)
            {
                if(rValues[i].size() != 3)
                    rValues[i].resize(3, false);

                SD_MathUtils<double>::StrainVectorToTensor(strain[i], epsilon);
                Vector eigenvalues = SD_MathUtils<double>::EigenValuesSym3x3(epsilon);

                noalias(rValues[i]) = SD_MathUtils<double>::OrganizeEigenvalues(eigenvalues);
            }
        }
        else if ( rVariable == NODAL_STRESS_VECTOR )
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

            Vector StressVector(strain_size);

            for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            {
                noalias(StressVector) = ZeroVector(strain_size);
                for (unsigned int j = 0; j < GetGeometry().size(); ++j)
                    noalias(StressVector) += Ncontainer(i, j)*GetGeometry()[j].GetSolutionStepValue(STRESSES);

                if (rValues[i].size() != strain_size)
                    rValues[i].resize(strain_size);
                noalias(rValues[i]) = StressVector;
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif
        }
        else if ( rVariable == STRESSES )
        {
            if ( rValues.size() != mConstitutiveLawVector.size() )
                rValues.resize( mConstitutiveLawVector.size() );

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                if (rValues[i].size() != strain_size)
                    rValues[i].resize(strain_size, false);

                rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
            }
        }
        else
        {
            if ( rValues.size() != mConstitutiveLawVector.size() )
                rValues.resize( mConstitutiveLawVector.size() );

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
            }
        }
    }

    /**
     * Calculate double Variables at each integration point, used for post-processing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the quadrature points, output of the method
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
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

        if( rVariable == INTEGRATION_POINT_GLOBAL || rVariable == INTEGRATION_POINT_GLOBAL_IN_CURRENT_CONFIGURATION )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                rValues[point] = GetGeometry().GlobalCoordinates(rValues[point], integration_points[point]);
            }

            return;
        }

        if( rVariable == INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points =
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod );

            Vector N( GetGeometry().size() );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                GetGeometry().ShapeFunctionsValues( N, integration_points[point] );

                noalias( rValues[point] ) = ZeroVector(3);
                for(std::size_t i = 0 ; i < GetGeometry().size() ; ++i)
                    noalias( rValues[point] ) += N[i] * GetGeometry()[i].GetInitialPosition();
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

    /**
     * Calculate double Variables at each integration point, used for post-processing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the quadrature points, output of the method
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        if( rVariable == K0 )
        {
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
            {
                rValues[Point] = GetValue( K0 );
            }
            return;
        }
        else if( rVariable == STRAIN_ENERGY )
        {
            std::vector<Vector> StrainList(rValues.size());
//            CalculateOnIntegrationPoints(STRAIN, StrainList, rCurrentProcessInfo);
            CalculateOnIntegrationPoints(THREED_STRAIN, StrainList, rCurrentProcessInfo);

            std::vector<Vector> StressList(rValues.size());
//            CalculateOnIntegrationPoints(STRESSES, StressList, rCurrentProcessInfo);
            CalculateOnIntegrationPoints(THREED_STRESSES, StressList, rCurrentProcessInfo);

            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                // calculate strain energy as C = 0.5 * <epsilon, sigma>
                rValues[i] = 0.5 * inner_prod(StrainList[i], StressList[i]);
            }
        }
        else if( rVariable == JACOBIAN_0 )
        {
            // initializing the Jacobian in the reference configuration
            GeometryType::JacobiansType J0;
            Matrix DeltaPosition(GetGeometry().size(), 3);

            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            {
                noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
            }

            J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

            // compute the Jacobian
            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                rValues[i] = MathUtils<double>::Det(J0[i]);
            }
        }
        else if( rVariable == INTEGRATION_WEIGHT )
        {
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                rValues[i] = integration_points[i].Weight();
            }
        }
        else if( rVariable == MATERIAL_DENSITY )
        {
            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                rValues[i] = this->GetValue(MATERIAL_DENSITY);
            }
        }
        else
        {
            //reading integration points and local gradients
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
            {
                rValues[Point] = mConstitutiveLawVector[Point]->GetValue( rVariable, rValues[Point] );
            }
        }
    }

    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if(rVariable == PARENT_ELEMENT_ID)
        {
            rValues.resize(mConstitutiveLawVector.size());
            std::fill(rValues.begin(), rValues.end(), Id());
        }

        if(rVariable == INTEGRATION_POINT_INDEX)
        {
            rValues.resize(mConstitutiveLawVector.size());
            for(unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            {
                rValues[i] = i;
            }
        }
    }

    #ifndef SD_APP_FORWARD_COMPATIBILITY
    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<std::string>& rVariable, std::vector<std::string>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
        {
            rValues[Point] = mConstitutiveLawVector[Point]->GetValue(rVariable, rValues[Point]);
        }
    }
    #endif

    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, std::vector<ConstitutiveLaw::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        if( rVariable == CONSTITUTIVE_LAW )
        {
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
            {
                rValues[Point] = mConstitutiveLawVector[Point];
            }
        }
    }

    /**
     * Set a Matrix Variable from outside
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector of the values on the quadrature points
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::SetValuesOnIntegrationPoints( const Variable<Matrix>& rVariable, const std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            std::stringstream ss;
            ss << "Error at KinematicLinear element " << Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl;
            ss << "rValues.size(): " << rValues.size() << std::endl;
            ss << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
        }
    }

    /**
     * Set a Vector Variable from outside
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector of the values on the quadrature points
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::SetValuesOnIntegrationPoints( const Variable<Vector>& rVariable, const std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            std::stringstream ss;
            ss << "Error at KinematicLinear element " << Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl;
            ss << "rValues.size(): " << rValues.size() << std::endl;
            ss << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }

        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); ++PointNumber )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
        }
    }

    /**
     * Set a Double Variable from outside
     * @param rVariable Global name of the variable to be calculated
     * @param rValue value on the quadrature points
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::SetValuesOnIntegrationPoints( const Variable<double>& rVariable,
            const std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rVariable == K0 )
        {
            SetValue( K0, rValues[0] );
        }
        else
        {
            if ( rValues.size() != mConstitutiveLawVector.size() )
            {
                std::stringstream ss;
                ss << "Error at KinematicLinear element " << Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl;
                ss << "rValues.size(): " << rValues.size() << std::endl;
                ss << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
            }
        }
    }
    /**
     * Set an Int Variable from outside
     * @param rVariable Global name of the variable to be calculated
     * @param rValue value on the quadrature points
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::SetValuesOnIntegrationPoints( const Variable<int>& rVariable,
            const std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            std::stringstream ss;
            ss << "Error at KinematicLinear element " << Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl;
            ss << "rValues.size(): " << rValues.size() << std::endl;
            ss << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
        }
    }

    void KinematicLinear::SetValuesOnIntegrationPoints( const Kratos::Variable<ConstitutiveLaw::Pointer>& rVariable,
            const std::vector< ConstitutiveLaw::Pointer >& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rVariable == CONSTITUTIVE_LAW )
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            if( mConstitutiveLawVector.size() != rValues.size() )
            {
                mConstitutiveLawVector.resize( rValues.size() );
            }

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i] = rValues[i];
                mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif
        }
        else if ( rVariable == CONSTITUTIVE_LAW_NO_INITIALIZE )
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            if( mConstitutiveLawVector.size() != rValues.size() )
            {
                mConstitutiveLawVector.resize( rValues.size() );
            }

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i] = rValues[i];
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif
        }
    }

    int KinematicLinear::Check( const Kratos::ProcessInfo& rCurrentProcessInfo ) const
    {
        KRATOS_TRY

//	    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

        if ( this->Id() < 1 )
        {
            KRATOS_THROW_ERROR( std::logic_error, "Element found with Id 0 or negative, Id() =", Id() );
        }

        //verify that the constitutive law exists
        if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property", this->GetProperties().Id() );
        }

        // verify the strain measure
        auto strain_mearure = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainMeasure();
        if ( strain_mearure != ConstitutiveLaw::StrainMeasure_Infinitesimal
          && strain_mearure != ConstitutiveLaw::StrainMeasure_GreenLagrange )
        {
            std::stringstream ss;
            ss << "The strain measure " << strain_mearure << " is not supported by this element";
            KRATOS_THROW_ERROR( std::logic_error, ss.str(), "" )
        }

        //Verify that the body force is defined
        if ( this->GetProperties().Has( BODY_FORCE ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "BODY_FORCE not provided for property", this->GetProperties().Id() )
        }

        return 0;

        KRATOS_CATCH( "" );
    }

} // Namespace Kratos

#ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
#undef ENABLE_DEBUG_CONSTITUTIVE_LAW
#endif
