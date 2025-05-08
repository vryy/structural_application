/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite BaseType Analysis
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
#include "custom_elements/kinematic_linear.h"

#include "includes/fnv_1a_hash.h"
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"
#include "custom_utilities/bathe_recover_stress_utility.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"


//TODO: there is a potential bug at the CalculateRightHandSide, which is used for calculating the reaction. In principal, it should not tell the material to update them-self, however, CalculateRightHandSide indirectly call CalculateMaterialResponse. THis should be fixed, by introducing another abstract layer to update the material in the input parameters for CalculateAll

namespace Kratos
{

    template<typename TNodeType>
    BaseKinematicLinear<TNodeType>::BaseKinematicLinear( IndexType NewId, typename GeometryType::Pointer pGeometry )
        : BaseType( NewId, pGeometry )
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
    template<typename TNodeType>
    BaseKinematicLinear<TNodeType>::BaseKinematicLinear( IndexType NewId,
            typename GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties )
        : BaseType( NewId, pGeometry, pProperties )
    {
    }

    template<typename TNodeType>
    typename BaseKinematicLinear<TNodeType>::BaseType::Pointer BaseKinematicLinear<TNodeType>::Create( IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties ) const
    {
        return typename BaseType::Pointer( new BaseKinematicLinear( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    template<typename TNodeType>
    typename BaseKinematicLinear<TNodeType>::BaseType::Pointer BaseKinematicLinear<TNodeType>::Create( IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties ) const
    {
        return typename BaseType::Pointer( new BaseKinematicLinear( NewId, pGeom, pProperties ) );
    }

    template<typename TNodeType>
    BaseKinematicLinear<TNodeType>::~BaseKinematicLinear()
    {
    }

    /**
     * Initialization of the element, called at the begin of each simulation.
     * Member variables and the Material law are initialized here
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        //dimension of the problem
        unsigned int dim = this->WorkingSpaceDimension();
        if (rCurrentProcessInfo[RESET_CONFIGURATION] == 0)
        {
            // integration rule
            const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

            // // use the default integration rule in the case of finite cell geometry
            // This is not necessary if the INTEGRATION_ORDER is set for finite cell geometry
            // std::string geo_name = typeid(this->GetGeometry()).name();
            // if ( geo_name.find("FiniteCellGeometry") != std::string::npos )
            //     mThisIntegrationMethod = this->GetGeometry().GetDefaultIntegrationMethod();

            // number of integration points used
            const typename GeometryType::IntegrationPointsArrayType& integration_points = this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

            // Initialization of the constitutive law vector and
            // declaration, definition and initialization of the material
            // laws at each integration point
            mConstitutiveLawVector.resize( integration_points.size() );
            InitializeMaterial(rCurrentProcessInfo);

            // initialize zero displacement
            mInitialDisp.resize( this->GetGeometry().size(), dim, false );
            noalias(mInitialDisp) = ZeroMatrixType(this->GetGeometry().size(), dim);

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
    //            std::cout << "BaseType " << Id() << " mInitialDisp is reinitialized to " << mInitialDisp << std::endl;
            #endif

            // initializing the Jacobian in the reference configuration
            typename GeometryType::JacobiansType J0;
            this->CalculateJacobian( J0, ThisIntegrationMethod );

            // calculating the domain size
            DataType TotalDomainInitialSize = 0.00;
            for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
            {
                //getting informations for integration
                ValueType IntegrationWeight = integration_points[PointNumber].Weight();
                //calculating the total domain size
                TotalDomainInitialSize += MathUtils<DataType>::Det(J0[PointNumber]) * IntegrationWeight;
            }
            // KRATOS_WATCH(integration_points.size())
            // KRATOS_WATCH(TotalDomainInitialSize)
            // // debugging quadrature h8 on cube
            // typedef GeometryType::IntegrationPointType IntegrationPointType;
            // double TotalDomainInitialSize_New = 0.0;
            // MatrixType J0n(dim, dim);
            // std::vector<IntegrationPointType> points;
            // // const double a = 1.00/std::sqrt(3.0);
            // const double a = 0.57735026918962576450914878050195745564760175127012687601860232648397767230293334569371539558574952522520871380513556767665664836499965082627055183736479121617603107730076852735599160670036155830775501;
            // points.push_back(IntegrationPointType(-a, -a, -a, 1.0));
            // points.push_back(IntegrationPointType(a, -a, -a, 1.0));
            // points.push_back(IntegrationPointType(a, a, -a, 1.0));
            // points.push_back(IntegrationPointType(-a, a, -a, 1.0));
            // points.push_back(IntegrationPointType(-a, -a, a, 1.0));
            // points.push_back(IntegrationPointType(a, -a, a, 1.0));
            // points.push_back(IntegrationPointType(a, a, a, 1.0));
            // points.push_back(IntegrationPointType(-a, a, a, 1.0));
            // for ( unsigned int PointNumber = 0; PointNumber < points.size(); ++PointNumber )
            // {
            //     this->GetGeometry().Jacobian(J0n, points[PointNumber]);
            //     TotalDomainInitialSize_New += MathUtils<DataType>::Det(J0n) * points[PointNumber].Weight();
            // }
            // KRATOS_WATCH(integration_points.size())
            // KRATOS_WATCH(TotalDomainInitialSize)
            // const double error = TotalDomainInitialSize - 1.0;
            // const double error_new = TotalDomainInitialSize_New - 1.0;
            // printf("error: %.10e\n", error);
            // printf("error_new: %.10e\n", error_new);
            // KRATOS_ERROR << "stop here";
            // // end of debugging
            // std::stringstream ss;
            // ss << "BaseType " << Id() << " domain size: " << TotalDomainInitialSize;
            // std::cout << ss.str() << std::endl;
            if constexpr (std::is_arithmetic<DataType>::value)
            {
                if ( TotalDomainInitialSize < 0.0 )
                {
                    if ( TotalDomainInitialSize < -1.0e-10 )
                    {
                        KRATOS_ERROR << "Error on element -> " << this->Id() << std::endl
                                     << "Properties " << this->GetProperties().Id() << ": " << this->GetProperties() << std::endl
                                     << "Domain size can not be less than 0, TotalDomainInitialSize = " << TotalDomainInitialSize;
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
            }
            this->SetValue(VARSEL(DataType, GEOMETRICAL_DOMAIN_SIZE), TotalDomainInitialSize);
        }
        else if (rCurrentProcessInfo[RESET_CONFIGURATION] == 1)
        {
            //Set Up Initial displacement for StressFreeActivation of Elements
            if (mInitialDisp.size1() != this->GetGeometry().size() || mInitialDisp.size2() != dim)
                mInitialDisp.resize( this->GetGeometry().size(), dim, false );

            for ( unsigned int node = 0; node < this->GetGeometry().size(); ++node )
                for ( unsigned int i = 0; i < dim; ++i )
                    mInitialDisp( node, i ) = this->GetGeometry()[node].GetSolutionStepValue( VARSEL( DataType, DISPLACEMENT ) )[i];

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
    //        std::cout << "BaseType " << Id() << " mInitialDisp is initialized to " << mInitialDisp << std::endl;
            #endif
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Initialization of the Material law at each integration point
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::InitializeMaterial(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        #pragma omp critical
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
                mConstitutiveLawVector[i] = this->GetProperties()[VARSEL(ConstitutiveLawType, CONSTITUTIVE_LAW)]->Clone();
        }

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
            const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

            #ifdef ENABLE_BEZIER_GEOMETRY
            this->GetGeometry().Initialize(ThisIntegrationMethod);
            #endif

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i]->InitializeMaterial( this->GetProperties(), this->GetGeometry(), row( this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod ), i ) );

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
                mConstitutiveLawVector[i]->Check( this->GetProperties(), this->GetGeometry(), rCurrentProcessInfo );
//                if( mConstitutiveLawVector[i]->IsIncremental() )
//                    KRATOS_THROW_ERROR( std::logic_error, "This element does not provide incremental strains!", "" );
//                if( mConstitutiveLawVector[i]->GetStrainMeasure() != typename ConstitutiveLawType::StrainMeasure_Linear )
//                    KRATOS_THROW_ERROR( std::logic_error, "This element formulated in linear strain measure", "" );
//                if( mConstitutiveLawVector[i]->GetStressMeasure() != typename ConstitutiveLawType::StressMeasure_PK1 )
//                    KRATOS_THROW_ERROR( std::logic_error, "This element is formulated in PK1 stresses", "" );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            this->GetGeometry().Clean();
            #endif
        }
        else
        {
            Vector dummy;
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                mConstitutiveLawVector[i]->InitializeMaterial( this->GetProperties(), this->GetGeometry(), dummy );
            }
        }

        KRATOS_CATCH( "" )
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

        #ifdef ENABLE_BEZIER_GEOMETRY
        this->GetGeometry().Initialize(ThisIntegrationMethod);
        #endif

        // ProcessInfo DummyProcessInfo;

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            // mConstitutiveLawVector[i]->SetValue( PARENT_ELEMENT_ID, this->Id(), DummyProcessInfo);
            // mConstitutiveLawVector[i]->SetValue( INTEGRATION_POINT_INDEX, i, DummyProcessInfo);
            mConstitutiveLawVector[i]->ResetMaterial( this->GetProperties(), this->GetGeometry(), row( this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod ), i ) );
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        this->GetGeometry().Clean();
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
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                       VectorType& rRightHandSideVector,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       bool CalculateStiffnessMatrixFlag,
                                       bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        // std::cout << "start computing element " << Id() << std::endl;

        unsigned int number_of_nodes = this->GetGeometry().size();
        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);

        //this is the size of the elements stiffness matrix/force vector
        unsigned int mat_size = this->GetGeometry().size() * dim;

        //Initialize local variables
        MatrixType B( strain_size, mat_size );
        MatrixType TanC( strain_size, strain_size );
        VectorType StrainVector( strain_size );
        VectorType StressVector( strain_size );
        MatrixType DN_DX( number_of_nodes, dim );
        Vector N( number_of_nodes );
        MatrixType CurrentDisp( number_of_nodes, dim );
        MatrixType InvJ0(dim, dim);
        DataType DetJ0;
        DataType IntToReferenceWeight;
        const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

        //constitutive law
        typename ConstitutiveLawType::Parameters const_params;
        const_params.SetStrainVector(StrainVector);
        if (CalculateResidualVectorFlag)
            const_params.SetStressVector(StressVector);
        if (CalculateStiffnessMatrixFlag)
            const_params.SetConstitutiveMatrix(TanC);
        const_params.SetProcessInfo(rCurrentProcessInfo);
        const_params.SetMaterialProperties(this->GetProperties());
        const_params.SetElementGeometry(this->GetGeometry());
        typename ConstitutiveLawType::StressMeasure stress_measure = ConstitutiveLawType::StressMeasure_Cauchy;

        //resize the LHS=StiffnessMatrix if its size is not correct
        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != mat_size )
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrixType( mat_size, mat_size ); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            //resize the RHS=force vector if its size is not correct
            if ( rRightHandSideVector.size() != mat_size )
                rRightHandSideVector.resize( mat_size, false );

            noalias( rRightHandSideVector ) = ZeroVectorType( mat_size ); //resetting RHS
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        this->GetGeometry().Initialize(ThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const typename GeometryType::IntegrationPointsArrayType& integration_points = this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

        // KRATOS_WATCH(integration_points.size())
        // std::cout << "quadrature listing:" << std::endl;
        // for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        // {
        //    std::cout << "integration point " << PointNumber << ": " << integration_points[PointNumber] << std::endl;
        //    VectorType shape_values;
        //    this->GetGeometry().ShapeFunctionsValues(shape_values, integration_points[PointNumber]);
        //    std::cout << "shape values: " << shape_values << std::endl;
        // }

        const typename GeometryType::ShapeFunctionsGradientsType& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients( ThisIntegrationMethod );

        const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        typename GeometryType::JacobiansType J0;
        this->CalculateJacobian( J0, ThisIntegrationMethod );

        //Current displacements
        for ( unsigned int node = 0; node < this->GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = subrange(this->GetGeometry()[node].GetSolutionStepValue( VARSEL( DataType, DISPLACEMENT ) ), 0, dim);

        // if (Id() == 749)
        // {
        //     KRATOS_WATCH(CurrentDisp)
        // }

//        for(std::size_t i = 0; i < this->GetGeometry().size(); ++i)
//            std::cout << " " << this->GetGeometry()[i].Id();
//        std::cout << std::endl;

//        std::vector<std::size_t> equation_ids;
//        this->EquationIdVector(equation_ids, rCurrentProcessInfo);
//        std::cout << "  equation_ids:";
//        for(std::size_t i = 0; i < equation_ids.size(); ++i)
//            std::cout << " " << equation_ids[i];
//        std::cout << std::endl;

        //auxiliary terms
        const VectorType& BodyForce = this->GetProperties()[BODY_FORCE];

        /////////////////////////////////////////////////////////////////////////
        //// Compute the B-dilatational operator
        //// Reference: Thomas Hughes, The Finite BaseType Method
        /////////////////////////////////////////////////////////////////////////
        bool is_bbar = false;
        if(this->GetProperties().Has(IS_BBAR))
            is_bbar = this->GetProperties()[IS_BBAR];
        MatrixType Bdil_bar;
        DataType TotalDomainInitialSize = this->GetValue( VARSEL(DataType, GEOMETRICAL_DOMAIN_SIZE) );
        if(is_bbar)
        {
            Bdil_bar.resize(number_of_nodes, dim, false);
            noalias(Bdil_bar) = ZeroMatrixType(number_of_nodes, dim);
            for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
            {
                noalias( N ) = row(Ncontainer, PointNumber);
                IntToReferenceWeight = this->GetIntegrationWeight(integration_points[PointNumber].Weight(), N);
                if ( dim == 2 ) IntToReferenceWeight *= this->GetProperties()[THICKNESS];
                MathUtils<DataType>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
                noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );
                noalias(Bdil_bar) += DN_DX * IntToReferenceWeight * DetJ0;
            }
            Bdil_bar /= TotalDomainInitialSize;
        }
    //    KRATOS_WATCH(Bdil_bar / dim)
    //    KRATOS_WATCH(TotalDomainInitialSize)

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            MathUtils<DataType>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );
            noalias( N ) = row(Ncontainer, PointNumber);

            //Initializing B_Operator at the current integration point
            if(is_bbar)
                CalculateBBaroperator( B, DN_DX, Bdil_bar );
            else
                CalculateBoperator( B, N, DN_DX );

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
            IntToReferenceWeight = this->GetIntegrationWeight(integration_points[PointNumber].Weight(), N);

            //modify integration weight in case of 2D
            if ( dim == 2 ) IntToReferenceWeight *= this->GetProperties()[THICKNESS];

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                //calculate stiffness matrix
                noalias( rLeftHandSideMatrix ) +=
                    prod( trans( B ), ( IntToReferenceWeight * DetJ0 ) * MatrixType( prod( TanC, B ) ) );
            }

            if ( CalculateResidualVectorFlag == true )
            {
                //contribution of external forces
                CalculateAndAdd_ExtForceContribution( N, rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight, DetJ0);

                //contribution of gravity (if there is)
                AddBodyForcesToRHS( rRightHandSideVector, N, IntToReferenceWeight, DetJ0 );

                //contribution of internal forces
                AddInternalForcesToRHS( rRightHandSideVector, B, StressVector, IntToReferenceWeight, DetJ0 );
            }
        }//loop over integration points

        // if (Id() == DEBUG_ELEMENT_ID)
        // {
        //     Vector rExternalForces = rRightHandSideVector;
        //     rExternalForces.clear();
        //     for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        //     {
        //         MathUtils<DataType>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
        //         noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );

        //         //calculating weights for integration on the reference configuration
        //         IntToReferenceWeight = this->GetIntegrationWeight(integration_points, PointNumber, Ncontainer);

        //         //modify integration weight in case of 2D
        //         if ( dim == 2 ) IntToReferenceWeight *= this->GetProperties()[THICKNESS];

        //         CalculateAndAdd_ExtForceContribution( N, rCurrentProcessInfo, BodyForce, rExternalForces, IntToReferenceWeight, DetJ0);
        //         AddBodyForcesToRHS( rExternalForces, N, IntToReferenceWeight, DetJ0 );
        //     }

        //     KRATOS_WATCH(rExternalForces)
            // KRATOS_WATCH(rRightHandSideVector)
            // KRATOS_WATCH(rLeftHandSideMatrix)
        // }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //clean the internal data of the geometry
        this->GetGeometry().Clean();
        #endif

        KRATOS_CATCH( "" )
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::ApplyPrescribedDofs(const MatrixType& LHS_Contribution, VectorType& RHS_Contribution, const ProcessInfo& CurrentProcessInfo) const
    {
        // modify the right hand side to account for prescribed displacement
        // according to the book of Bazant & Jirasek, this scheme is more stable than the total displacement scheme for prescribing displacement.
        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int mat_size = dim * this->GetGeometry().size();
        for ( unsigned int node = 0; node < this->GetGeometry().size(); ++node )
        {
            if(this->GetGeometry()[node].IsFixed(VARSELC(DataType, DISPLACEMENT, X)))
            {
                DataType temp = this->GetGeometry()[node].GetSolutionStepValue(VARSELC(DataType, PRESCRIBED_DELTA_DISPLACEMENT, X));
                if (temp != 0.0)
                    for( unsigned int i = 0; i < mat_size; ++i )
                        RHS_Contribution[i] -= LHS_Contribution(i, node * dim) * temp;
            }

            if(this->GetGeometry()[node].IsFixed(VARSELC(DataType, DISPLACEMENT, Y)))
            {
                DataType temp = this->GetGeometry()[node].GetSolutionStepValue(VARSELC(DataType, PRESCRIBED_DELTA_DISPLACEMENT, Y));
                if (temp != 0.0)
                    for( unsigned int i = 0; i < mat_size; ++i )
                        RHS_Contribution[i] -= LHS_Contribution(i, node * dim + 1) * temp;
            }

            if (dim > 2)
            {
                if(this->GetGeometry()[node].IsFixed(VARSELC(DataType, DISPLACEMENT, Z)))
                {
                    DataType temp = this->GetGeometry()[node].GetSolutionStepValue(VARSELC(DataType, PRESCRIBED_DELTA_DISPLACEMENT, Z));
                    if (temp != 0.0)
                        for( unsigned int i = 0; i < mat_size; ++i )
                            RHS_Contribution[i] -= LHS_Contribution(i, node * dim + 2) * temp;
                }
            }
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::ComputePrescribedForces(const MatrixType& LHS_Contribution, VectorType& Force, const ProcessInfo& CurrentProcessInfo) const
    {
        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int mat_size = dim * this->GetGeometry().size();
        for ( unsigned int node = 0; node < this->GetGeometry().size(); ++node )
        {
            if(this->GetGeometry()[node].IsFixed(VARSELC(DataType, DISPLACEMENT, X)))
            {
                DataType temp = this->GetGeometry()[node].GetSolutionStepValue(VARSELC(DataType, PRESCRIBED_DISPLACEMENT, X));
                if (temp != 0.0)
                    for( unsigned int i = 0; i < mat_size; ++i )
                        Force[i] -= LHS_Contribution(i, node * dim) * temp;
            }

            if(this->GetGeometry()[node].IsFixed(VARSELC(DataType, DISPLACEMENT, Y)))
            {
                DataType temp = this->GetGeometry()[node].GetSolutionStepValue(VARSELC(DataType, PRESCRIBED_DISPLACEMENT, Y));
                if (temp != 0.0)
                    for( unsigned int i = 0; i < mat_size; ++i )
                        Force[i] -= LHS_Contribution(i, node * dim + 1) * temp;
            }

            if (dim > 2)
            {
                if(this->GetGeometry()[node].IsFixed(VARSELC(DataType, DISPLACEMENT, Z)))
                {
                    DataType temp = this->GetGeometry()[node].GetSolutionStepValue(VARSELC(DataType, PRESCRIBED_DISPLACEMENT, Z));
                    if (temp != 0.0)
                        for( unsigned int i = 0; i < mat_size; ++i )
                            Force[i] -= LHS_Contribution(i, node * dim + 2) * temp;
                }
            }
        }
    }

    /**
     * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
     * method with CalculateStiffnessMatrixFlag = true and CalculateResidualVectorFlag = false
     * @param rLeftHandSideMatrix (inner and outer) stiffness matrix, size (number_of_nodes*dim x number_of_nodes*dim)
     * @param rCurrentProcessInfo
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = false;
        VectorType temp = VectorType();

        CalculateAll( rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
    }

    /**
     * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
     * method with CalculateStiffnessMatrixFlag = false and CalculateResidualVectorFlag = true
     * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
     * @param rCurrentProcessInfo
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = MatrixType();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );

//        //calculation flags
//        bool CalculateStiffnessMatrixFlag = true;
//        bool CalculateResidualVectorFlag = true;
//        MatrixType Ke = MatrixType();

//        CalculateAll( Ke, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );

//        const unsigned int Dim = this->WorkingSpaceDimension();
//        VectorType u(Ke.size1());
//        std::size_t cnt = 0;
//        for ( unsigned int node = 0; node < this->GetGeometry().size(); ++node )
//        {
//            u(cnt++) = this->GetGeometry()[node].GetSolutionStepValue(VARSELC( DataType, DISPLACEMENT, X ));
//            u(cnt++) = this->GetGeometry()[node].GetSolutionStepValue(VARSELC( DataType, DISPLACEMENT, Y ));
//            if (Dim == 3)
//                u(cnt++) = this->GetGeometry()[node].GetSolutionStepValue(VARSELC( DataType, DISPLACEMENT, Z ));
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
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

        #ifdef COMPUTE_NUMERICAL_DERIVATIVE
        if (Id() == DEBUG_ELEMENT_ID)
        {
            KRATOS_WATCH(rLeftHandSideMatrix)

            const ValueType epsilon = this->GetProperties()[LOCAL_ERROR_TOLERANCE];
            MatrixType NumericalStiffness(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2());
            CalculateNumericalStiffness(NumericalStiffness, rCurrentProcessInfo, epsilon);
            KRATOS_WATCH(NumericalStiffness)

            const MatrixType DiffStiff = rLeftHandSideMatrix - NumericalStiffness;
            const DataType norm_diff_stiff = norm_frobenius(DiffStiff);
            const DataType norm_stiff = norm_frobenius(rLeftHandSideMatrix);
            KRATOS_WATCH(norm_diff_stiff / norm_stiff)
        }
        #endif
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateNumericalStiffness(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo, const ValueType epsilon)
    {
        unsigned int dim_disp = ( this->WorkingSpaceDimension() );
        unsigned int mat_size = this->GetGeometry().size() * dim_disp;

        VectorType RefRhs(mat_size), NewRhs(mat_size);
        MatrixType dummy;
        CalculateAll( dummy, RefRhs, rCurrentProcessInfo, false, true );

        for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
        {
            const auto& disp = this->GetGeometry()[i].GetSolutionStepValue(VARSEL( DataType, DISPLACEMENT ));
            auto new_disp = disp;

            for (unsigned int j = 0; j < dim_disp; ++j)
            {
                noalias(new_disp) = disp;
                new_disp[j] += epsilon;
                noalias(this->GetGeometry()[i].GetSolutionStepValue(VARSEL( DataType, DISPLACEMENT ))) = new_disp;

                CalculateAll( dummy, NewRhs, rCurrentProcessInfo, false, true );

                column(rLeftHandSideMatrix, i*dim_disp + j) = (RefRhs - NewRhs) / epsilon;

                noalias(this->GetGeometry()[i].GetSolutionStepValue(VARSEL( DataType, DISPLACEMENT ))) = disp;
            }
        }

        // calculate one more time to restore the element state
        CalculateAll( dummy, NewRhs, rCurrentProcessInfo, false, true );
    }

    /**
     * THIS method is called from the scheme at the start of each solution step
     * @param rCurrentProcessInfo
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::InitializeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
        // check if the constitutive law need current strain
        bool need_current_strain_vector = false;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            if (mConstitutiveLawVector[Point]->Has(VARSEL(DataType, CURRENT_STRAIN_VECTOR)))
            {
                need_current_strain_vector = true;
                break;
            }
        }

        if ( need_current_strain_vector )
        {
            std::vector<VectorType> Values;
            this->CalculateOnIntegrationPoints( VARSEL(DataType, CURRENT_STRAIN_VECTOR), Values, CurrentProcessInfo );
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->SetValue( VARSEL(DataType, CURRENT_STRAIN_VECTOR), Values[Point], CurrentProcessInfo );
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
            const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            this->GetGeometry().Initialize(ThisIntegrationMethod);
            #endif

            const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->InitializeSolutionStep( this->GetProperties(), this->GetGeometry(), row(Ncontainer, Point), CurrentProcessInfo );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            this->GetGeometry().Clean();
            #endif
        }
        else
        {
            Vector dummy;
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->InitializeSolutionStep( this->GetProperties(), this->GetGeometry(), dummy, CurrentProcessInfo );
            }
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::InitializeNonLinearIteration(const ProcessInfo& CurrentProcessInfo)
    {
        //reset all resistant forces at node
        for ( unsigned int i = 0; i < this->GetGeometry().size(); ++i )
        {
            this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, REACTION, X ) ) = 0.0;
            this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, REACTION, Y ) ) = 0.0;
            this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, REACTION, Z ) ) = 0.0;
        }

//        if( mConstitutiveLawVector[0]->Has( CURRENT_STRAIN_VECTOR ) )
//        {
//            std::vector<VectorType> Values;
//            this->CalculateOnIntegrationPoints( CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo );
//            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
//            {
//                mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
//            }
//        }

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
            const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            this->GetGeometry().Initialize(ThisIntegrationMethod);
            #endif

            const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->InitializeNonLinearIteration( this->GetProperties(), this->GetGeometry(), row(Ncontainer, Point), CurrentProcessInfo );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            this->GetGeometry().Clean();
            #endif
        }
        else
        {
            Vector dummy;
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->InitializeNonLinearIteration( this->GetProperties(), this->GetGeometry(), dummy, CurrentProcessInfo );
            }
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::FinalizeNonLinearIteration(const ProcessInfo& CurrentProcessInfo)
    {
        // check if the constitutive law need current strain
        bool need_current_strain_vector = false;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            if (mConstitutiveLawVector[Point]->Has(VARSEL(DataType, CURRENT_STRAIN_VECTOR)))
            {
                need_current_strain_vector = true;
                break;
            }
        }

        if ( need_current_strain_vector )
        {
            std::vector<VectorType> Values;
            this->CalculateOnIntegrationPoints( VARSEL(DataType, CURRENT_STRAIN_VECTOR), Values, CurrentProcessInfo );
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->SetValue( VARSEL(DataType, CURRENT_STRAIN_VECTOR), Values[Point], CurrentProcessInfo );
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
            const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            this->GetGeometry().Initialize(ThisIntegrationMethod);
            #endif

            const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->FinalizeNonLinearIteration( this->GetProperties(), this->GetGeometry(), row(Ncontainer, Point), CurrentProcessInfo );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            this->GetGeometry().Clean();
            #endif
        }
        else
        {
            Vector dummy;
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->FinalizeNonLinearIteration( this->GetProperties(), this->GetGeometry(), dummy, CurrentProcessInfo );
            }
        }
    }

    /**
     * THIS method is called from the scheme after each solution step, here the time step
     * start and end point variables can be transferred n --> n+1
     * @param rCurrentProcessInfo
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
        // check if the constitutive law need current strain
        bool need_current_strain_vector = false;
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            if (mConstitutiveLawVector[Point]->Has(VARSEL(DataType, CURRENT_STRAIN_VECTOR)))
            {
                need_current_strain_vector = true;
                break;
            }
        }

        if ( need_current_strain_vector )
        {
            std::vector<VectorType> Values;
            this->CalculateOnIntegrationPoints( VARSEL(DataType, CURRENT_STRAIN_VECTOR), Values, CurrentProcessInfo );
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->SetValue( VARSEL(DataType, CURRENT_STRAIN_VECTOR), Values[Point], CurrentProcessInfo );
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
            const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            this->GetGeometry().Initialize(ThisIntegrationMethod);
            #endif

            const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->FinalizeSolutionStep( this->GetProperties(), this->GetGeometry(), row(Ncontainer, Point), CurrentProcessInfo );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            this->GetGeometry().Clean();
            #endif
        }
        else
        {
            Vector dummy;
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->FinalizeSolutionStep( this->GetProperties(), this->GetGeometry(), dummy, CurrentProcessInfo );
            }
        }
    }

    //************************************************************************************
    //************************************************************************************
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::MassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_ERROR << "Deprecated method";
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = this->WorkingSpaceDimension();
        unsigned int number_of_nodes = this->GetGeometry().size();
        unsigned int mat_size = dimension * number_of_nodes;

        if ( rMassMatrix.size1() != mat_size )
            rMassMatrix.resize( mat_size, mat_size, false );

        noalias( rMassMatrix ) = ZeroMatrixType( mat_size, mat_size );

        ValueType density = 0.0;
        if( this->GetValue(USE_DISTRIBUTED_PROPERTIES) )
            density = this->GetValue( DENSITY );
        else
            density = this->GetProperties()[ DENSITY ];

        /// Lumped mass

        // double TotalDomainInitialSize = this->GetValue(GEOMETRICAL_DOMAIN_SIZE);
        // double TotalMass = density * TotalDomainInitialSize;
        // if ( dimension == 2 ) TotalMass *= this->GetProperties()[THICKNESS];

        // VectorType LumpFact;

        // LumpFact = this->GetGeometry().LumpingFactors( LumpFact );

        // for ( unsigned int i = 0; i < number_of_nodes; ++i )
        // {
        //     double temp = LumpFact[i] * TotalMass;

        //     for ( unsigned int j = 0; j < dimension; ++j )
        //     {
        //         unsigned int index = i * dimension + j;
        //         rMassMatrix( index, index ) = temp;
        //     }
        // }

        /// Consistent mass

        const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        this->GetGeometry().Initialize(ThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const typename GeometryType::IntegrationPointsArrayType& integration_points =
            this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

        const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        typename GeometryType::JacobiansType J0;
        this->CalculateJacobian( J0, ThisIntegrationMethod );

        Vector N(number_of_nodes);
        DataType DetJ0;
        DataType IntToReferenceWeight;

        for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            DetJ0 = MathUtils<DataType>::Det(J0[PointNumber]);
            noalias( N ) = row( Ncontainer, PointNumber );

            //calculating weights for integration on the reference configuration
            IntToReferenceWeight = this->GetIntegrationWeight(integration_points[PointNumber].Weight(), N) * density;

            //modify integration weight in case of 2D
            if ( dimension == 2 ) IntToReferenceWeight *= this->GetProperties()[THICKNESS];

            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                for ( unsigned int j = 0; j < number_of_nodes; ++j )
                {
                    for ( unsigned int k = 0; k < dimension; ++k )
                        rMassMatrix(dimension*i + k, dimension*j + k) += N(i) * N(j) * IntToReferenceWeight * DetJ0;
                }
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //clean the internal data of the geometry
        this->GetGeometry().Clean();
        #endif

        KRATOS_CATCH( "" )
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::DampMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_ERROR << "Deprecated method";
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        unsigned int number_of_nodes = this->GetGeometry().size();
        unsigned int dim = this->WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int mat_size = number_of_nodes * dim;

        if ( rDampMatrix.size1() != mat_size )
            rDampMatrix.resize( mat_size, mat_size, false );

        noalias( rDampMatrix ) = ZeroMatrixType( mat_size, mat_size );

        // Rayleigh damping

        ValueType alpha = 0.0, beta = 0.0;

        if(this->GetProperties().Has(RAYLEIGH_DAMPING_ALPHA))
        {
            alpha = this->GetProperties()[RAYLEIGH_DAMPING_ALPHA];
        }

        if(this->GetProperties().Has(RAYLEIGH_DAMPING_BETA))
        {
            beta = this->GetProperties()[RAYLEIGH_DAMPING_BETA];
        }

        if (alpha > 0.0)
        {
            CalculateMassMatrix( rDampMatrix, rCurrentProcessInfo );

            rDampMatrix *= alpha;
        }

        if (beta > 0.0)
        {
            MatrixType StiffnessMatrix = ZeroMatrixType( mat_size, mat_size );

            VectorType RHS_Vector = ZeroVectorType( mat_size );

            CalculateAll( StiffnessMatrix, RHS_Vector, rCurrentProcessInfo, true, false );

            noalias( rDampMatrix ) += beta * StiffnessMatrix;
        }

        // KRATOS_WATCH(Id())
        // KRATOS_WATCH(rDampMatrix)

        KRATOS_CATCH( "" )
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::AddInertiaForces(VectorType& rRightHandSideVector, DataType coeff, const ProcessInfo& rCurrentProcessInfo)
    {
        VectorType Acceleration;
        this->GetSecondDerivativesVector(Acceleration, 0);

        MatrixType MassMatrix;
        this->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

        if (coeff == 1.0)
            noalias(rRightHandSideVector) -= prod(MassMatrix, Acceleration);
        else
            noalias(rRightHandSideVector) -= coeff * prod(MassMatrix, Acceleration);
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::AddDampingForces(VectorType& rRightHandSideVector, DataType coeff, const ProcessInfo& rCurrentProcessInfo)
    {
        VectorType Velocity;
        this->GetFirstDerivativesVector(Velocity, 0);

        MatrixType DampMatrix;
        this->CalculateDampingMatrix(DampMatrix, rCurrentProcessInfo);

        if (coeff == 1.0)
            noalias(rRightHandSideVector) -= prod(DampMatrix, Velocity);
        else
            noalias(rRightHandSideVector) -= coeff * prod(DampMatrix, Velocity);
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateLocalAccelerationContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        VectorType Acceleration;
        this->GetSecondDerivativesVector(Acceleration, 0);

        if (rCurrentProcessInfo[TIME_INTEGRATION_SCHEME] == 0) // default behavior for linear structural dynamics
        {
            MatrixType& rMassMatrix = rLeftHandSideMatrix;
            this->CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);

            noalias(rRightHandSideVector) -= prod(rMassMatrix, Acceleration);
        }
        else if (rCurrentProcessInfo[TIME_INTEGRATION_SCHEME] == FNV1a32Hash::CalculateHash("ResidualBasedNewmarkScheme<1>"))
        {
            const ValueType alpha_m = rCurrentProcessInfo[NEWMARK_ALPHAM];
            const ValueType beta = rCurrentProcessInfo[NEWMARK_BETA];
            const ValueType Dt = rCurrentProcessInfo[DELTA_TIME];

            ///

            MatrixType MassMatrix;
            this->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

            noalias(rRightHandSideVector) -= prod(MassMatrix, Acceleration);

            const ValueType aux = (1 - alpha_m) / (beta*pow(Dt, 2));
            noalias(rLeftHandSideMatrix) += aux * MassMatrix;
        }
        else
        {
            KRATOS_ERROR << "BaseKinematicLinear::" << __FUNCTION__ << " is not yet implemented for time integration scheme "
                         << rCurrentProcessInfo[TIME_INTEGRATION_SCHEME];
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        VectorType Velocity;
        this->GetFirstDerivativesVector(Velocity, 0);

        if (rCurrentProcessInfo[TIME_INTEGRATION_SCHEME] == 0) // default behavior for linear structural dynamics
        {
            MatrixType& rDampMatrix = rLeftHandSideMatrix;
            this->CalculateDampingMatrix(rDampMatrix, rCurrentProcessInfo);

            noalias(rRightHandSideVector) -= prod(rDampMatrix, Velocity);
        }
        else if (rCurrentProcessInfo[TIME_INTEGRATION_SCHEME] == FNV1a32Hash::CalculateHash("ResidualBasedNewmarkScheme<1>"))
        {
            const ValueType alpha_f = rCurrentProcessInfo[NEWMARK_ALPHAF];
            const ValueType gamma = rCurrentProcessInfo[NEWMARK_GAMMA];
            const ValueType beta = rCurrentProcessInfo[NEWMARK_BETA];
            const ValueType Dt = rCurrentProcessInfo[DELTA_TIME];

            ///

            MatrixType DampMatrix;
            this->CalculateDampingMatrix(DampMatrix, rCurrentProcessInfo);

            noalias(rRightHandSideVector) -= prod(DampMatrix, Velocity);

            const ValueType aux = (1 - alpha_f) * gamma/(beta*Dt);
            noalias(rLeftHandSideMatrix) += aux * DampMatrix;
        }
        else
        {
            KRATOS_ERROR << "BaseKinematicLinear::" << __FUNCTION__ << " is not yet implemented for time integration scheme "
                         << rCurrentProcessInfo[TIME_INTEGRATION_SCHEME];
        }
    }

    //************************************************************************************
    //************************************************************************************
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::GetValuesVector( VectorType& values, int Step ) const
    {
        const unsigned int number_of_nodes = this->GetGeometry().size();
        const unsigned int dim = this->WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size )
            values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = i * dim;
            values[index] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, DISPLACEMENT, X ), Step );
            values[index + 1] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, DISPLACEMENT, Y ), Step );
            if ( dim == 3 )
                values[index + 2] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, DISPLACEMENT, Z ), Step );
        }
    }

    //************************************************************************************
    //************************************************************************************
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::GetFirstDerivativesVector( VectorType& values, int Step ) const
    {
        const unsigned int number_of_nodes = this->GetGeometry().size();
        const unsigned int dim = this->WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size )
            values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = i * dim;
            values[index] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, DISPLACEMENT_DT, X ), Step );
            values[index + 1] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, DISPLACEMENT_DT, Y ), Step );
            if ( dim == 3 )
                values[index + 2] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, DISPLACEMENT_DT, Z ), Step );
        }
    }

    //************************************************************************************
    //************************************************************************************
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::GetSecondDerivativesVector( VectorType& values, int Step ) const
    {
        const unsigned int number_of_nodes = this->GetGeometry().size();
        const unsigned int dim = this->WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size )
            values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = i * dim;
            values[index] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, ACCELERATION, X ), Step );
            values[index + 1] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, ACCELERATION, Y ), Step );
            if ( dim == 3 )
                values[index + 2] = this->GetGeometry()[i].GetSolutionStepValue( VARSELC( DataType, ACCELERATION, Z ), Step );
        }
    }

    /**
     * returns the used integration method
     */
    template<typename TNodeType>
    typename BaseKinematicLinear<TNodeType>::IntegrationMethod BaseKinematicLinear<TNodeType>::GetIntegrationMethod() const
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
        else if(this->GetProperties().Has( INTEGRATION_ORDER ))
        {
            if(this->GetProperties()[INTEGRATION_ORDER] == 1)
            {
                return GeometryData::IntegrationMethod::GI_GAUSS_1;
            }
            else if(this->GetProperties()[INTEGRATION_ORDER] == 2)
            {
                return GeometryData::IntegrationMethod::GI_GAUSS_2;
            }
            else if(this->GetProperties()[INTEGRATION_ORDER] == 3)
            {
                return GeometryData::IntegrationMethod::GI_GAUSS_3;
            }
            else if(this->GetProperties()[INTEGRATION_ORDER] == 4)
            {
                return GeometryData::IntegrationMethod::GI_GAUSS_4;
            }
            else if(this->GetProperties()[INTEGRATION_ORDER] == 5)
            {
                return GeometryData::IntegrationMethod::GI_GAUSS_5;
            }
            else
                KRATOS_ERROR << Info() << " does not support for integration order " << this->GetProperties()[INTEGRATION_ORDER];
        }
        else
            return this->GetGeometry().GetDefaultIntegrationMethod(); // default method
    }

    /**
     * Informations for the builder and solver to assemble the global vectors and matrices.
     * Here a VectorType containing the EquationIds of the differnt Dofs is created
     * @param rResult VectorType of the EquationIds
     * @param rCurrentProcessInfo
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::EquationIdVector( EquationIdVectorType& rResult,
            const ProcessInfo& CurrentProcessInfo ) const
    {
        unsigned int dim = ( this->WorkingSpaceDimension() );
        unsigned int mat_size = this->GetGeometry().size() * dim;

        if ( rResult.size() != mat_size )
            rResult.resize( mat_size, false );

        for ( unsigned int i = 0 ; i < this->GetGeometry().size() ; ++i )
        {
            int index = i * dim;
            rResult[index] = this->GetGeometry()[i].GetDof( VARSELC( DataType, DISPLACEMENT, X ) ).EquationId();
            rResult[index+1] = this->GetGeometry()[i].GetDof( VARSELC( DataType, DISPLACEMENT, Y ) ).EquationId();
            if(dim == 3)
                rResult[index+2] = this->GetGeometry()[i].GetDof( VARSELC( DataType, DISPLACEMENT, Z ) ).EquationId();
        }
    }

    /**
     * Informations for the builder and solver to assemble the global vectors and matrices.
     * Here a Container containing the pointers rto the DOFs of this element is created
     * @param ElementalDofList Container with of the DOFs associated with the nodes
     *                           of this element
     * @param rCurrentProcessInfo
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo&
            CurrentProcessInfo ) const
    {
        unsigned int dim = this->WorkingSpaceDimension();

        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0 ; i < this->GetGeometry().size() ; ++i )
        {
            ElementalDofList.push_back( this->GetGeometry()[i].pGetDof( VARSELC( DataType, DISPLACEMENT, X ) ) );
            ElementalDofList.push_back( this->GetGeometry()[i].pGetDof( VARSELC( DataType, DISPLACEMENT, Y ) ) );
            if(dim == 3)
                ElementalDofList.push_back( this->GetGeometry()[i].pGetDof( VARSELC( DataType, DISPLACEMENT, Z ) ) );
        }
    }

    /**
     * Adds the Body Forces to the load vector
     * @param R RHS VectorType
     * @param N_DISP shape function values at the current integration points
     * @param Weight current integration weight
     * @param detJ current Determinant of the Jacobian
     */
    template<typename TNodeType>
    inline void BaseKinematicLinear<TNodeType>::AddBodyForcesToRHS( VectorType& R, const Vector& N_DISP, DataType Weight, DataType detJ ) const
    {
        KRATOS_TRY

        unsigned int dim = this->WorkingSpaceDimension();

        array_1d<ValueType, 3> gravity;
        noalias( gravity ) = this->GetProperties()[GRAVITY];

        ValueType density = 0.0;
        if( this->GetValue( USE_DISTRIBUTED_PROPERTIES ) )
        {
            density = this->GetValue(DENSITY);
        }
        else
        {
            density = this->GetProperties()[DENSITY];
        }

        for ( unsigned int prim = 0; prim < this->GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                R( prim*dim + i ) += N_DISP( prim ) * density * gravity( i ) * detJ * Weight;
            }
        }

        KRATOS_CATCH( "" )
    }

    template<typename TNodeType>
    inline void BaseKinematicLinear<TNodeType>::CalculateAndAdd_ExtForceContribution(
            const Vector& N,
            const ProcessInfo& CurrentProcessInfo,
            const VectorType& BodyForce,
            VectorType& rRightHandSideVector,
            DataType Weight,
            DataType detJ) const
    {
        KRATOS_TRY

        unsigned int number_of_nodes = this->GetGeometry().size();
        unsigned int dimension = this->WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            int index = dimension * i;

            for ( unsigned int j = 0; j < dimension; ++j )
                rRightHandSideVector[index + j] += Weight * detJ * N[i] * BodyForce[j];
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Adds the Internal Forces to the load vector
     * @param R RHS VectorType
     * @param B_Operator B-Operator at the current integration point
     * @param StressVector current stress vector
     * @param Weight current integration weight
     * @param detJ current Determinant of the Jacobian
     */
//    void BaseKinematicLinear<TNodeType>::AddInternalForcesToRHS( VectorType& R, const MatrixType& B_Operator, VectorType& StressVector, ValueType Weight, DataType detJ )
//    {
//        KRATOS_TRY

//        unsigned int dim = this->WorkingSpaceDimension();
//        unsigned int strain_size = dim * (dim + 1) / 2;
//
//        for ( unsigned int prim = 0; prim < this->GetGeometry().size(); ++prim )
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

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::AddInternalForcesToRHS( VectorType& R, const MatrixType& B_Operator, VectorType& StressVector, DataType Weight, DataType detJ ) const
    {
        KRATOS_TRY

        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);
        VectorType InternalForces(3);

        for ( unsigned int prim = 0; prim < this->GetGeometry().size(); ++prim )
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
//            this->GetGeometry()[prim].GetSolutionStepValue( REACTION ) += InternalForces;
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Adds the Contribution of the current quadrature point to the load vector
     * @param K LHS MatrixType
     * @param tan_C 6*6 algorithmic tangent of the materia law (derivation of stresses
     *               regarding strains
     * @param B_Operator B-Operator at the current integration point
     * @param Weight current integration weight
     * @param detJ current Determinant of the Jacobian
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateStiffnesMatrix( MatrixType& K, const MatrixType& tan_C, const MatrixType& B_Operator, DataType Weight, DataType detJ ) const
    {
        KRATOS_TRY

        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);

        for ( unsigned int prim = 0; prim < this->GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                for ( unsigned int sec = 0; sec < this->GetGeometry().size(); ++sec )
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

//        noalias(K) -= prod( trans(B_Operator), (Weight * detJ) * MatrixType(prod(C, B_Operator)) );

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
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateStressAndTangentialStiffness( VectorType& StressVector, MatrixType& tanC_U,
            VectorType& StrainVector, const MatrixType& B_Operator, int PointNumber,
            const ProcessInfo& CurrentProcessInfo ) const
    {
        KRATOS_TRY

        KRATOS_CATCH( "" )
    }

    /**
     * Computes the strain vector
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateStrain( const MatrixType& B, const MatrixType& Displacements, VectorType& StrainVector ) const
    {
        KRATOS_TRY

        unsigned int Dim = this->WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(Dim);
        noalias( StrainVector ) = ZeroVectorType( strain_size );

        for ( unsigned int node = 0; node < this->GetGeometry().size(); ++node )
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
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateBoperator( MatrixType& B_Operator, const Vector& N, const MatrixType& DN_DX ) const
    {
        KRATOS_TRY

        unsigned int dim = this->WorkingSpaceDimension();
        SD_MathUtils<DataType>::CalculateB( dim, B_Operator, DN_DX );

        KRATOS_CATCH( "" )
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateBBaroperator( MatrixType& B_Operator, const MatrixType& DN_DX, const MatrixType& Bdil_bar ) const
    {
        KRATOS_TRY

        unsigned int number_of_nodes = this->GetGeometry().size();
        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);

        noalias( B_Operator ) = ZeroMatrixType( strain_size, number_of_nodes * dim );

        if(dim == 2)
        {
            KRATOS_ERROR << "Bbar formulation for 2D is not supported";
            // const ValueType R1R3 = (1.0 / 3);
            // ValueType tmp1;
            // ValueType tmp2;

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
            constexpr ValueType R1R3 = (1.0 / 3);
            DataType tmp1;
            DataType tmp2;
            DataType tmp3;
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

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateJacobian( typename GeometryType::JacobiansType& J, const IntegrationMethod ThisIntegrationMethod ) const
    {
        typename GeometryType::MatrixType DeltaPosition(this->GetGeometry().size(), 3);

        for ( unsigned int node = 0; node < this->GetGeometry().size(); ++node )
            noalias( row( DeltaPosition, node ) ) = this->GetGeometry()[node].Coordinates()
                                            - this->GetGeometry()[node].GetInitialPosition();

        J = this->GetGeometry().Jacobian( J, ThisIntegrationMethod, DeltaPosition );
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateOnIntegrationPoints( const Variable<MatrixType>& rVariable, std::vector<MatrixType>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        // TODO: This needs to be reviewed (BUI)

        unsigned int number_of_nodes = this->GetGeometry().size();
        unsigned int dim = this->WorkingSpaceDimension();;
        unsigned int strain_size = this->GetStrainSize(dim);

        //Initialize local variables
        MatrixType B( strain_size, number_of_nodes*dim );
        MatrixType TanC( strain_size, strain_size );
        VectorType StrainVector( strain_size );
        VectorType StressVector( strain_size );
        Vector N( number_of_nodes );
        MatrixType DN_DX( number_of_nodes, dim );
        MatrixType CurrentDisp( number_of_nodes, dim );
        MatrixType InvJ0(dim, dim);
        DataType DetJ0;

        //constitutive law
        typename ConstitutiveLawType::Parameters const_params;
        const_params.SetStrainVector(StrainVector);
        const_params.SetStressVector(StressVector);
        const_params.SetConstitutiveMatrix(TanC);
        const_params.SetProcessInfo(rCurrentProcessInfo);
        const_params.SetMaterialProperties(this->GetProperties());
        const_params.SetElementGeometry(this->GetGeometry());
        typename ConstitutiveLawType::StressMeasure stress_measure = ConstitutiveLawType::StressMeasure_Cauchy;
        const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        this->GetGeometry().Initialize(ThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const typename GeometryType::IntegrationPointsArrayType& integration_points =
            this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );
        const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

        if ( rValues.size() != integration_points.size() )
            rValues.resize( integration_points.size() );

        const typename GeometryType::ShapeFunctionsGradientsType& DN_De =
            this->GetGeometry().ShapeFunctionsLocalGradients( ThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        typename GeometryType::JacobiansType J0;
        this->CalculateJacobian( J0, ThisIntegrationMethod );

        //Current displacements
        for ( unsigned int node = 0; node < this->GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = subrange(this->GetGeometry()[node].GetSolutionStepValue( VARSEL( DataType, DISPLACEMENT ) ), 0, dim);

        //Declaration of the integration weight
        //    ValueType Weight;

        //loop over all integration points
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            MathUtils<DataType>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );
            noalias( N ) = row(Ncontainer, PointNumber);

            //Initializing B_Operator at the current integration point
            CalculateBoperator( B, N, DN_DX );

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
        this->GetGeometry().Clean();
        #endif

        KRATOS_CATCH( "" )
    }

    /**
     * Calculate MatrixType Variables at each integration point, used for post-processing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues VectorType to store the values on the quadrature points, output of the method
     * @param rCurrentProcessInfo
     *
     * Calculate VectorType Variables at each integration point, used for post-processing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues VectorType to store the values on the quadrature points, output of the method
     * @param rCurrentProcessInfo
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateOnIntegrationPoints( const Variable<VectorType>& rVariable, std::vector<VectorType>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        unsigned int dim = this->WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);
        unsigned int number_of_nodes = this->GetGeometry().size();
        unsigned int mat_size = dim * number_of_nodes;
        const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

        if ( rVariable == VARSEL(DataType, RECOVERY_STRESSES) )
        {
            /////////////////////////////////////////////////////////////////////////
            //// Calculate recover stresses
            /////////////////////////////////////////////////////////////////////////
            if( this->GetProperties().Has( STRESS_RECOVERY_TYPE ) == true )
            {
                int Type = this->GetProperties()[ STRESS_RECOVERY_TYPE ];

                if(Type == 0)
                {
                    // no recovery
                    this->CalculateOnIntegrationPoints( VARSEL(DataType, STRESSES), rValues, rCurrentProcessInfo );
                }
                else if(Type == 1)
                {
                    // new recovery method from Bathe

                    int ExpansionLevel = this->GetProperties()[ NEIGHBOUR_EXPANSION_LEVEL ];

                    if constexpr (std::is_arithmetic<DataType>::value)
                    {
                        BatheRecoverStressUtility StressUtils(ExpansionLevel);
                        StressUtils.CalculateImprovedStressOnIntegrationPoints( *this, rValues, rCurrentProcessInfo );
                    }
                    else
                        KRATOS_ERROR << "Bathe recovery stress is not implemented for complex number";
                }
                else
                    KRATOS_ERROR << "The stress recovery type " << Type << " is not supported";

            }
            else
                KRATOS_ERROR << "The stress recovery method is not defined for element " << this->Id();
        }
        else if ( rVariable == VARSEL(DataType, STRAIN) || rVariable == VARSEL(DataType, CURRENT_STRAIN_VECTOR) )
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            this->GetGeometry().Initialize(ThisIntegrationMethod);
            #endif

            // calculate shape function values and local gradients
            MatrixType B(strain_size, mat_size);
            VectorType StrainVector(strain_size);
            Vector N(number_of_nodes);
            MatrixType DN_DX(number_of_nodes, dim);
            MatrixType CurrentDisp(number_of_nodes, dim);
            MatrixType InvJ0(dim, dim);
            DataType DetJ0;

            const typename GeometryType::ShapeFunctionsGradientsType& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients( ThisIntegrationMethod );
            const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

            //initializing the Jacobian in the reference configuration
            typename GeometryType::JacobiansType J0;
            this->CalculateJacobian( J0, ThisIntegrationMethod );

            // extract current displacements
            for (unsigned int node = 0; node < this->GetGeometry().size(); ++node)
                noalias(row(CurrentDisp, node)) =
                    subrange(this->GetGeometry()[node].GetSolutionStepValue( VARSEL( DataType, DISPLACEMENT ) ), 0, dim);

            for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            {
                MathUtils<DataType>::InvertMatrix( J0[i], InvJ0, DetJ0 );
                noalias(N) = row(Ncontainer, i);
                noalias(DN_DX) = prod(DN_De[i], InvJ0);

                // compute B_Operator at the current integration point
                CalculateBoperator(B, N, DN_DX);

                // compute the strain at integration point
                CalculateStrain(B, CurrentDisp, StrainVector);

                if (rValues[i].size() != strain_size)
                    rValues[i].resize(strain_size, false);
                noalias(rValues[i]) = StrainVector;
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            this->GetGeometry().Clean();
            #endif
        }
        else if ( rVariable == VARSEL(DataType, PRINCIPAL_STRESS) )
        {
            std::vector<VectorType> stresses;
            this->CalculateOnIntegrationPoints( VARSEL(DataType, STRESSES), stresses, rCurrentProcessInfo );

            MatrixType sigma(3, 3);
            const ValueType conv = 1.e-8;
            const ValueType zero = 1.0e-12;
            for(std::size_t i = 0; i < stresses.size(); ++i)
            {
                if(rValues[i].size() != 3)
                    rValues[i].resize(3, false);

                SD_MathUtils<DataType>::StressVectorToTensor(stresses[i], sigma);

                if constexpr (std::is_arithmetic<DataType>::value)
                {
                    VectorType eigenvalues = SD_MathUtils<DataType>::EigenValuesSym3x3(sigma);
                    noalias(rValues[i]) = SD_MathUtils<DataType>::OrganizeEigenvalues(eigenvalues);
                }
                else
                {
                    KRATOS_ERROR << "Complex eigenvalues are not yet implemented";
                }
            }
        }
        else if ( rVariable == VARSEL(DataType, PRINCIPAL_STRAIN) )
        {
            std::vector<VectorType> strain;
            this->CalculateOnIntegrationPoints( VARSEL(DataType, STRAIN), strain, rCurrentProcessInfo );

            MatrixType epsilon(3, 3);
            const ValueType conv = 1.e-8;
            const ValueType zero = 1.0e-12;
            for(std::size_t i = 0; i < strain.size(); ++i)
            {
                if(rValues[i].size() != 3)
                    rValues[i].resize(3, false);

                SD_MathUtils<DataType>::StrainVectorToTensor(strain[i], epsilon);

                if constexpr (std::is_arithmetic<DataType>::value)
                {
                    VectorType eigenvalues = SD_MathUtils<DataType>::EigenValuesSym3x3(epsilon);
                    noalias(rValues[i]) = SD_MathUtils<DataType>::OrganizeEigenvalues(eigenvalues);
                }
                else
                {
                    KRATOS_ERROR << "Complex eigenvalues are not yet implemented";
                }
            }
        }
        else if ( rVariable == VARSEL(DataType, NODAL_STRESS_VECTOR) )
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            this->GetGeometry().Initialize(ThisIntegrationMethod);
            #endif

            const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

            VectorType StressVector(strain_size);

            for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            {
                noalias(StressVector) = ZeroVectorType(strain_size);
                for (unsigned int j = 0; j < this->GetGeometry().size(); ++j)
                    noalias(StressVector) += Ncontainer(i, j)*this->GetGeometry()[j].GetSolutionStepValue(STRESSES);

                if (rValues[i].size() != strain_size)
                    rValues[i].resize(strain_size);
                noalias(rValues[i]) = StressVector;
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            this->GetGeometry().Clean();
            #endif
        }
        else if ( rVariable == VARSEL(DataType, STRESSES) || rVariable == VARSEL(DataType, PRESTRESS)
               || rVariable == VARSEL(DataType, ELASTIC_STRAIN_VECTOR) || rVariable == VARSEL(DataType, PLASTIC_STRAIN_VECTOR) )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                if (rValues[i].size() != strain_size)
                    rValues[i].resize(strain_size, false);

                rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
            }
        }
        else
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
            }
        }
    }

    /**
     * Calculate double Variables at each integration point, used for post-processing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues VectorType to store the values on the quadrature points, output of the method
     * @param rCurrentProcessInfo
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateOnIntegrationPoints( const Variable<array_1d<DataType, 3> >& rVariable, std::vector<array_1d<DataType, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        const unsigned int dim = this->WorkingSpaceDimension();
        const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        this->GetGeometry().Initialize(ThisIntegrationMethod);
        #endif

        if( rVariable == VARSEL(DataType, DISPLACEMENT) )
        {
            const typename GeometryType::IntegrationPointsArrayType& integration_points =
                    this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

            const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                noalias(rValues[point]) = ZeroVectorType(3);
                for(std::size_t i = 0; i < this->GetGeometry().size(); ++i)
                {
                    const auto& displacement = this->GetGeometry()[i].GetSolutionStepValue(VARSEL( DataType, DISPLACEMENT ));
                    noalias(rValues[point]) += Ncontainer(point, i) * displacement;
                }
            }
        }
        else if( rVariable == VARSEL(DataType, INITIAL_DISPLACEMENT) )
        {
            const typename GeometryType::IntegrationPointsArrayType& integration_points =
                    this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

            const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                noalias(rValues[point]) = ZeroVectorType(3);
                for(std::size_t i = 0; i < this->GetGeometry().size(); ++i)
                {
                    for (unsigned int j = 0; j < dim; ++j)
                        rValues[point][j] += Ncontainer(point, i) * mInitialDisp(i, j);
                }
            }
        }
        else if( rVariable == VARSEL(DataType, INTEGRATION_POINT_GLOBAL) || rVariable == VARSEL(DataType, INTEGRATION_POINT_GLOBAL_IN_CURRENT_CONFIGURATION) )
        {
            const typename GeometryType::IntegrationPointsArrayType& integration_points =
                    this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

            if constexpr (std::is_same<DataType, ValueType>::value)
            {
                for(std::size_t point = 0; point < integration_points.size(); ++point)
                {
                    rValues[point] = this->GetGeometry().GlobalCoordinates(rValues[point], integration_points[point]);
                }
            }
            else
            {
                typename GeometryType::CoordinatesArrayType tmp;
                for(std::size_t point = 0; point < integration_points.size(); ++point)
                {
                    tmp = this->GetGeometry().GlobalCoordinates(tmp, integration_points[point]);
                    rValues[point] = tmp;
                }
            }
        }
        else if( rVariable == VARSEL(DataType, INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION) )
        {
            const typename GeometryType::IntegrationPointsArrayType& integration_points =
                    this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

            Vector N( this->GetGeometry().size() );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                this->GetGeometry().ShapeFunctionsValues( N, integration_points[point] );

                noalias( rValues[point] ) = ZeroVectorType(3);
                for(std::size_t i = 0 ; i < this->GetGeometry().size() ; ++i)
                    noalias( rValues[point] ) += N[i] * this->GetGeometry()[i].GetInitialPosition();
            }
        }
        else if( rVariable == VARSEL(DataType, INTEGRATION_POINT_LOCAL) )
        {
            const typename GeometryType::IntegrationPointsArrayType& integration_points =
                    this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                noalias(rValues[point]) = integration_points[point];
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        this->GetGeometry().Clean();
        #endif
    }

    /**
     * Calculate double Variables at each integration point, used for post-processing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues VectorType to store the values on the quadrature points, output of the method
     * @param rCurrentProcessInfo
     */
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateOnIntegrationPoints( const Variable<DataType>& rVariable, std::vector<DataType>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

        if( rVariable == K0 )
        {
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
            {
                rValues[Point] = this->GetValue( K0 );
            }
            return;
        }
        else if( rVariable == VARSEL(DataType, STRAIN_ENERGY) )
        {
            std::vector<VectorType> StrainList(rValues.size());
//            CalculateOnIntegrationPoints(STRAIN, StrainList, rCurrentProcessInfo);
            CalculateOnIntegrationPoints(VARSEL(DataType, THREED_STRAIN), StrainList, rCurrentProcessInfo);

            std::vector<VectorType> StressList(rValues.size());
//            CalculateOnIntegrationPoints(STRESSES, StressList, rCurrentProcessInfo);
            CalculateOnIntegrationPoints(VARSEL(DataType, THREED_STRESSES), StressList, rCurrentProcessInfo);

            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                // calculate strain energy as C = 0.5 * <epsilon, sigma>
                rValues[i] = 0.5 * inner_prod(StrainList[i], StressList[i]);
            }
        }
        else if( rVariable == VARSEL(DataType, JACOBIAN_0) )
        {
            // initializing the Jacobian in the reference configuration
            typename GeometryType::JacobiansType J0;
            this->CalculateJacobian( J0, ThisIntegrationMethod );

            // compute the Jacobian determinant
            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                rValues[i] = MathUtils<DataType>::Det(J0[i]);
            }
        }
        else if( rVariable == VARSEL(DataType, INTEGRATION_WEIGHT) )
        {
            const typename GeometryType::IntegrationPointsArrayType& integration_points = this->GetGeometry().IntegrationPoints( ThisIntegrationMethod );
            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                rValues[i] = integration_points[i].Weight();
            }
        }
        else if( rVariable == VARSEL(DataType, MATERIAL_DENSITY) )
        {
            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                rValues[i] = this->GetValue(rVariable);
            }
        }
        else
        {
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
            {
                rValues[Point] = mConstitutiveLawVector[Point]->GetValue( rVariable, rValues[Point] );
            }
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if (rValues.size() != mConstitutiveLawVector.size())
            rValues.resize(mConstitutiveLawVector.size());

        if(rVariable == PARENT_ELEMENT_ID)
        {
            std::fill(rValues.begin(), rValues.end(), this->Id());
        }
        else if(rVariable == INTEGRATION_POINT_INDEX)
        {
            for(unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
                rValues[i] = i;
        }
        else
        {
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
                rValues[Point] = mConstitutiveLawVector[Point]->GetValue(rVariable, rValues[Point]);
        }
    }

    #ifndef SD_APP_FORWARD_COMPATIBILITY
    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateOnIntegrationPoints( const Variable<std::string>& rVariable, std::vector<std::string>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
        {
            rValues[Point] = mConstitutiveLawVector[Point]->GetValue(rVariable, rValues[Point]);
        }
    }
    #endif

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::CalculateOnIntegrationPoints( const Variable<typename ConstitutiveLawType::Pointer>& rVariable, std::vector<typename ConstitutiveLawType::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        if( rVariable == VARSEL(ConstitutiveLawType, CONSTITUTIVE_LAW) )
        {
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
            {
                rValues[Point] = mConstitutiveLawVector[Point];
            }
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::SetValuesOnIntegrationPoints( const Variable<MatrixType>& rVariable, const std::vector<MatrixType>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            KRATOS_ERROR << "Error at BaseKinematicLinear element " << this->Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl
                         << "rValues.size(): " << rValues.size() << std::endl
                         << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
        }

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::SetValuesOnIntegrationPoints( const Variable<VectorType>& rVariable, const std::vector<VectorType>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            KRATOS_ERROR << "Error at BaseKinematicLinear element " << this->Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl
                         << "rValues.size(): " << rValues.size() << std::endl
                         << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
        }

        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); ++PointNumber )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::SetValuesOnIntegrationPoints( const Variable<array_1d<DataType, 3> >& rVariable, const std::vector<array_1d<DataType, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if (rVariable == INTEGRATION_POINT_LOCAL)
            KRATOS_ERROR << "BaseKinematicLinear element does not support for given natural coordinates of the integration point" << std::endl;

        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            KRATOS_ERROR << "Error at BaseKinematicLinear element " << this->Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl
                         << "rValues.size(): " << rValues.size() << std::endl
                         << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
        }

        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); ++PointNumber )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::SetValuesOnIntegrationPoints( const Variable<DataType>& rVariable, const std::vector<DataType>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rVariable == VARSEL(DataType, K0) )
        {
            this->SetValue( VARSEL(DataType, K0), rValues[0] );
        }
        else
        {
            if ( rValues.size() != mConstitutiveLawVector.size() )
            {
                KRATOS_ERROR << "Error at BaseKinematicLinear element " << this->Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl
                             << "rValues.size(): " << rValues.size() << std::endl
                             << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
            }

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
            }
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::SetValuesOnIntegrationPoints( const Variable<int>& rVariable, const std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            KRATOS_ERROR << "Error at BaseKinematicLinear element " << this->Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl
                         << "rValues.size(): " << rValues.size() << std::endl
                         << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
        }

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::SetValuesOnIntegrationPoints( const Variable<bool>& rVariable, const std::vector<bool>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            KRATOS_ERROR << "Error at BaseKinematicLinear element " << this->Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl
                         << "rValues.size(): " << rValues.size() << std::endl
                         << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
        }

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
        }
    }

    template<typename TNodeType>
    void BaseKinematicLinear<TNodeType>::SetValuesOnIntegrationPoints( const Variable<typename ConstitutiveLawType::Pointer>& rVariable, const std::vector< typename ConstitutiveLawType::Pointer >& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

        if ( rVariable == VARSEL(ConstitutiveLawType, CONSTITUTIVE_LAW) )
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            this->GetGeometry().Initialize(ThisIntegrationMethod);
            #endif

            if( mConstitutiveLawVector.size() != rValues.size() )
            {
                mConstitutiveLawVector.resize( rValues.size() );
            }

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i] = rValues[i];
                mConstitutiveLawVector[i]->InitializeMaterial( this->GetProperties(), this->GetGeometry(), row( this->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod ), i ) );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            this->GetGeometry().Clean();
            #endif
        }
        else if ( rVariable == VARSEL(ConstitutiveLawType, CONSTITUTIVE_LAW_NO_INITIALIZE) )
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            this->GetGeometry().Initialize(ThisIntegrationMethod);
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
            this->GetGeometry().Clean();
            #endif
        }
    }

    template<typename TNodeType>
    int BaseKinematicLinear<TNodeType>::Check( const Kratos::ProcessInfo& rCurrentProcessInfo ) const
    {
        KRATOS_TRY

        if ( this->Id() < 1 )
        {
            KRATOS_ERROR << "BaseType found with Id 0 or negative";
        }

        //verify that the constitutive law exists
        if ( this->GetProperties().Has( VARSEL(ConstitutiveLawType, CONSTITUTIVE_LAW) ) == false )
        {
            KRATOS_ERROR << "constitutive law not provided for property " << this->GetProperties().Id();
        }

        // verify the strain measure
        // auto strain_mearure = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainMeasure();
        // if ( strain_mearure != typename ConstitutiveLawType::StrainMeasure_Infinitesimal
        //   && strain_mearure != typename ConstitutiveLawType::StrainMeasure_GreenLagrange )
        // {
        //     std::stringstream ss;
        //     ss << "The strain measure " << strain_mearure << " is not supported by this element";
        //     KRATOS_THROW_ERROR( std::logic_error, ss.str(), "" )
        // }
        typename ConstitutiveLawType::Features features;
        this->GetProperties().GetValue( VARSEL(ConstitutiveLawType, CONSTITUTIVE_LAW) )->GetLawFeatures(features);
        if ( std::find(features.GetStrainMeasures().begin(), features.GetStrainMeasures().end(), ConstitutiveLawType::StrainMeasure_Infinitesimal) == features.GetStrainMeasures().end()
          // && std::find(features.GetStrainMeasures().begin(), features.GetStrainMeasures().end(), typename ConstitutiveLawType::StrainMeasure_GreenLagrange) == features.GetStrainMeasures().end()
           )
        {
            KRATOS_ERROR << "The constitutive law strain measures are not supported by this element";
        }

        //Verify that the body force is defined
        if ( this->GetProperties().Has( BODY_FORCE ) == false )
        {
            KRATOS_ERROR << "BODY_FORCE not provided for property " << this->GetProperties().Id();
        }

        return 0;

        KRATOS_CATCH( "" );
    }

    /// template class instantiation
    template class BaseKinematicLinear<RealNode>;
    template class BaseKinematicLinear<ComplexNode>;
    template class BaseKinematicLinear<GComplexNode>;

} // Namespace Kratos

#undef ENABLE_DEBUG_CONSTITUTIVE_LAW
#undef COMPUTE_NUMERICAL_DERIVATIVE
#undef DEBUG_ELEMENT_ID
