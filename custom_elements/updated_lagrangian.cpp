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
 *   Date:                $Date: 13 Sep 2023 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

// System includes


// External includes

// Project includes
#include "custom_elements/updated_lagrangian.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    UpdatedLagrangian::UpdatedLagrangian( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseType( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
        //THIS IS THE DEFAULT CONSTRUCTOR
    }

    UpdatedLagrangian::UpdatedLagrangian( IndexType NewId,
            GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseType( NewId, pGeometry, pProperties )
    {}

    Element::Pointer UpdatedLagrangian::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new UpdatedLagrangian( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    Element::Pointer UpdatedLagrangian::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new UpdatedLagrangian( NewId, pGeom, pProperties ) );
    }

    UpdatedLagrangian::~UpdatedLagrangian()
    {
    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangian::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();

        BaseType::Initialize(rCurrentProcessInfo);

        mLastF.resize(mConstitutiveLawVector.size());
        mCurrentF.resize(mConstitutiveLawVector.size());

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        {
            mLastF[i].resize(dim, dim, false);
            noalias(mLastF[i]) = IdentityMatrix(dim);
            mCurrentF[i].resize(dim, dim, false);
            noalias(mCurrentF[i]) = IdentityMatrix(dim);
        }
    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangian::CalculateAll( MatrixType& rLeftHandSideMatrix,
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

        Matrix F( f_size, f_size ), Fd( f_size, f_size );

        Matrix InvF( f_size, f_size );

        Matrix A( g_size, g_size );

        Matrix C( dim, dim );

        Vector StrainVector( strain_size );

        Vector StressVector( strain_size );

        Matrix DN_DX( number_of_nodes, dim );
        Vector N( number_of_nodes );
        Matrix DN_Dx( number_of_nodes, dim );

        Matrix CurrentDisp( number_of_nodes, 3 );

        Matrix InvJn(dim, dim), InvJ(dim, dim);

        double DetJn, DetJ;

        double DetF;

        Matrix InvFd(f_size, f_size);

        double DetFd;

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

        // current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        // calculating Jacobian in previous state
        GeometryType::JacobiansType Jn;
        Jn = GetGeometry().Jacobian( Jn, mThisIntegrationMethod );

        // calculating Jacobian in current configuration
        GeometryType::JacobiansType J;
        Matrix ShiftedDisp = -CurrentDisp;
        J = GetGeometry().Jacobian( J, mThisIntegrationMethod, ShiftedDisp );

        // /////////////////////////////////////////////////////////////////////////
        // //// Compute FO, GO for Fbar formulation
        // //// Reference: Souze de Neto, Computational Plasticity, Box 15.1
        // /////////////////////////////////////////////////////////////////////////
        // const bool is_fbar = GetProperties().Has(IS_FBAR) ? GetProperties()[IS_FBAR] : false;
        // Matrix FO, GOx, Q, A3d, eye, stress_tensor;
        // SD_MathUtils<double>::Fourth_Order_Tensor Atensor, Qtensor, ItimesI;
        // double DetFO;
        // if (is_fbar)
        // {
        //     CoordinatesArrayType origin;
        //     GeometryUtility::ComputeOrigin(GetGeometry().GetGeometryFamily(), origin);

        //     Matrix JnO, JO;
        //     JnO = GetGeometry().Jacobian( JnO, origin, DeltaPosition );
        //     JO = GetGeometry().Jacobian( JO, origin, ShiftedDisp );

        //     Vector NO( number_of_nodes );
        //     NO = GetGeometry().ShapeFunctionsValues(NO, origin);

        //     Matrix DN_De_O( number_of_nodes, dim );
        //     DN_De_O = GetGeometry().ShapeFunctionsLocalGradients(DN_De_O, origin);

        //     Matrix InvJnO(dim, dim), InvJO(dim, dim);
        //     double DetJnO, DetJO;
        //     Matrix DN_Dx_O( number_of_nodes, dim );

        //     MathUtils<double>::InvertMatrix( JnO, InvJnO, DetJnO );
        //     MathUtils<double>::InvertMatrix( JO, InvJO, DetJO );
        //     noalias( DN_Dx_O ) = prod( DN_De_O, InvJO );

        //     // FO
        //     Matrix DN_DX_O( number_of_nodes, dim );
        //     noalias( DN_DX_O ) = prod( DN_De_O, InvJnO );

        //     Matrix GOX(g_size, number_of_nodes * dim);
        //     this->CalculateG( GOX, NO, DN_DX_O );

        //     FO.resize(f_size, f_size, false);
        //     this->CalculateF( FO, GOX, CurrentDisp );
        //     DetFO = MathUtils<double>::Det(FO);

        //     #ifdef CHECK_DEFORMATION_GRADIENT
        //     if (DetFO < 0.0)
        //     {
        //         KRATOS_WATCH(origin)
        //         KRATOS_WATCH(JnO)
        //         KRATOS_WATCH(InvJnO)
        //         KRATOS_WATCH(JO)
        //         KRATOS_WATCH(FO)
        //         KRATOS_WATCH(DetFO)
        //         KRATOS_ERROR << "Deformation gradient is negative at origin of element " << Id();
        //     }
        //     #endif

        //     // GOx
        //     GOx.resize( g_size, number_of_nodes * dim, false );
        //     this->CalculateG( GOx, NO, DN_Dx_O, CurrentDisp );

        //     #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
        //     if (Id() == DEBUG_ELEMENT_ID)
        //     {
        //         KRATOS_WATCH(GOx)
        //     }
        //     #endif

        //     // resize Q as needed
        //     Q.resize(g_size, g_size, false);

        //     // resize A3d as needed
        //     A3d.resize(9, 9, false);

        //     // initialize the tensors
        //     SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Atensor);
        //     SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Qtensor);

        //     eye = IdentityMatrix(3);
        //     if (f_size == 2) eye(2, 2) = 0.0; // this modification is required because we don't want to involve the zz component
        //     SD_MathUtils<double>::CalculateFourthOrderZeroTensor(ItimesI);
        //     SD_MathUtils<double>::OuterProductFourthOrderTensor(1.0, eye, eye, ItimesI);

        //     stress_tensor.resize(3, 3);
        // }

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
            MathUtils<double>::InvertMatrix( Jn[PointNumber], InvJn, DetJn );
            MathUtils<double>::InvertMatrix( J[PointNumber], InvJ, DetJ );
            noalias( N ) = row( Ncontainer, PointNumber );
            noalias( DN_Dx ) = prod( DN_De[PointNumber], InvJ );
            noalias( DN_DX ) = prod( DN_De[PointNumber], InvJn );

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
            if (Id() == DEBUG_ELEMENT_ID)
            {
                KRATOS_WATCH(integration_points[PointNumber])
                KRATOS_WATCH(DN_Dx)
                // KRATOS_WATCH(DN_DX)
            }
            #endif

            //deformation gradient
            this->CalculateG( GX, N, DN_DX );
            this->CalculateF( Fd, GX, CurrentDisp );
            MathUtils<double>::InvertMatrix( Fd, InvFd, DetFd );

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

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
            // Matrix Finv(dim, dim);
            // for (unsigned int i = 0; i < dim; ++i)
            // {
            //     for (unsigned int j = 0; j < dim; ++j)
            //     {
            //         Finv(i, j) = 0.0;
            //         for (unsigned int n = 0; n < GetGeometry().size(); ++n)
            //             Finv(i, j) -= DN_Dx(n, j) * GetGeometry()[n].GetSolutionStepValue(DISPLACEMENT)[i];
            //     }
            //     Finv(i, i) += 1.0;
            // }
            // KRATOS_WATCH(prod(F, Finv)) // should produce identity matrix
            #endif

            noalias(F) = prod(Fd, mLastF[PointNumber]);
            DetF = MathUtils<double>::Det(F);
            noalias(mCurrentF[PointNumber]) = F;

            #ifdef CHECK_DEFORMATION_GRADIENT
            if (DetF < 0.0)
            {
                KRATOS_WATCH(F)
                KRATOS_ERROR << "Deformation gradient is negative at integration point " << PointNumber
                             << " of element " << Id();
            }
            #endif

            // if (Id() == 20)
            // {
            //     KRATOS_WATCH(F)
            //     KRATOS_WATCH(DetFO)
            //     KRATOS_WATCH(DetF)
            // }

            // if (is_fbar)
            // {
            //     if (f_size == 2) // plane strain
            //     {
            //         F *= std::sqrt(DetFO/DetF);
            //         #ifdef USE_DETERMINANT_FO_FOR_FBAR
            //         const_params.SetDeterminantF(DetFO); // for the reason why to do this, see iffba2.f
            //         #else
            //         const_params.SetDeterminantF(std::sqrt(DetFO/DetF)*DetF);
            //         #endif
            //     }
            //     else if (f_size == 3) // plane stress, axisymmetric, 3D
            //     {
            //         F *= std::cbrt(DetFO/DetF);
            //         #ifdef USE_DETERMINANT_FO_FOR_FBAR
            //         const_params.SetDeterminantF(DetFO); // for the reason why to do this, see iffba2.f
            //         #else
            //         const_params.SetDeterminantF(std::cbrt(DetFO/DetF)*DetF);
            //         #endif
            //     }
            //     else
            //     {
            //         KRATOS_ERROR << "Invalid size " << f_size << " of deformation gradient";
            //     }
            // }
            // else
                const_params.SetDeterminantF(DetF);

            //strain calculation
            noalias( C ) = prod( trans( F ), F );

            BaseType::CalculateStrain( C, StrainVector );

            //integrate the material
            BaseType::mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse( const_params, stress_measure );

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight = this->GetIntegrationWeight(integration_points[PointNumber].Weight(), N, CurrentDisp) * DetJn;

            if ( dim == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

            if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
            {
                //calculating operator B
                BaseType::CalculateB( B, N, prod(DN_DX, InvFd) * DetFd, CurrentDisp );

                //contribution to external forces
                const Vector& BodyForce = GetProperties()[BODY_FORCE];
                BaseType::CalculateAndAdd_ExtForceContribution( N, rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight );

                //contribution of gravity (if there is)
                BaseType::AddBodyForcesToRHS( rRightHandSideVector, N, IntToReferenceWeight );

                //contribution of internal forces
                noalias( rRightHandSideVector ) -= IntToReferenceWeight * prod( trans( B ), StressVector );

                #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
                if (Id() == DEBUG_ELEMENT_ID)
                {
                    KRATOS_WATCH(B)
                    KRATOS_WATCH(StressVector)
                    KRATOS_WATCH(rRightHandSideVector)
                }
                #endif
            }

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                //calculating operator G
                this->CalculateG( Gx, N, DN_Dx, CurrentDisp );

                //contributions to stiffness
                noalias( rLeftHandSideMatrix ) += prod( trans( Gx ), ( IntToReferenceWeight ) * Matrix( prod( A, Gx ) ) ); //to be optimized to remove the temporary

                // if (is_fbar)
                // {
                //     //compute Q matrix
                //     mConstitutiveLawVector[PointNumber]->GetValue(CAUCHY_STRESS_TENSOR, stress_tensor);
                //     if ( dim == 2)
                //     {
                //         mConstitutiveLawVector[PointNumber]->GetValue(THREED_ALGORITHMIC_TANGENT, A3d);
                //         SD_MathUtils<double>::UnsymmetricMatrixToTensor(A3d, Atensor);
                //     }
                //     else if (dim == 3)
                //     {
                //         SD_MathUtils<double>::UnsymmetricMatrixToTensor(A, Atensor);
                //     }
                //     #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
                //     if (Id() == DEBUG_ELEMENT_ID)
                //     {
                //         KRATOS_WATCH(stress_tensor)
                //     //     KRATOS_WATCH(A3d)
                //     //     KRATOS_WATCH(Atensor)

                //     //     Matrix ItimesImat(g_size, g_size);
                //     //     SD_MathUtils<double>::TensorToUnsymmetricMatrix(ItimesI, ItimesImat);
                //     //     KRATOS_WATCH(ItimesImat)
                //     //     Matrix AxItimesI = prod(A, ItimesImat);
                //     //     KRATOS_WATCH(0.5*AxItimesI)
                //     }
                //     #endif

                //     SD_MathUtils<double>::ZeroFourthOrderTensor(Qtensor);
                //     if (f_size == 2) // plane strain
                //     {
                //         SD_MathUtils<double>::ProductFourthOrderTensor(0.5, Atensor, ItimesI, Qtensor);
                //         SD_MathUtils<double>::OuterProductFourthOrderTensor(-0.5, stress_tensor, eye, Qtensor);
                //     }
                //     else if (f_size == 3) // plane stress, axisymmetric, 3D
                //     {
                //         SD_MathUtils<double>::ProductFourthOrderTensor(1.0/3, Atensor, ItimesI, Qtensor);
                //         SD_MathUtils<double>::OuterProductFourthOrderTensor(-2.0/3, stress_tensor, eye, Qtensor);
                //     }
                //     else
                //         KRATOS_ERROR << "Invalid size " << f_size << " of deformation gradient";

                //     SD_MathUtils<double>::TensorToUnsymmetricMatrix(Qtensor, Q);

                //     //contributions to stiffness
                //     noalias( rLeftHandSideMatrix ) += prod( trans( Gx ), ( IntToReferenceWeight ) * Matrix( prod( Q, GOx - Gx ) ) );
                // }

                #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
                if (Id() == DEBUG_ELEMENT_ID)
                {
                    KRATOS_WATCH(A)
                    if (is_fbar)
                    {
                        // KRATOS_WATCH(Qtensor)
                        KRATOS_WATCH(Q)
                        KRATOS_WATCH(GOx)
                    }
                    KRATOS_WATCH(Gx)
                    KRATOS_WATCH(rLeftHandSideMatrix)
                }
                #endif

                #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
                if (Id() == DEBUG_ELEMENT_ID)
                {
                    KRATOS_WATCH("-----------")
                }
                #endif
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //clean the internal data of the geometry
        GetGeometry().Clean();
        #endif

        #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
        if (Id() == DEBUG_ELEMENT_ID)
        {
            // KRATOS_WATCH(rRightHandSideVector)
            // KRATOS_WATCH(rLeftHandSideMatrix)
            KRATOS_WATCH("----------------------")
        }
        #endif

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangian::FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo )
    {
        BaseType::FinalizeSolutionStep( CurrentProcessInfo );

        for (unsigned int i = 0; i < mLastF.size(); ++i)
            noalias(mLastF[i]) = mCurrentF[i];
    }

} // Namespace Kratos
