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
 *   Last Modified by:    $Author: vladislav $
 *   Date:                $Date: 2015-07-28 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

// System includes

// External includes
#include "boost/math/constants/constants.hpp"

// Project includes
#include "timoshenko_beam_element.h"
#include "structural_application_variables.h"

#define PI boost::math::constants::pi<double>()

namespace Kratos
{

    //*****************************************************************************
    //*****************************************************************************

    // Constructors
    TimoshenkoBeamElement::TimoshenkoBeamElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
    {
    }

    TimoshenkoBeamElement::TimoshenkoBeamElement( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
    {
        KRATOS_TRY

        unsigned int dimension = GetGeometry().WorkingSpaceDimension();  // checks dimension of problem

        if (dimension != 3)
            KRATOS_ERROR << "This element only works in 3D space";

        KRATOS_CATCH("")
    }

    Element::Pointer TimoshenkoBeamElement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new TimoshenkoBeamElement( NewId, pGeom, pProperties ) );
    }

    Element::Pointer TimoshenkoBeamElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer( new TimoshenkoBeamElement( NewId, GetGeometry().Create(ThisNodes), pProperties ) );
    }

    // Destructor
    TimoshenkoBeamElement::~TimoshenkoBeamElement()
    {
    }

    //************************************************************************************
    //THIS IS THE INITIALIZATION OF THE ELEMENT (CALLED AT THE BEGIN OF EACH CALCULATION)
    //************************************************************************************

    void TimoshenkoBeamElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        CalculateSectionProperties();

        unsigned int number_of_nodes = GetGeometry().size();

        if (number_of_nodes == 2)
        {
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
        }
        else if (number_of_nodes == 3)
        {
            mThisIntegrationMethod = IntegrationMethod::GI_GAUSS_3;
        }

        if (rCurrentProcessInfo[RESET_CONFIGURATION] == 0)
        {
            mInitialDisp.resize(number_of_nodes, 3, false);
            noalias(mInitialDisp) = ZeroMatrix(number_of_nodes, 3);
            mInitialRot.resize(number_of_nodes, 3, false);
            noalias(mInitialRot) = ZeroMatrix(number_of_nodes, 3);
        }
        else if (rCurrentProcessInfo[RESET_CONFIGURATION] == 1)
        {
            if (mInitialDisp.size1() != number_of_nodes || mInitialDisp.size2() != 3)
                mInitialDisp.resize( number_of_nodes, 3, false );

            if (mInitialRot.size1() != number_of_nodes || mInitialRot.size2() != 3)
                mInitialRot.resize( number_of_nodes, 3, false );

            for ( unsigned int node = 0; node < number_of_nodes; ++node )
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

    void TimoshenkoBeamElement::InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo)
    {

    }

    //************************************************************************************
    //************************************************************************************

    void TimoshenkoBeamElement::FinalizeSolutionStep(const ProcessInfo& CurrentProcessInfo)
    {

    }

    //*************************************************************************************
    //*************************************************************************************

    void TimoshenkoBeamElement::CalculateSectionProperties()
    {
        KRATOS_TRY

        //**********************************************************************************
        // Initializing area
        //**********************************************************************************
        mArea = GetProperties()[AREA];
        //**********************************************************************************
        // Initializing effective shear area for rectangular cross section
        //**********************************************************************************
        mArea_y = GetProperties()[AREA_Y];
        mArea_z = GetProperties()[AREA_Z];
        //**********************************************************************************
        // Calculating moments of inertia
        //**********************************************************************************
        mInertia_x = GetProperties()[INERTIA_X];
        mInertia_y = GetProperties()[INERTIA_Y];
        mInertia_z = GetProperties()[INERTIA_Z];

        KRATOS_CATCH( "" )
    }

    //*************************************************************************************
    //*************************************************************************************

    void TimoshenkoBeamElement::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo&
            CurrentProcessInfo) const
    {
        unsigned int number_of_nodes = GetGeometry().size();
        ElementalDofList.resize(0);

        if (number_of_nodes == 2)
        {
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
        else if (number_of_nodes == 3)
        {
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
            ElementalDofList.push_back(GetGeometry()[2].pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(GetGeometry()[2].pGetDof(DISPLACEMENT_Y));
            ElementalDofList.push_back(GetGeometry()[2].pGetDof(DISPLACEMENT_Z));
            ElementalDofList.push_back(GetGeometry()[2].pGetDof(ROTATION_X));
            ElementalDofList.push_back(GetGeometry()[2].pGetDof(ROTATION_Y));
            ElementalDofList.push_back(GetGeometry()[2].pGetDof(ROTATION_Z));
        }
    }

    //************************************************************************************
    //************************************************************************************

    void TimoshenkoBeamElement::EquationIdVector(EquationIdVectorType& rResult,
            const ProcessInfo& CurrentProcessInfo) const
    {
        unsigned int number_of_nodes = GetGeometry().size();
        if (number_of_nodes == 2)
        {
            if(rResult.size() != 12)
                rResult.resize(12, false);

            rResult[0]    = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
            rResult[1]    = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[2]    = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[3]    = GetGeometry()[0].GetDof(ROTATION_X).EquationId();
            rResult[4]    = GetGeometry()[0].GetDof(ROTATION_Y).EquationId();
            rResult[5]    = GetGeometry()[0].GetDof(ROTATION_Z).EquationId();
            rResult[6]    = GetGeometry()[1].GetDof(DISPLACEMENT_X).EquationId();
            rResult[7]    = GetGeometry()[1].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[8]    = GetGeometry()[1].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[9]    = GetGeometry()[1].GetDof(ROTATION_X).EquationId();
            rResult[10]   = GetGeometry()[1].GetDof(ROTATION_Y).EquationId();
            rResult[11]   = GetGeometry()[1].GetDof(ROTATION_Z).EquationId();
        }
        else if (number_of_nodes == 3)
        {
            if(rResult.size() != 18)
                rResult.resize(18, false);

            rResult[0]    = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
            rResult[1]    = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[2]    = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[3]    = GetGeometry()[0].GetDof(ROTATION_X).EquationId();
            rResult[4]    = GetGeometry()[0].GetDof(ROTATION_Y).EquationId();
            rResult[5]    = GetGeometry()[0].GetDof(ROTATION_Z).EquationId();
            rResult[6]    = GetGeometry()[1].GetDof(DISPLACEMENT_X).EquationId();
            rResult[7]    = GetGeometry()[1].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[8]    = GetGeometry()[1].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[9]    = GetGeometry()[1].GetDof(ROTATION_X).EquationId();
            rResult[10]   = GetGeometry()[1].GetDof(ROTATION_Y).EquationId();
            rResult[11]   = GetGeometry()[1].GetDof(ROTATION_Z).EquationId();
            rResult[12]   = GetGeometry()[2].GetDof(DISPLACEMENT_X).EquationId();
            rResult[13]   = GetGeometry()[2].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[14]   = GetGeometry()[2].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[15]   = GetGeometry()[2].GetDof(ROTATION_X).EquationId();
            rResult[16]   = GetGeometry()[2].GetDof(ROTATION_Y).EquationId();
            rResult[17]   = GetGeometry()[2].GetDof(ROTATION_Z).EquationId();
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
    int  TimoshenkoBeamElement::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        // verify that the area is given by properties
        if (this->GetProperties().Has(AREA) == false)
        {
            // if( this->GetProperties()[AREA] == 0.0 )
                KRATOS_ERROR << "AREA not provided for this element " << this->Id();
        }

        // verify that the inertia is given by properties
        if (this->GetProperties().Has(INERTIA_X) == false)
        {
            // if( GetProperties()[INERTIA_X] == 0.0 )
                KRATOS_ERROR << "INERTIA_X not provided for this element " << this->Id();
        }

        if (this->GetProperties().Has(INERTIA_Y) == false)
        {
            // if( GetProperties()[INERTIA_Y] == 0.0 )
                KRATOS_ERROR << "INERTIA_Y not provided for this element " << this->Id();
        }
        if (this->GetProperties().Has(INERTIA_Z) == false)
        {
            // if( GetProperties()[INERTIA_Z] == 0.0 )
                KRATOS_ERROR << "INERTIA_Z not provided for this element " << this->Id();
        }

        return 0;

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void TimoshenkoBeamElement::CalculateExternalLoadVector(Matrix& Rotation, Vector& LocalBody, Vector& GlobalBody) const
    {
        KRATOS_TRY

        //************************************************************************************
        //Calculates external load vector
        //
        //Works only for constant distributed load - should be extended for arbitrary distibution
        // Also extend for integration over Gauss points
        //************************************************************************************

        double alpha =  0.00;
        double sign  =  1.00;
        const double Length = GetGeometry().Length();
        double  sine;
        double  cosine;

        array_1d<double, 3> Weight;
        const Vector& Gravity = GetProperties()[GRAVITY];
        const double density = GetProperties()[DENSITY];
        // Vector Gravity(3);
        // double density = 0.0;
        noalias(Weight) = density * Gravity;

        array_1d<double, 3> Distr_load;
        noalias(Distr_load) = mArea * Weight;

        //KRATOS_WATCH(Distr_load)

        array_1d<double, 12> Load_X = ZeroVector(12);
        array_1d<double, 12> Load_Y = ZeroVector(12);
        array_1d<double, 12> Load_Z = ZeroVector(12);

        array_1d<double, 2> Load;
        array_1d<double, 6> x_zero;

        Vector Normal_Loads;
        Vector Vector_zero;

        Normal_Loads.resize(3, false);
        Vector_zero.resize(3, false);

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

        //Load in X direction
        //***********************************
        if(Distr_load[0]!=0.00)
        {
            Normal_Loads[0]   = 0.00;
            Normal_Loads[1]   = Vector_zero[1] ;
            Normal_Loads[2]   = Vector_zero[2] ;

            if (Vector_zero[0] < 0)
            {
                sign = -1.00;
            }
            if( norm_2(Normal_Loads)==0 || norm_2(Vector_zero)==0  )
            {
                alpha = sign*PI/2;
            }
            else
            {
                alpha = inner_prod(Normal_Loads, Vector_zero)/(norm_2(Vector_zero)*norm_2(Normal_Loads));
                alpha = sign*acos(alpha);
            }

            sine = sin(alpha);
            cosine = cos(alpha);
            if(fabs(sine) < 1E-7) sine = 0.00;
            if(fabs(cosine) < 1E-7) cosine = 0.00;

            Load[0]= Distr_load[0]*sine;
            Load[1]= Distr_load[0]*cosine;

            Load_X[0]= Load[0]*Length/2.00;
            Load_X[1]= -(Load[1]*Length)/2.00;
            Load_X[2]= 0.00;

            Load_X[3]= 0.00;
            Load_X[4]= 0.00;
            Load_X[5]= -(Load[1])*Length*Length/12.00;;
            Load_X[6]= Load[0]*Length/2.00;
            Load_X[7]= -(Load[1])*Length/2.00;
            Load_X[8]= 0.00;
            Load_X[9]= 0.00;
            Load_X[10]= 0.00;
            Load_X[11]= (Load[1])*Length*Length/12.00;

            noalias(GlobalBody) -= prod(Rotation,Load_X);
            noalias(LocalBody)  -= Load_X;
        }

        //Load in Y direction
        //***********************************
        if(Distr_load[1]!=0.00)
        {
            Normal_Loads[0]    = Vector_zero[0] ;
            Normal_Loads[1]    = 0.00 ;
            Normal_Loads[2]    = Vector_zero[2];

            if (Vector_zero[1]<0)
            {
                sign =-1.00;
            }
            if( norm_2(Normal_Loads)==0 || norm_2( Vector_zero)==0  )
            {
                alpha = sign*PI/2;
            }
            else
            {
                alpha = inner_prod(Normal_Loads,Vector_zero)/(norm_2(Vector_zero)*norm_2( Normal_Loads));
                alpha = sign*acos(alpha);
            }

            sine = sin(alpha);
            cosine = cos(alpha);

            if(fabs(sine) < 1E-7) sine = 0.00;
            if(fabs(cosine) < 1E-7) cosine = 0.00;

            Load[0] = Distr_load[1]*sine;
            Load[1] = Distr_load[1]*cosine;

            Load_Y[0] = -Load[0]*Length/2.00;
            Load_Y[1] = -(Load[1]*Length)/2.00;
            Load_Y[2] = 0.00;

            Load_Y[3] = 0.00;

            Load_Y[4] = 0.00;

            Load_Y[5] = -(Load[1])*Length*Length/12.00;
            Load_Y[6] = -(Load[0])*Length/2.00;
            Load_Y[7] = -(Load[1])*Length/2.00;
            Load_Y[8] = 0.00;
            Load_Y[9] = 0.00;
            Load_Y[10] = 0.00;
            Load_Y[11] = (Load[1])*Length*Length/12.00;

            noalias(GlobalBody) -= prod(Rotation, Load_Y);
            noalias(LocalBody)  -= Load_Y;

            //std::cout << "load[0] " << sine << " load[1] " << cosine << std::endl;
        }

        //Load in Z direction
        //***********************************
        if(Distr_load[2]!=0.00)
        {

            Normal_Loads[0]    = Vector_zero[0] ;
            Normal_Loads[1]    = Vector_zero[1] ;
            Normal_Loads[2]    = 0.00;

            if (Vector_zero[2]<0)
            {
                sign =-1.00;
            }
            if( norm_2(Normal_Loads)==0 || norm_2( Vector_zero)==0  )
            {
                alpha = sign*PI/2;
            }
            else
            {
                alpha = inner_prod(Normal_Loads,Vector_zero)/(norm_2(Vector_zero)*norm_2( Normal_Loads));
                alpha    = sign*acos(alpha);
            }

            sine = sin(alpha);
            cosine = cos(alpha);

            if(fabs(sine) < 1E-7) sine = 0.00;
            if(fabs(cosine) < 1E-7) cosine = 0.00;


            Load[0]= Distr_load[2]*sine;  // load in axial direction
            Load[1]= Distr_load[2]*cosine; //

            Load_Z[0]= -Load[0]*Length/2.00;
            Load_Z[1]= 0.00;
            Load_Z[2]= -(Load[1]*Length)/2.00;
            Load_Z[3]= 0.00;
            Load_Z[4]= -Load[1]*Length*Length/12.00;

            Load_Z[5]= 0.00;
            Load_Z[6]= -Load[0]*Length/2.00;
            Load_Z[7]= 0.00;
            Load_Z[8]= -(Load[1])*Length/2.00;
            Load_Z[9]=  0.00;
            Load_Z[10]= (Load[1])*Length*Length/12.00;
            Load_Z[11]=  0.00;

            noalias(GlobalBody) -= prod(Rotation, Load_Z); //External load vector in global coordinates
            noalias(LocalBody)  -= Load_Z;
        }

        //KRATOS_WATCH(GlobalBody)
        KRATOS_CATCH("")
    }

    //*****************************************************************************
    //*****************************************************************************

    void TimoshenkoBeamElement::CalculateTransformationMatrix(Matrix& Rotation, const Vector& TangentialVector) const
    {
        KRATOS_TRY

        unsigned int number_of_nodes = GetGeometry().size();
        Vector Normal_zero(9); // vector containing directional cosines

        noalias(Normal_zero) =    zero_vector<double>(9);
        noalias(Rotation)    =    zero_matrix<double>(6*number_of_nodes, 6*number_of_nodes);

        double nx, ny, nz, theta;

        for ( unsigned int i = 0; i < 3; ++i )
        {
            Normal_zero[i] = TangentialVector[i];
        }

        nx = Normal_zero[0];
        ny = Normal_zero[1];
        nz = Normal_zero[2];

        if (nx ==0.0)
        {
            theta = PI/2;
            if (ny == 0.0)
            {
                theta = 0.0;
            }
        }
        else
        {
            theta = atan(ny/nx);
        }

        if(nx < 0.0)
            theta = theta + PI;

        Normal_zero[3] = -sin(theta);
        Normal_zero[4] =  cos(theta);
        Normal_zero[5] =  0.0;
        Normal_zero[6] = -nz*cos(theta);
        Normal_zero[7] = -nz*sin(theta);
        Normal_zero[8] =  nx*cos(theta) + ny*sin(theta);

        // Generation of transformation matrix
        for (unsigned int kk=0; kk < 6*number_of_nodes; kk += 3)
        {
            for (unsigned int i=0; i<3; i++)
            {
                for(unsigned int j=0; j<3; j++)
                {
                    Rotation(i+kk,j+kk) = Normal_zero(3*j+i);
                }
            }
        }

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void TimoshenkoBeamElement::CalculateLocalMatrix(Matrix& LocalMatrix) const
    {
        KRATOS_TRY

        unsigned int number_of_nodes = GetGeometry().size();

        if(LocalMatrix.size1() != 6*number_of_nodes || LocalMatrix.size2() != 6*number_of_nodes)
        {
            LocalMatrix.resize(6*number_of_nodes, 6*number_of_nodes, false);
        }
        noalias(LocalMatrix) = ZeroMatrix(6*number_of_nodes, 6*number_of_nodes);

        //Initialization of local stiffness matrix
        const double Poisson = GetProperties()[POISSON_RATIO];
        const double Youngs  = GetProperties()[YOUNG_MODULUS];
        const double ShearModulus = Youngs / (2.0*(1.0 + Poisson));

        double const EA   =  mArea * Youngs;
        double const EIy  =  mInertia_y * Youngs;
        double const EIz  =  mInertia_z * Youngs;
        double const GJ   =  mInertia_x * ShearModulus;
        double const GAy  = mArea_y * ShearModulus;
        double const GAz  = mArea_z * ShearModulus;

        //this is the size of the elements stiffness matrix/force vector
        unsigned int mat_size = number_of_nodes * 6;

        //Initialize local variables
        Matrix B( 6, mat_size );
        Matrix C = zero_matrix<double>(6, 6);
        C(0, 0) = EA;
        C(1, 1) = GAy;
        C(2, 2) = GAz;
        C(3, 3) = GJ;
        C(4, 4) = EIy;
        C(5, 5) = EIz;

        Matrix DN_DX( number_of_nodes, 1 );

        //reading integration points, shape function values and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        Matrix DeltaPosition(number_of_nodes, 3);

        for ( unsigned int node = 0; node < number_of_nodes; ++node )
            noalias( row( DeltaPosition, node ) ) = this->GetGeometry()[node].Coordinates()
                            - this->GetGeometry()[node].GetInitialPosition();

        typename GeometryType::JacobiansType J0;
        J0 = this->GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            double DetJ0 = std::sqrt(MathUtils<double>::Det(Matrix(prod(trans(J0[PointNumber]), J0[PointNumber]))));

            noalias( DN_DX ) = DN_De[PointNumber] / DetJ0;

            //Initializing B_Operator at the current integration point
            CalculateBoperator( B, DN_DX, Ncontainer);

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight = integration_points[PointNumber].Weight();

            //calculate stiffness matrix
            noalias( LocalMatrix ) +=
                prod( trans(B), (IntToReferenceWeight * DetJ0) * Matrix(prod(C, B)));
        }

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void TimoshenkoBeamElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int mat_size = number_of_nodes * 6;

        if (rLeftHandSideMatrix.size1() != mat_size
            || rLeftHandSideMatrix.size2() != mat_size)
        {
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        }
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

        //Initialization of local stiffness matrix
        const double Poisson = GetProperties()[POISSON_RATIO];
        const double Youngs  = GetProperties()[YOUNG_MODULUS];
        const double ShearModulus = Youngs / (2.0*(1.0 + Poisson));

        double const EA   = mArea * Youngs;
        double const EIy  = mInertia_y * Youngs;
        double const EIz  = mInertia_z * Youngs;
        double const GJ   = mInertia_x * ShearModulus;
        double const GAy  = mArea_y * ShearModulus;
        double const GAz  = mArea_z * ShearModulus;

        //Initialize local variables
        Matrix B( 6, mat_size );
        Matrix C = zero_matrix<double>(6, 6);
        C(0, 0) = EA;
        C(1, 1) = GAy;
        C(2, 2) = GAz;
        C(3, 3) = GJ;
        C(4, 4) = EIy;
        C(5, 5) = EIz;

        Matrix DN_DX( number_of_nodes, 1 );
        Matrix LocalMatrix(mat_size, mat_size);
        Matrix Rotation(mat_size, mat_size);
        Matrix AuxMatrix(mat_size, mat_size);
        Vector TangentialVector(3);

        //reading integration points, shape function values and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        Matrix DeltaPosition(number_of_nodes, 3);

        for ( unsigned int node = 0; node < number_of_nodes; ++node )
            noalias( row( DeltaPosition, node ) ) = this->GetGeometry()[node].Coordinates()
                            - this->GetGeometry()[node].GetInitialPosition();

        typename GeometryType::JacobiansType J0;
        J0 = this->GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            double DetJ0 = std::sqrt(MathUtils<double>::Det(Matrix(prod(trans(J0[PointNumber]), J0[PointNumber]))));

            noalias( DN_DX ) = DN_De[PointNumber] / DetJ0;

            // Initializing B_Operator at the current integration point
            CalculateBoperator( B, DN_DX, Ncontainer);

            // calculating weights for integration on the reference configuration
            double IntToReferenceWeight = integration_points[PointNumber].Weight();

            // calculate local stiffness matrix
            noalias( LocalMatrix ) =
                prod( trans(B), (IntToReferenceWeight * DetJ0) * Matrix(prod(C, B)));

            // calculate tangential vector
            noalias(TangentialVector) = ZeroVector( 3 );
            for ( unsigned int n = 0; n < number_of_nodes; ++n )
            {
                TangentialVector[0] += GetGeometry()[n].X0() * DN_De[PointNumber]( n, 0 );
                TangentialVector[1] += GetGeometry()[n].Y0() * DN_De[PointNumber]( n, 0 );
                TangentialVector[2] += GetGeometry()[n].Z0() * DN_De[PointNumber]( n, 0 );
            }
            TangentialVector *= (1.0 / norm_2(TangentialVector));

            // calculate the transformation matrix
            CalculateTransformationMatrix(Rotation, TangentialVector);

            // contribute to the elemental stiffness matrix
            noalias(AuxMatrix) = prod(Rotation, LocalMatrix);
            noalias(rLeftHandSideMatrix) += prod(AuxMatrix, Matrix(trans(Rotation)));
        }

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void TimoshenkoBeamElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        unsigned int number_of_nodes = GetGeometry().size();

        //Initialization of variables
        Matrix LocalMatrix;
        Matrix Rotation;
        Matrix aux_matrix;
        Vector LocalBody;
        Vector CurrentDisplacement;

        CurrentDisplacement.resize(6*number_of_nodes);
        LocalMatrix.resize(6*number_of_nodes, 6*number_of_nodes, false);
        Rotation.resize(6*number_of_nodes, 6*number_of_nodes, false);
        aux_matrix.resize(6*number_of_nodes, 6*number_of_nodes, false);
        rLeftHandSideMatrix.resize(6*number_of_nodes, 6*number_of_nodes, false);
        rRightHandSideVector = ZeroVector(6*number_of_nodes);
        LocalBody = ZeroVector(6*number_of_nodes);

        CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

        // //Calculating RHS
        // CalculateExternalLoadVector(Rotation, LocalBody, rRightHandSideVector);

        for(unsigned int i = 0; i < number_of_nodes; ++i)
        {
            CurrentDisplacement(6*i    ) = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X) - mInitialDisp(i, 0);
            CurrentDisplacement(6*i + 1) = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y) - mInitialDisp(i, 1);
            CurrentDisplacement(6*i + 2) = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z) - mInitialDisp(i, 2);
            CurrentDisplacement(6*i + 3) = GetGeometry()[i].GetSolutionStepValue(ROTATION_X) - mInitialRot(i, 0);
            CurrentDisplacement(6*i + 4) = GetGeometry()[i].GetSolutionStepValue(ROTATION_Y) - mInitialRot(i, 1);
            CurrentDisplacement(6*i + 5) = GetGeometry()[i].GetSolutionStepValue(ROTATION_Z) - mInitialRot(i, 2);
        }

        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, CurrentDisplacement);

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void TimoshenkoBeamElement::CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX, const Matrix& Ncontainer ) const
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();

        noalias( B_Operator ) = ZeroMatrix( 6, 6*number_of_nodes );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            B_Operator( 0, i*6 ) = DN_DX( i, 0 );

            B_Operator( 1, i*6 + 1 ) = DN_DX( i, 0 );
            B_Operator( 1, i*6 + 5 ) = -Ncontainer(0, i);

            B_Operator( 2, i*6 + 2 ) = DN_DX( i, 0 );
            B_Operator( 2, i*6 + 4 ) = Ncontainer(0, i);

            B_Operator( 3, i*6 + 3 ) = DN_DX( i, 0 );

            B_Operator( 4, i*6 + 4 ) = DN_DX( i, 0 );

            B_Operator( 5, i*6 + 5 ) = DN_DX( i, 0 );
        }

        KRATOS_CATCH( "" )
    }

} // Namespace Kratos
