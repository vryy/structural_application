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
 *   Date:                $Date: 21 Nov 2023 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

// System includes


// External includes

// Project includes
#include "custom_elements/total_lagrangian_axisymmetric.h"
#include "custom_utilities/sd_math_utils.h"


namespace Kratos
{
    TotalLagrangianAxisymmetric::TotalLagrangianAxisymmetric( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseType( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
        //THIS IS THE DEFAULT CONSTRUCTOR
    }

    TotalLagrangianAxisymmetric::TotalLagrangianAxisymmetric( IndexType NewId,
            GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseType( NewId, pGeometry, pProperties )
    {}

    Element::Pointer TotalLagrangianAxisymmetric::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new TotalLagrangianAxisymmetric( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    Element::Pointer TotalLagrangianAxisymmetric::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new TotalLagrangianAxisymmetric( NewId, pGeom, pProperties ) );
    }

    TotalLagrangianAxisymmetric::~TotalLagrangianAxisymmetric()
    {
    }

    void TotalLagrangianAxisymmetric::CalculateF( Matrix& F, const Vector& N, const Matrix& DN_DX, const Matrix& CurrentDisp ) const
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = CurrentDisp.size1();

        double r = 0.0;
        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            r += N[i] * GetGeometry()[i].X0();
        }

        F.clear();

        for (unsigned int i = 0; i < 2; ++i)
        {
            for (unsigned int j = 0; j < 2; ++j)
            {
                for (unsigned int n = 0; n < number_of_nodes; ++n)
                {
                    F(i, j) += DN_DX(n, j) * CurrentDisp(n, i);
                }
            }
            F(i, i) += 1.0;
        }

        F(2, 2) = 1.0;
        for (unsigned int n = 0; n < number_of_nodes; ++n)
            F(2, 2) += N(n) / r * CurrentDisp(n, 0);

        KRATOS_CATCH( "" )
    }

    void TotalLagrangianAxisymmetric::CalculateB( Matrix& B, const Matrix& F, const Vector& N, const Matrix& DN_DX ) const
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();

        B.clear();

        double r = 0.0;
        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            r += N[i] * GetGeometry()[i].X0();
        }

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            const unsigned int index = i*2;
            B( 0, index     ) = F( 0, 0 ) * DN_DX( i, 0 );
            B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
            B( 1, index     ) = F( 0, 1 ) * DN_DX( i, 1 );
            B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
            B( 2, index     ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
            B( 2, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
            B( 3, index     ) = F( 2, 2 ) * N( i ) / r;
        }

        KRATOS_CATCH( "" )
    }

    void TotalLagrangianAxisymmetric::CalculateStrain( const Matrix& C, Vector& StrainVector ) const
    {
        KRATOS_TRY

        StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        StrainVector[2] = C( 0, 1 );

        StrainVector[3] = 0.5 * ( C( 2, 2 ) - 1.00 );

        KRATOS_CATCH( "" )
    }

    double TotalLagrangianAxisymmetric::GetIntegrationWeight( double Weight, const Vector& N ) const
    {
        // compute r
        double r = 0.0;
        for (std::size_t i = 0; i < GetGeometry().size(); ++i)
        {
            r += N(i) * GetGeometry()[i].X0();
        }

        return Weight * 2*SD_MathUtils<double>::Pi()*r;
    }

    void TotalLagrangianAxisymmetric::CalculateAndAddKg(
        MatrixType& K,
        const Vector& N,
        const Matrix& DN_DX,
        const Vector& StressVector,
        double Weight ) const
    {
        KRATOS_TRY

        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_nodes = GetGeometry().size();
        Matrix DN_DX_Axi(number_of_nodes, 3);

        double r = 0.0;
        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            r += N[i] * GetGeometry()[i].X0();
        }

        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            for (unsigned int j = 0; j < 2; ++j)
                DN_DX_Axi(i, j) = DN_DX(i, j);
            DN_DX_Axi(i, 2) = N(i)/r;
        }

        Matrix StressTensor(3, 3);
        SD_MathUtils<double>::StressVectorToTensor( StressVector, StressTensor );
        Matrix ReducedKg = prod( DN_DX_Axi, Weight * Matrix( prod( StressTensor, trans( DN_DX_Axi ) ) ) ); //to be optimized
        MathUtils<double>::ExpandAndAddReducedMatrix( K, ReducedKg, dimension );

        KRATOS_CATCH( "" )
    }

} // Namespace Kratos
