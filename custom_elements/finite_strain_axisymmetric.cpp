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
 *   Date:                $Date: 11 Oct 2023 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

// System includes


// External includes

// Project includes
#include "custom_elements/finite_strain_axisymmetric.h"
#include "custom_utilities/sd_math_utils.h"


namespace Kratos
{
    FiniteStrainAxisymmetric::FiniteStrainAxisymmetric( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseType( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
        //THIS IS THE DEFAULT CONSTRUCTOR
    }

    FiniteStrainAxisymmetric::FiniteStrainAxisymmetric( IndexType NewId,
            GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseType( NewId, pGeometry, pProperties )
    {}

    Element::Pointer FiniteStrainAxisymmetric::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new FiniteStrainAxisymmetric( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    Element::Pointer FiniteStrainAxisymmetric::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new FiniteStrainAxisymmetric( NewId, pGeom, pProperties ) );
    }

    FiniteStrainAxisymmetric::~FiniteStrainAxisymmetric()
    {
    }

    void FiniteStrainAxisymmetric::CalculateF( Matrix& F, const Matrix& G_Operator, const Matrix& CurrentDisp ) const
    {
        SD_MathUtils<double>::CalculateFaxi(F, G_Operator, CurrentDisp);
    }

    void FiniteStrainAxisymmetric::CalculateB( Matrix& B_Operator, const Vector& N, const Matrix& DN_DX, const Matrix& CurrentDisp ) const
    {
        KRATOS_TRY

        SD_MathUtils<double>::CalculateBaxi( B_Operator, GetGeometry(), N, DN_DX, CurrentDisp );

        KRATOS_CATCH( "" )
    }

    void FiniteStrainAxisymmetric::CalculateG( Matrix& G_Operator, const VectorType& N, const Matrix& DN_DX ) const
    {
        KRATOS_TRY

        SD_MathUtils<double>::CalculateGaxi( G_Operator, GetGeometry(), N, DN_DX );

        KRATOS_CATCH( "" )
    }

    void FiniteStrainAxisymmetric::CalculateG( Matrix& G_Operator, const VectorType& N, const Matrix& DN_DX, const Matrix& CurrentDisp ) const
    {
        KRATOS_TRY

        SD_MathUtils<double>::CalculateGaxi( G_Operator, GetGeometry(), N, DN_DX, CurrentDisp );

        KRATOS_CATCH( "" )
    }

    void FiniteStrainAxisymmetric::CalculateStrain( const Matrix& B, Vector& StrainVector ) const
    {
        KRATOS_TRY

        Matrix InvB(3, 3);
        double DetB;
        MathUtils<double>::InvertMatrix(B, InvB, DetB);

        StrainVector[0] = 0.5 * ( 1.00 - InvB( 0, 0 ) );

        StrainVector[1] = 0.5 * ( 1.00 - InvB( 1, 1 ) );

        StrainVector[2] = -InvB( 0, 1 );

        StrainVector[3] = 0.5 * ( 1.00 - InvB( 2, 2 ) );

        KRATOS_CATCH( "" )
    }

    double FiniteStrainAxisymmetric::GetIntegrationWeight( double Weight, const VectorType& N, const Matrix& CurrentDisp ) const
    {
        // compute r
        double r = 0.0;
        for (std::size_t i = 0; i < GetGeometry().size(); ++i)
        {
            r += N(i) * (GetGeometry()[i].X0() + CurrentDisp(i, 0));
        }

        return Weight * 2*SD_MathUtils<double>::Pi()*r;
    }

} // Namespace Kratos
