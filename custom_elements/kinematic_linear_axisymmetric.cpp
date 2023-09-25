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
 *   Date:                $Date: 14 Feb 2022 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

// System includes


// External includes

// Project includes
#include "custom_elements/kinematic_linear_axisymmetric.h"
#include "custom_utilities/sd_math_utils.h"

//#define ENABLE_DEBUG_CONSTITUTIVE_LAW

namespace Kratos
{
    KinematicLinearAxisymmetric::KinematicLinearAxisymmetric( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseType( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
        //THIS IS THE DEFAULT CONSTRUCTOR
    }

    KinematicLinearAxisymmetric::KinematicLinearAxisymmetric( IndexType NewId,
            GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseType( NewId, pGeometry, pProperties )
    {}

    Element::Pointer KinematicLinearAxisymmetric::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new KinematicLinearAxisymmetric( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    Element::Pointer KinematicLinearAxisymmetric::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new KinematicLinearAxisymmetric( NewId, pGeom, pProperties ) );
    }

    KinematicLinearAxisymmetric::~KinematicLinearAxisymmetric()
    {
    }

    void KinematicLinearAxisymmetric::CalculateBoperator( Matrix& B_Operator, const Vector& N, const Matrix& DN_DX ) const
    {
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = this->GetStrainSize(dim);

        noalias( B_Operator ) = ZeroMatrix( strain_size, number_of_nodes * dim );

        double r = 0.0;

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            r += N[i] * GetGeometry()[i].X0();
        }

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            B_Operator( 0, i*2     ) = DN_DX( i, 0 );
            B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 );
            B_Operator( 2, i*2     ) = DN_DX( i, 1 );
            B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
            B_Operator( 3, i*2     ) = N( i ) / r;
        }
    }

    void KinematicLinearAxisymmetric::CalculateBBaroperator( Matrix& B_Operator, const Matrix& DN_DX, const Matrix& Bdil_bar ) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
    }

    double KinematicLinearAxisymmetric::GetIntegrationWeight( const GeometryType::IntegrationPointsArrayType& integration_points,
            const unsigned int& PointNumber, const Matrix& Ncontainer ) const
    {
        // compute r
        double r = 0.0;
        for (std::size_t i = 0; i < GetGeometry().size(); ++i)
        {
            r += Ncontainer(PointNumber, i) * GetGeometry()[i].X0();
        }

        return integration_points[PointNumber].Weight() * 2*SD_MathUtils<double>::Pi()*r;
    }

} // Namespace Kratos

#ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
#undef ENABLE_DEBUG_CONSTITUTIVE_LAW
#endif
