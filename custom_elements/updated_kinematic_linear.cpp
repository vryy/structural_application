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
#include "custom_elements/updated_kinematic_linear.h"

namespace Kratos
{
    UpdatedKinematicLinear::UpdatedKinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseType( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
        //THIS IS THE DEFAULT CONSTRUCTOR
    }

    UpdatedKinematicLinear::UpdatedKinematicLinear( IndexType NewId,
            GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseType( NewId, pGeometry, pProperties )
    {}

    Element::Pointer UpdatedKinematicLinear::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new UpdatedKinematicLinear( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    Element::Pointer UpdatedKinematicLinear::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new UpdatedKinematicLinear( NewId, pGeom, pProperties ) );
    }

    UpdatedKinematicLinear::~UpdatedKinematicLinear()
    {
    }

    void UpdatedKinematicLinear::CalculateJacobian( GeometryType::JacobiansType& J, const IntegrationMethod ThisIntegrationMethod ) const
    {
        J = GetGeometry().Jacobian( J, ThisIntegrationMethod );
    }

} // Namespace Kratos
