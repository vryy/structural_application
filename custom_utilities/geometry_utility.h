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
*   Last Modified by:    $Author: Hoang-Giang Bui $
*   Date:                $Date: 19 Sep 2023 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_STRUCTURAL_GEOMETRY_UTILITY_H_INCLUDED )
#define  KRATOS_STRUCTURAL_GEOMETRY_UTILITY_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "geometries/geometry.h"

namespace Kratos
{

/**
 * Utility for geometric routines
 */
class GeometryUtility
{
public:

    /**
     * Type Definitions
     */

    /**
     * Counted pointer of ContactUtility
     */
    KRATOS_CLASS_POINTER_DEFINITION( GeometryUtility );

    /**
     * Life Cycle
     */

    /**
     * Constructor.
     */
    GeometryUtility( )
    {
    }

    /**
     * Destructor.
     */
    virtual ~GeometryUtility() {}

    /**
     * Operations
     */

    /// Obtain the origin in the local coordinates
    template<typename CoordinatesArrayType>
    static void ComputeOrigin(GeometryData::KratosGeometryFamily geom_family, CoordinatesArrayType& rCoordinates)
    {
        if ( geom_family == GeometryData::KratosGeometryFamily::Kratos_Hexahedra
          || geom_family == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral
          || geom_family == GeometryData::KratosGeometryFamily::Kratos_Linear )
        {
            rCoordinates[0] = 0.0;
            rCoordinates[1] = 0.0;
            rCoordinates[2] = 0.0;
        }
        else if ( geom_family == GeometryData::KratosGeometryFamily::Kratos_Tetrahedra )
        {
            rCoordinates[0] = 0.25;
            rCoordinates[1] = 0.25;
            rCoordinates[2] = 0.25;
        }
        else if ( geom_family == GeometryData::KratosGeometryFamily::Kratos_Triangle )
        {
            rCoordinates[0] = 1.0 / 3;
            rCoordinates[1] = 1.0 / 3;
            rCoordinates[2] = 0.0;
        }
        else if ( geom_family == GeometryData::KratosGeometryFamily::Kratos_NURBS )
        {
            rCoordinates[0] = 0.5;
            rCoordinates[1] = 0.5;
            rCoordinates[2] = 0.5;
        }
        else
        {
            KRATOS_ERROR << "Geometrily family of type " << static_cast<int>(geom_family)
                         << " is not yet supported";
        }
    }

};//class GeometryUtility
}  /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_GEOMETRY_UTILITY_H_INCLUDED  defined */
