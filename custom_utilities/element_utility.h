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
*   Date:                $Date: 5 June 2025 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_ELEMENT_UTILITY_H_INCLUDED )
#define  KRATOS_ELEMENT_UTILITY_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "includes/element.h"

namespace Kratos
{

/**
 * Utility for elemental routines
 */
class ElementUtility
{
public:

    /**
     * Type Definitions
     */

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /**
     * Counted pointer of ContactUtility
     */
    KRATOS_CLASS_POINTER_DEFINITION( ElementUtility );

    /**
     * Life Cycle
     */

    /**
     * Constructor.
     */
    ElementUtility( )
    {
    }

    /**
     * Destructor.
     */
    virtual ~ElementUtility() {}

    /**
     * Operations
     */

    /// Interpolate on the element with the specified integration method
    template<typename TEntityType, typename TDataType>
    static void Interpolate(std::vector<TDataType>& rValues, TEntityType& rElement,
            const Variable<TDataType>& rVariable, const IntegrationMethod ThisIntegrationMethod)
    {
        typedef typename TEntityType::GeometryType GeometryType;

        GeometryType& rGeometry = rElement.GetGeometry();

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        rGeometry.Initialize(ThisIntegrationMethod);
        #endif

        const typename GeometryType::IntegrationPointsArrayType& integration_points =
                rGeometry.IntegrationPoints( ThisIntegrationMethod );

        const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( ThisIntegrationMethod );

        if (rValues.size() != integration_points.size())
            rValues.resize(integration_points.size());

        for (std::size_t point = 0; point < integration_points.size(); ++point)
        {
            TDataType value = 0;
            for(std::size_t i = 0; i < rGeometry.size(); ++i)
            {
                value += Ncontainer(point, i) * rGeometry[i].GetSolutionStepValue(rVariable);
            }
            rValues[point] = value;
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        rGeometry.Clean();
        #endif
    }

    /// Interpolate on the element with the specified integration method
    template<typename TEntityType, typename TDataType>
    static void Interpolate(std::vector<array_1d<TDataType, 3> >& rValues, TEntityType& rElement,
            const Variable<array_1d<TDataType, 3> >& rVariable, const IntegrationMethod ThisIntegrationMethod)
    {
        typedef typename TEntityType::GeometryType GeometryType;

        GeometryType& rGeometry = rElement.GetGeometry();

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        rGeometry.Initialize(ThisIntegrationMethod);
        #endif

        const typename GeometryType::IntegrationPointsArrayType& integration_points =
                rGeometry.IntegrationPoints( ThisIntegrationMethod );

        const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( ThisIntegrationMethod );

        if (rValues.size() != integration_points.size())
            rValues.resize(integration_points.size());

        for (std::size_t point = 0; point < integration_points.size(); ++point)
        {
            array_1d<TDataType, 3> value;
            for (unsigned int i = 0; i < 3; ++i) value[i] = 0;
            for(std::size_t i = 0; i < rGeometry.size(); ++i)
            {
                noalias(value) += Ncontainer(point, i) * rGeometry[i].GetSolutionStepValue(rVariable);
            }
            rValues[point] = value;
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        rGeometry.Clean();
        #endif
    }

    /// Interpolate on the element with the specified integration method
    template<typename TEntityType, typename TVectorType>
    static void Interpolate(std::vector<TVectorType>& rValues, TEntityType& rElement,
            const Variable<TVectorType>& rVariable, const IntegrationMethod ThisIntegrationMethod,
            const unsigned int ncomponents)
    {
        typedef typename TEntityType::GeometryType GeometryType;

        GeometryType& rGeometry = rElement.GetGeometry();

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        rGeometry.Initialize(ThisIntegrationMethod);
        #endif

        const typename GeometryType::IntegrationPointsArrayType& integration_points =
                rGeometry.IntegrationPoints( ThisIntegrationMethod );

        const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( ThisIntegrationMethod );

        if (rValues.size() != integration_points.size())
            rValues.resize(integration_points.size());

        for (std::size_t point = 0; point < integration_points.size(); ++point)
        {
            TVectorType value(ncomponents);
            for (unsigned int i = 0; i < ncomponents; ++i) value[i] = 0;
            for(std::size_t i = 0; i < rGeometry.size(); ++i)
            {
                noalias(value) += Ncontainer(point, i) * rGeometry[i].GetSolutionStepValue(rVariable);
            }
            rValues[point] = value;
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        rGeometry.Clean();
        #endif
    }

}; // class ElementUtility

}  /* namespace Kratos.*/

#endif /* KRATOS_ELEMENT_UTILITY_H_INCLUDED  defined */
