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
*   Date:                $Date: 31 May 2022 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_SURFACE_UTILITY_H_INCLUDED )
#define  KRATOS_SURFACE_UTILITY_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "structural_application_variables.h"

namespace Kratos
{

/**
 * auxiliary functions for calculating on surface
 */
class SurfaceUtility
{
public:

    /**
     * Type Definitions
     */

    /**
     * Counted pointer of ContactUtility
     */
    KRATOS_CLASS_POINTER_DEFINITION( SurfaceUtility );

    typedef std::size_t IndexType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    typedef Condition::GeometryType GeometryType;

    typedef GeometryType::IntegrationPointsArrayType            IntegrationPointsArrayType;

    typedef GeometryType::ShapeFunctionsSecondDerivativesType   ShapeFunctionsSecondDerivativesType;

    typedef Condition::IntegrationMethod                        IntegrationMethod;

    /**
     * Life Cycle
     */

    /**
     * Constructor.
     */
    SurfaceUtility( int echo_level )
    {
        mEchoLevel = echo_level;
    }

    /**
     * Destructor.
     */
    virtual ~SurfaceUtility() {}

    /**
     * Operations
     */

    /// Integrate the strain on the surface:
    ///     epsilon = int_S 1/2 * (u \otimes n + n \otimes u) dA
    /// If TFrame == 0, the normal vector is computed in the reference configuration, as well as the Jacobian
    /// If TFrame == 1, the integration is performed in the current configuration
    /// On output, strain vector of the form [e_xx e_yy e_zz e_xy e_yz e_xz] is returned
    template<int TDim, int TFrame>
    Vector CalculateStrain(const ConditionsContainerType& rConditions) const
    {
        Vector N(3), T1(3), T2(3), u(3);
        Matrix epsilon(3, 3);
        double DetJ;
        double Area = 0.0;

        epsilon.clear();

        for (auto it = rConditions.begin(); it != rConditions.end(); ++it)
        {
            const GeometryType& rGeometry = it->GetGeometry();

            const IntegrationMethod& ThisIntegrationMethod = rGeometry.GetDefaultIntegrationMethod();

            const IntegrationPointsArrayType& integration_points = rGeometry.IntegrationPoints( ThisIntegrationMethod );

            const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( ThisIntegrationMethod );

            const GeometryType::ShapeFunctionsGradientsType& DN_De = rGeometry.ShapeFunctionsLocalGradients( ThisIntegrationMethod );

            for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
            {
                const Matrix& DN = DN_De[PointNumber];

                // compute the tangential and normal vectors
                if (TDim == 2)
                {
                    T1.clear();

                    for(std::size_t n = 0; n < rGeometry.PointsNumber(); ++n)
                    {
                        if (TFrame == 0)
                        {
                            // contribution to tangential vectors
                            T1[0] += rGeometry[n].X0() * DN(n, 0);
                            T1[1] += rGeometry[n].Y0() * DN(n, 0);
                        }
                        else if (TFrame == 1)
                        {
                            // contribution to tangential vectors
                            T1[0] += (rGeometry[n].X0() + rGeometry[n].GetSolutionStepValue(DISPLACEMENT_X)) * DN(n, 0);
                            T1[1] += (rGeometry[n].Y0() + rGeometry[n].GetSolutionStepValue(DISPLACEMENT_Y)) * DN(n, 0);
                        }
                    }

                    N.clear();
                    N[0] = T1[1];
                    N[1] = -T1[0]; // be careful here. It's only correct for GiD mesh when normal in 2D mesh always pointing inside

                    DetJ = norm_2(N);
                    N *= (1.0 / DetJ);
                }
                else if (TDim == 3)
                {
                    T1.clear();
                    T2.clear();

                    for(std::size_t n = 0; n < rGeometry.PointsNumber(); ++n)
                    {
                        if (TFrame == 0)
                        {
                            // contribution to tangential vectors
                            T1[0] += rGeometry[n].X0() * DN(n, 0);
                            T1[1] += rGeometry[n].Y0() * DN(n, 0);
                            T1[2] += rGeometry[n].Z0() * DN(n, 0);

                            T2[0] += rGeometry[n].X0() * DN(n, 1);
                            T2[1] += rGeometry[n].Y0() * DN(n, 1);
                            T2[2] += rGeometry[n].Z0() * DN(n, 1);
                        }
                        else if (TFrame == 1)
                        {
                            // contribution to tangential vectors
                            T1[0] += (rGeometry[n].X0() + rGeometry[n].GetSolutionStepValue(DISPLACEMENT_X)) * DN(n, 0);
                            T1[1] += (rGeometry[n].Y0() + rGeometry[n].GetSolutionStepValue(DISPLACEMENT_Y)) * DN(n, 0);
                            T1[2] += (rGeometry[n].Z0() + rGeometry[n].GetSolutionStepValue(DISPLACEMENT_Z)) * DN(n, 0);

                            T2[0] += (rGeometry[n].X0() + rGeometry[n].GetSolutionStepValue(DISPLACEMENT_X)) * DN(n, 1);
                            T2[1] += (rGeometry[n].Y0() + rGeometry[n].GetSolutionStepValue(DISPLACEMENT_Y)) * DN(n, 1);
                            T2[2] += (rGeometry[n].Z0() + rGeometry[n].GetSolutionStepValue(DISPLACEMENT_Z)) * DN(n, 1);
                        }
                    }

                    N[0] = T1[1] * T2[2] - T1[2] * T2[1];
                    N[1] = T1[2] * T2[0] - T1[0] * T2[2];
                    N[2] = T1[0] * T2[1] - T1[1] * T2[0];

                    DetJ = norm_2(N);
                    N *= (1.0 / DetJ);
                }

                if (mEchoLevel > 1)
                {
                    Vector P(3);
                    P.clear();
                    for(std::size_t n = 0; n < rGeometry.PointsNumber(); ++n)
                    {
                        P[0] += rGeometry[n].X0() * Ncontainer(PointNumber, n);
                        P[1] += rGeometry[n].Y0() * Ncontainer(PointNumber, n);
                        P[2] += rGeometry[n].Z0() * Ncontainer(PointNumber, n);
                    }

                    std::cout << "Normal vector at point (" << P[0] << ", " << P[1] << ", " << P[2] << "): " << N << std::endl;
                }

                Area += DetJ * integration_points[PointNumber].Weight();

                // interpolate the displacement
                u.clear();
                for(std::size_t n = 0; n < rGeometry.PointsNumber(); ++n)
                {
                    u[0] += rGeometry[n].GetSolutionStepValue(DISPLACEMENT_X) * Ncontainer(PointNumber, n);
                    u[1] += rGeometry[n].GetSolutionStepValue(DISPLACEMENT_Y) * Ncontainer(PointNumber, n);
                    u[2] += rGeometry[n].GetSolutionStepValue(DISPLACEMENT_Z) * Ncontainer(PointNumber, n);
                }

                // contribute to the strain tensor
                noalias(epsilon) += 0.5 * (outer_prod(N, u) + outer_prod(u, N)) * DetJ * integration_points[PointNumber].Weight();
            }
        }

        if (mEchoLevel > 0)
        {
            KRATOS_WATCH(epsilon)
            KRATOS_WATCH(Area)
        }

        // epsilon /= Area;

        Vector strain_vector(6);
        strain_vector[0] = epsilon(0, 0);
        strain_vector[1] = epsilon(1, 1);
        strain_vector[2] = epsilon(2, 2);
        strain_vector[3] = 2.0*epsilon(0, 1);
        strain_vector[4] = 2.0*epsilon(1, 2);
        strain_vector[5] = 2.0*epsilon(0, 2);

        return strain_vector;
    }

private:

    int mEchoLevel;

};//class SurfaceUtility
}  /* namespace Kratos.*/

#endif /* KRATOS_SURFACE_UTILITY_H_INCLUDED  defined */
