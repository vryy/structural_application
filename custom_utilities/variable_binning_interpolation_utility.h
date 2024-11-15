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
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 29 Mar 2017 $
*   Revision:            $Revision: 1.13 $
*
* ***********************************************************/

#if !defined(KRATOS_VARIABLE_BINNING_INTERPOLATION_UTILITY_INCLUDED )
#define  KRATOS_VARIABLE_BINNING_INTERPOLATION_UTILITY_INCLUDED

//System includes
#ifdef _OPENMP
#include <omp.h>
#endif

//External includes

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "utilities/progress.h"
#include "custom_utilities/variable_interpolation_utility.h"

namespace Kratos
{

/**
 * Utility to transfer the variables by interpolation.
 * It uses spatial binning to quickly search for corresponding element
 */
class VariableBinningInterpolationUtility : public VariableInterpolationUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( VariableBinningInterpolationUtility );

    typedef VariableInterpolationUtility BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef BaseType::PointType PointType;
    typedef BaseType::NodesContainerType NodesContainerType;
    typedef BaseType::ElementsContainerType ElementsContainerType;

    /**
     * Constructor.
     */
    VariableBinningInterpolationUtility(ElementsContainerType& pElements,
            const double& Dx, const double& Dy, const double& Dz)
    : BaseType(pElements), mDx(Dx), mDy(Dy), mDz(Dz)
    {
        this->Initialize(pElements);
        std::cout << "VariableBinningInterpolationUtility created" << std::endl;
    }

    /**
     * Destructor.
     */
    virtual ~VariableBinningInterpolationUtility()
    {
    }

    /**
     * Operations
     */

protected:

    /// Initialize the elements binning
    void Initialize( const ElementsContainerType& pElements ) final
    {
        std::cout << "Initialize the spatial binning" << std::endl;

#ifdef _OPENMP
        double start_init = omp_get_wtime();
#endif

        Kratos::progress_display show_progress( pElements.size() );
        std::vector<double> vmin(3);
        std::vector<double> vmax(3);
        for ( ElementsContainerType::const_iterator it = pElements.begin(); it != pElements.end(); ++it )
        {
            this->FindBoundingBox(vmin, vmax, it->GetGeometry());

            if (GetEchoLevel() > 4)
            {
                std::cout << "vmin: " << vmin[0] << " " << vmin[1] << " " << vmin[2] << std::endl;
                std::cout << "vmax: " << vmax[0] << " " << vmax[1] << " " << vmax[2] << std::endl;
            }

            int i_min = (int) floor(vmin[0] / mDx), i_max = (int) floor(vmax[0] / mDx);
            int j_min = (int) floor(vmin[1] / mDy), j_max = (int) floor(vmax[1] / mDy);
            int k_min = (int) floor(vmin[2] / mDz), k_max = (int) floor(vmax[2] / mDz);

            if (GetEchoLevel() > 4)
            {
                std::cout << " [" << i_min << "," << i_max << "]";
                std::cout << " [" << j_min << "," << j_max << "]";
                std::cout << " [" << k_min << "," << k_max << "]";
                std::cout << std::endl;
            }

            for (int i = i_min; i <= i_max; ++i)
            {
                for (int j = j_min; j <= j_max; ++j)
                {
                    for (int k = k_min; k <= k_max; ++k)
                    {
                        SpatialKey key(i, j, k);
                        mBinElements[key].insert(it->Id());

                        if (GetEchoLevel() > 4)
                            std::cout << "element " << it->Id() << " is added to the bin with key (" << i << "," << j << "," << k << ")" << std::endl;
                    }
                }
            }

            ++show_progress;
        }

        KRATOS_WATCH(mBinElements.size())

#ifdef _OPENMP
        double stop_init = omp_get_wtime();
        std::cout << "Initialize binning completed, time = " << (stop_init-start_init) << "s" << std::endl;
#endif
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************

    void FindPotentialPartners( const PointType& rSourcePoint, ElementsContainerType& pMasterElements ) const final
    {
        // get the containing elements from the bin
        int ix = (int) floor(rSourcePoint.X() / mDx);
        int iy = (int) floor(rSourcePoint.Y() / mDy);
        int iz = (int) floor(rSourcePoint.Z() / mDz);

        if (GetEchoLevel() > 4)
            std::cout << "point " << rSourcePoint << " has key (" << ix << "," << iy << "," << iz << ")" << std::endl;

        SpatialKey key(ix, iy, iz);
        std::map<SpatialKey, std::set<std::size_t> >::const_iterator it_bin_elements = mBinElements.find(key);

        if(it_bin_elements != mBinElements.end())
        {
            for(std::set<std::size_t>::const_iterator it = it_bin_elements->second.begin(); it != it_bin_elements->second.end(); ++it )
            {
                auto it_elem = BaseType::mpElements.find(*it);
                if (it_elem != BaseType::mpElements.end())
                    pMasterElements.push_back(*(it_elem.base()));
            }
        }
    }

private:

    struct SpatialKey
    {
        public:
            SpatialKey(int ix, int iy, int iz) : x(ix), y(iy), z(iz) {}
            bool operator<(const SpatialKey& rOther) const
            {
                if(x == rOther.x)
                {
                    if(y == rOther.y)
                    {
                        return z < rOther.z;
                    }
                    else
                        return y < rOther.y;
                }
                else
                    return x < rOther.x;
            }
            int kx() const {return x;}
            int ky() const {return y;}
            int kz() const {return z;}
        private:
            int x, y, z;
    };

    double mDx, mDy, mDz;
    std::map<SpatialKey, std::set<std::size_t> > mBinElements;

    void FindBoundingBox(std::vector<double>& vmin, std::vector<double>& vmax, GeometryType& rGeometry) const
    {
        vmin[0] = rGeometry[0].X(); vmin[1] = rGeometry[0].Y(); vmin[2] = rGeometry[0].Z();
        vmax[0] = rGeometry[0].X(); vmax[1] = rGeometry[0].Y(); vmax[2] = rGeometry[0].Z();
        for (std::size_t i = 1; i < rGeometry.size(); ++i)
        {
            if (rGeometry[i].X() < vmin[0]) vmin[0] = rGeometry[i].X();
            if (rGeometry[i].X() > vmax[0]) vmax[0] = rGeometry[i].X();
            if (rGeometry[i].Y() < vmin[1]) vmin[1] = rGeometry[i].Y();
            if (rGeometry[i].Y() > vmax[1]) vmax[1] = rGeometry[i].Y();
            if (rGeometry[i].Z() < vmin[2]) vmin[2] = rGeometry[i].Z();
            if (rGeometry[i].Z() > vmax[2]) vmax[2] = rGeometry[i].Z();
        }
    }

};//Class VariableBinningInterpolationUtility

}//namespace Kratos.

#endif /* KRATOS_VARIABLE_BINNING_INTERPOLATION_UTILITY_INCLUDED  defined */
