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
*   Date:                $Date: 9 Apr 2021 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_VARIABLE_BVH_INTERPOLATION_UTILITY_INCLUDED )
#define  KRATOS_VARIABLE_BVH_INTERPOLATION_UTILITY_INCLUDED

//System includes
#ifdef _OPENMP
#include <omp.h>
#endif

//External includes

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "spatial_containers/bounding_volume_tree.h"
#include "custom_utilities/variable_interpolation_utility.h"

namespace Kratos
{

/**
 * Utility to transfer the variables by interpolation.
 * It uses Bounding Volume Hierarchy to quickly search for corresponding element
 */
class VariableBVHInterpolationUtility : public VariableInterpolationUtility
{
public:

    typedef VariableInterpolationUtility BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef BaseType::PointType PointType;
    typedef BaseType::NodesContainerType NodesContainerType;
    typedef BaseType::ElementsContainerType ElementsContainerType;
    typedef BoundingVolumeTree<0, ElementsContainerType> BoundingVolumeTreeType;

    /**
     * Constructor.
     */
    VariableBVHInterpolationUtility(const ElementsContainerType& pElements, const int bv_type)
    : BaseType(pElements)
    {
        mpBVTree = typename BoundingVolumeTreeType::Pointer(new BoundingVolumeTreeType(bv_type));
        this->Initialize(pElements);
        std::cout << "VariableBVHInterpolationUtility created, BVH depth = " << mpBVTree->Depth() << std::endl;
    }

    /**
     * Destructor.
     */
    ~VariableBVHInterpolationUtility() override
    {
    }

    /**
     * Operations
     */

protected:

    /// Initialize the elements binning
    void Initialize( const ElementsContainerType& pElements ) final
    {
        if (GetEchoLevel() > 0)
            std::cout << "Initialize the Bounding Volume Hierarchy" << std::endl;

#ifdef _OPENMP
        double start_init = omp_get_wtime();
#endif

        SimpleBoundingVolumePartitioner<0, ElementsContainerType> partitioner;
        mpBVTree->BuildTreeTopDown(pElements, partitioner);

#ifdef _OPENMP
        double stop_init = omp_get_wtime();
        if (GetEchoLevel() > 0)
            std::cout << "Initialize BVH completed, time = " << (stop_init-start_init) << "s" << std::endl;
#endif
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************

    void FindPotentialPartners( const PointType& rSourcePoint, ElementsContainerType& pMasterElements ) const final
    {
        const double TOL = 1.0e-8;
        std::set<std::size_t> elem_ids;
        mpBVTree->GetContainingGeometries(elem_ids, rSourcePoint, TOL);

        for( auto it = elem_ids.begin(); it != elem_ids.end(); ++it )
        {
            auto it_elem = BaseType::mpElements.find(*it);
            if (it_elem != BaseType::mpElements.end())
                pMasterElements.push_back(*(it_elem.base()));
        }
    }

private:

    typename BoundingVolumeTreeType::Pointer mpBVTree;

}; // Class VariableBVHInterpolationUtility

} // namespace Kratos.

#endif /* KRATOS_VARIABLE_BVH_INTERPOLATION_UTILITY_INCLUDED  defined */
