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

#if !defined(KRATOS_VARIABLE_ADVANCED_TRANSFER_UTILITY_INCLUDED )
#define  KRATOS_VARIABLE_ADVANCED_TRANSFER_UTILITY_INCLUDED

//System includes
#ifdef _OPENMP
#include <omp.h>
#endif

//External includes
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp"
#include "boost/progress.hpp"

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "spaces/ublas_space.h"
#include "structural_application.h"

namespace Kratos
{

/**
 * Utility to transfer the variables using L2-projection.
 * In this utility, a set of elements, e.g. a sub-set of elements in the model_part, have to be given. The variables can then be transferred to the nodes of the elements.
 * Subsequently the nodal variables can be transferred to other mesh by using interpolation.
 */
class VariableAdvancedTransferUtility
{
public:
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef LinearSolver<SparseSpaceType, DenseSpaceType> LinearSolverType;
    typedef Element::GeometryType GeometryType;
    typedef GeometryType::PointType NodeType;
    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef NodeType::PointType PointType;
    typedef ModelPart::NodesContainerType NodesContainerType;
    typedef ModelPart::ElementsContainerType ElementsContainerType;

    /**
     * Constructor.
     */
    VariableAdvancedTransferUtility(ModelPart::ElementsContainerType& pElements,
            const double& Dx, const double& Dy, const double& Dz)
    : mpElements(pElements), mEchoLevel(-1), mDx(Dx), mDy(Dy), mDz(Dz)
    {
        this->InitializeBinning(mpElements);
        std::cout << "VariableAdvancedTransferUtility created" << std::endl;
    }

    /**
     * Destructor.
     */
    virtual ~VariableAdvancedTransferUtility()
    {
    }

    /**
     * Access
     */

    void SetEchoLevel(const int& Level)
    {
        mEchoLevel = Level;
    }

    const int& GetEchoLevel() const
    {
        return mEchoLevel;
    }

    /// Initialize the elements binning
    void InitializeBinning(ElementsContainerType& pElements)
    {
        for ( ElementsContainerType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            std::pair<std::vector<double>, std::vector<double> > box = FindBoundingBox((*it)->GetGeometry());

            const std::vector<double>& vmin = box.first;
            const std::vector<double>& vmax = box.second;

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
                        mBinElements[key].insert((*it)->Id());

                        if (GetEchoLevel() > 4)
                            std::cout << "element " << (*it)->Id() << " is added to the bin with key (" << i << "," << j << "," << k << ")" << std::endl;
                    }
                }
            }
        }
        KRATOS_WATCH(mBinElements.size())
    }

    /// Transfer the variable at node from source mesh to target mesh
    /// Using a simple scheme, where the node in target mesh is located in source mesh, and then the value is determined
    template<class TVariableType>
    void TransferVariablesFromNodeToNode(NodesContainerType& rTargetNodes, TVariableType& rThisVariable)
    {
        for(NodesContainerType::iterator it = rTargetNodes.begin(); it != rTargetNodes.end(); ++it)
        {
            //Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
            Point<3> sourceLocalPoint;
            Element::Pointer sourceElement;

            bool found = this->SearchPartnerWithBin( *it, mpElements, sourceElement, sourceLocalPoint );

            if (found)
            {
                Vector shape_functions_values;
                shape_functions_values = sourceElement->GetGeometry().ShapeFunctionsValues(shape_functions_values, sourceLocalPoint);

                it->GetSolutionStepValue(rThisVariable) =
                    shape_functions_values[0] * sourceElement->GetGeometry()[0].GetSolutionStepValue(rThisVariable);

                for(unsigned int i = 1; i < sourceElement->GetGeometry().size(); ++i)
                {
                    it->GetSolutionStepValue(rThisVariable) +=
                        shape_functions_values[i] * sourceElement->GetGeometry()[i].GetSolutionStepValue(rThisVariable);
                }
            }
            else
            {
                std::cout<<"###### NO PARTNER FOUND IN OLD MESH : TransferVariablesBetweenMeshes(..." << rThisVariable.Name() << "...)#####"<<std::endl;
                continue;
            }
        }
    }

protected:

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    /// Find an element in pMasterElements contains rSourcePoint and assign it to pTargetElement.
    /// The rLocalTargetPoint is the local point in pTargetElement of rSourcePoint
    /// REMARK: we should disable the move mesh flag if we want to search in the reference configuration
    bool SearchPartnerWithBin( const PointType& rSourcePoint, ModelPart::ElementsContainerType& pMasterElements,
            Element::Pointer& pTargetElement, PointType& rLocalTargetPoint ) const
    {
        ModelPart::ElementsContainerType pMasterElementsCandidates;

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
                std::size_t elem_id = *it;
                Element::GeometryType& r_geom = pMasterElements[elem_id].GetGeometry();

                bool is_inside = r_geom.IsInside( rSourcePoint, rLocalTargetPoint );
                if( is_inside )
                {
                    pTargetElement = pMasterElements(elem_id);
                    return true;
                }
            }
        }

        if (GetEchoLevel() > 4)
            std::cout << " !!!! WARNING: NO ELEMENT FOUND TO CONTAIN " << rSourcePoint << " !!!! " << std::endl;

        return false;
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

    std::pair<std::vector<double>, std::vector<double> > FindBoundingBox(GeometryType& rGeometry) const
    {
        std::vector<double> vmin(3);
        std::vector<double> vmax(3);

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

        return std::make_pair(vmin, vmax);
    }

    ElementsContainerType mpElements;
    int mEchoLevel;
    double mDx, mDy, mDz;
    std::map<SpatialKey, std::set<std::size_t> > mBinElements;

};//Class VariableAdvancedTransferUtility

}//namespace Kratos.

#endif /* KRATOS_VARIABLE_ADVANCED_TRANSFER_UTILITY_INCLUDED  defined */
