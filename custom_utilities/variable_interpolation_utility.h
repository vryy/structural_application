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

#if !defined(KRATOS_VARIABLE_INTERPOLATION_UTILITY_INCLUDED )
#define  KRATOS_VARIABLE_INTERPOLATION_UTILITY_INCLUDED

//System includes
#ifdef _OPENMP
#include <omp.h>
#endif

//External includes

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "utilities/openmp_utils.h"
#include "utilities/progress.h"
#include "custom_utilities/variable_utility.h"

namespace Kratos
{

/**
 * Utility to transfer the variables with the efficient search functionality
 */
class VariableInterpolationUtility : public VariableUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(VariableInterpolationUtility);

    typedef VariableUtility BaseType;
    typedef Element::GeometryType GeometryType;
    typedef GeometryType::PointType NodeType;
    typedef NodeType::PointType PointType;
    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef BaseType::NodesContainerType NodesContainerType;
    typedef BaseType::ElementsContainerType ElementsContainerType;

    /**
     * Constructor.
     */
    VariableInterpolationUtility(const ElementsContainerType& pElements)
    : BaseType(pElements)
    {
        std::cout << "VariableInterpolationUtility created" << std::endl;
    }

    VariableInterpolationUtility(const ElementsContainerType& pElements, const int EchoLevel)
    : BaseType(pElements, EchoLevel)
    {
        if (GetEchoLevel() > 0)
            std::cout << "VariableInterpolationUtility created" << std::endl;
    }

    /**
     * Destructor.
     */
    ~VariableInterpolationUtility() override
    {
    }

    /**
     * Operations
     */

    /// Get the elements of which the BV contains the point
    ElementsContainerType FindPotentialPartners( const PointType& rSourcePoint ) const
    {
        ElementsContainerType pMasterElements;
        this->FindPotentialPartners( rSourcePoint, pMasterElements );
        return pMasterElements;
    }

    /// Get the element containing the point
    Element::Pointer SearchPartner( const PointType& rSourcePoint, ElementsContainerType& pMasterElements ) const
    {
        Element::Pointer pElement;

        PointType localPoint;
        this->SearchPartner( rSourcePoint, pMasterElements, pElement, localPoint);

        return pElement;
    }

    /// Transfer the double variable to node of the target model_part
    void TransferVariablesToNodes(ModelPart& rTarget, const Variable<double>& rThisVariable)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << rTarget.Name() << ", Variable<double> " << rThisVariable.Name() << std::endl;

        TransferVariablesToNodes(rTarget.Nodes(), rThisVariable);
    }

    /// Transfer the double variable to node of the target node mesh
    void TransferVariablesToNodes(NodesContainerType& rTargetNodes, const Variable<double>& rThisVariable)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << " Variable<double> " << rThisVariable.Name() << std::endl;

        TransferVariablesToNodesImpl<DoubleVariableInitializer>(rTargetNodes, rThisVariable);
    }

    /// Transfer the double variable to node of the target model_part
    void TransferVariablesToNodes(ModelPart& rTarget, const Variable<array_1d<double, 3> >& rThisVariable)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << rTarget.Name() << ", Variable<array_1d<double, 3> > " << rThisVariable.Name() << std::endl;

        TransferVariablesToNodes(rTarget.Nodes(), rThisVariable);
    }

    /// Transfer the double variable to node of the target node mesh
    void TransferVariablesToNodes(NodesContainerType& rTargetNodes, const Variable<array_1d<double, 3> >& rThisVariable)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << " Variable<array_1d<double, 3> > " << rThisVariable.Name() << std::endl;

        TransferVariablesToNodesImpl<Array1DVariableInitializer>(rTargetNodes, rThisVariable);
    }

    /// Transfer the double variable to node of the target model_part
    void TransferVariablesToNodes(ModelPart& rTarget, const Variable<Vector>& rThisVariable, const std::size_t& ncomponents)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << rTarget.Name() << ", Variable<Vector> " << rThisVariable.Name() << std::endl;

        TransferVariablesToNodes(rTarget.Nodes(), rThisVariable, ncomponents);
    }

    /// Transfer the double variable to node of the target node mesh
    void TransferVariablesToNodes(NodesContainerType& rTargetNodes, const Variable<Vector>& rThisVariable, const std::size_t& ncomponents)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << " Variable<Vector> " << rThisVariable.Name() << std::endl;

        if (ncomponents == 3)
            TransferVariablesToNodesImpl<VectorVariableInitializer<3> >(rTargetNodes, rThisVariable);
        else if (ncomponents == 6)
            TransferVariablesToNodesImpl<VectorVariableInitializer<6> >(rTargetNodes, rThisVariable);
    }

    /// Transfer the double variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(ModelPart& rTarget, const Variable<double>& rThisVariable)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << rTarget.Name() << ", Variable<double> " << rThisVariable.Name() << std::endl;

        TransferVariablesToGaussPoints(rTarget.Elements(), rThisVariable, rTarget.GetProcessInfo());
    }

    /// Transfer the double variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(ElementsContainerType& TargetMeshElementsArray, const Variable<double>& rThisVariable, const ProcessInfo& CurrentProcessInfo)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << " Variable<double> " << rThisVariable.Name() << std::endl;

        TransferVariablesToGaussPointsImpl<DoubleVariableInitializer>(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
    }

    /// Transfer the array_1d variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(ModelPart& rTarget, const Variable<array_1d<double, 3> >& rThisVariable)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << rTarget.Name() << ", Variable<array_1d<double, 3> > " << rThisVariable.Name() << std::endl;

        TransferVariablesToGaussPoints(rTarget.Elements(), rThisVariable, rTarget.GetProcessInfo());
    }

    /// Transfer the double variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(ElementsContainerType& TargetMeshElementsArray, const Variable<array_1d<double, 3> >& rThisVariable, const ProcessInfo& CurrentProcessInfo)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << " Variable<array_1d<double, 3> > " << rThisVariable.Name() << std::endl;

        TransferVariablesToGaussPointsImpl<Array1DVariableInitializer>(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
    }

    /// Transfer the vector variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(ModelPart& rTarget, const Variable<Vector>& rThisVariable, std::size_t ncomponents = 6)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << rTarget.Name() << ", Variable<Vector> " << rThisVariable.Name() << std::endl;

        TransferVariablesToGaussPoints(rTarget.Elements(), rThisVariable, rTarget.GetProcessInfo(), ncomponents);
    }

    /// Transfer the vector variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(ElementsContainerType& TargetMeshElementsArray, const Variable<Vector>& rThisVariable, const ProcessInfo& CurrentProcessInfo, std::size_t ncomponents = 6)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << " Variable<Vector> " << rThisVariable.Name() << std::endl;

        if (ncomponents == 3)
            TransferVariablesToGaussPointsImpl<VectorVariableInitializer<3> >(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
        else if (ncomponents == 6)
            TransferVariablesToGaussPointsImpl<VectorVariableInitializer<6> >(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
    }

protected:

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************

    /// Find the master element candidates that contains the point.
    /// REMARK: we should disable the move mesh flag if we want to search in the reference configuration
    virtual void FindPotentialPartners( const PointType& rSourcePoint, ElementsContainerType& pMasterElements ) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Find an element in pMasterElements contains rSourcePoint and assign it to pTargetElement.
    /// The rLocalTargetPoint is the local point in pTargetElement of rSourcePoint
    /// REMARK: we should disable the move mesh flag if we want to search in the reference configuration
    bool SearchPartner( const PointType& rSourcePoint, ElementsContainerType& pMasterElements,
            Element::Pointer& pTargetElement, PointType& rLocalTargetPoint ) const
    {
        for( ElementsContainerType::ptr_iterator it = pMasterElements.ptr_begin(); it != pMasterElements.ptr_end(); ++it )
        {
            GeometryType& r_geom = (*it)->GetGeometry();

            bool is_inside = r_geom.IsInside( rSourcePoint, rLocalTargetPoint );
            if( is_inside )
            {
                pTargetElement = *it;
                return true;
            }
        }

        if (GetEchoLevel() > 4)
            std::cout << " !!!! WARNING: NO ELEMENT FOUND TO CONTAIN " << rSourcePoint << " !!!! " << std::endl;

        return false;
    }

    /// Interpolate the double value in the element
    void ValueVectorInOldMesh( double& newValue, const Element& oldElement, const PointType& localPoint,
                               const Variable<double>& rThisVariable )
    {
        Vector shape_functions_values;
        shape_functions_values = oldElement.GetGeometry().ShapeFunctionsValues(shape_functions_values, localPoint);

        newValue = 0.0;
        for(unsigned int i = 0; i < oldElement.GetGeometry().size(); ++i)
        {
            const double temp = oldElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);
            newValue += shape_functions_values[i] * temp;
        }
    }

    /// Interpolate the array_1d value in the element
    void ValueVectorInOldMesh( array_1d<double, 3>& newValue, const Element& oldElement, const PointType& localPoint,
                               const Variable<array_1d<double, 3> >& rThisVariable )
    {
        Vector shape_functions_values;
        shape_functions_values = oldElement.GetGeometry().ShapeFunctionsValues(shape_functions_values, localPoint);

        array_1d<double, 3> temp;
        noalias(newValue) = ZeroVector(3);
        for(unsigned int i = 0; i < oldElement.GetGeometry().size(); ++i)
        {
            noalias(temp) = oldElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);
            noalias(newValue) += shape_functions_values[i] * temp;
        }
    }

    /// Interpolate the vector value in the element
    void ValueVectorInOldMesh( Vector& newValue, const Element& oldElement, const PointType& localPoint,
                               const Variable<Vector>& rThisVariable )
    {
        Vector shape_functions_values;
        shape_functions_values = oldElement.GetGeometry().ShapeFunctionsValues(shape_functions_values, localPoint);

        const std::size_t ncomponents = newValue.size();
        noalias(newValue) = ZeroVector(ncomponents);
        Vector temp(ncomponents);
        for(unsigned int i = 0; i < oldElement.GetGeometry().size(); ++i)
        {
            noalias(temp) = oldElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);
            noalias(newValue) += shape_functions_values[i] * temp;
        }
    }

    /// Transfer the variable to Gauss points of the target mesh
    template<class TVariableInitializer>
    void TransferVariablesToGaussPointsImpl( ElementsContainerType& TargetMeshElementsArray,
                                             const ProcessInfo& CurrentProcessInfo,
                                             const typename TVariableInitializer::VariableType& rThisVariable)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints, Variable " << rThisVariable.Name() << std::endl;

        int number_of_threads = 1;
        std::vector<unsigned int> element_partition;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
        double start_transfer = omp_get_wtime();
#endif
        OpenMPUtils::CreatePartition(number_of_threads, TargetMeshElementsArray.size(), element_partition);
        KRATOS_WATCH( number_of_threads );
        std::cout << "element_partition:";
        for (std::size_t i = 0; i < element_partition.size(); ++i)
            std::cout << " " << element_partition[i];
        std::cout << std::endl;
        // KRATOS_WATCH( element_partition );
        Kratos::progress_display show_progress( TargetMeshElementsArray.size() );
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for(int k = 0; k < number_of_threads; ++k)
        {
            ElementsContainerType::ptr_iterator it_begin =
                TargetMeshElementsArray.ptr_begin() + element_partition[k];
            ElementsContainerType::ptr_iterator it_end =
                TargetMeshElementsArray.ptr_begin() + element_partition[k+1];
            for (ElementsContainerType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                if( ((*it)->GetValue(IS_INACTIVE) == true) && !(*it)->Is(ACTIVE) )
                    continue;

                // KRATOS_WATCH((*it)->Id())
                // KRATOS_WATCH(typeid((*it)->GetGeometry()).name())
                const IntegrationPointsArrayType& integration_points
                    = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());
                // KRATOS_WATCH(integration_points.size())
                std::vector<typename TVariableInitializer::DataType> ValuesOnIntPoint(integration_points.size());
                TVariableInitializer::Initialize(ValuesOnIntPoint);
                for(unsigned int point = 0; point< integration_points.size(); ++point)
                {
                    PointType sourceLocalPoint;
                    PointType targetLocalPoint;
                    noalias(targetLocalPoint) = integration_points[point];
                    PointType targetGlobalPoint;
                    (*it)->GetGeometry().GlobalCoordinates(targetGlobalPoint, targetLocalPoint);
//                    KRATOS_WATCH(targetGlobalPoint)
                    ElementsContainerType pMasterElements;
                    this->FindPotentialPartners(targetGlobalPoint, pMasterElements);
                    Element::Pointer sourceElement;
                    //Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
                    bool found = this->SearchPartner( targetGlobalPoint, pMasterElements, sourceElement, sourceLocalPoint );
                    if(found)
                    {
                        // KRATOS_WATCH(sourceElement->Id())
                        // KRATOS_WATCH(typeid(*sourceElement).name())
                        // KRATOS_WATCH(sourceElement->Is(ACTIVE))
                        // KRATOS_WATCH(sourceLocalPoint)
                        ValueVectorInOldMesh( ValuesOnIntPoint[point], *sourceElement, sourceLocalPoint, rThisVariable );
                        // KRATOS_WATCH(ValuesOnIntPoint[point])
                    }
                    else
                    {
                        std::cout << "###### NO PARTNER FOUND IN OLD MESH : TransferVariablesToGaussPoints(..."
                                  << rThisVariable.Name() << "...) at point " << targetGlobalPoint << "#####" << std::endl;
                        continue;
                    }
                }

                (*it)->SetValuesOnIntegrationPoints( rThisVariable, ValuesOnIntPoint, CurrentProcessInfo );

                ++show_progress;
            }
        }

#ifdef _OPENMP
        double stop_transfer = omp_get_wtime();
        std::cout << "TransferVariablesToGaussPoints time: " << stop_transfer - start_transfer << std::endl;
#endif
    }

    /// Transfer the variable at node target mesh
    template<class TVariableInitializer>
    void TransferVariablesToNodesImpl(NodesContainerType& rTargetNodes, const typename TVariableInitializer::VariableType& rThisVariable)
    {
        if (GetEchoLevel() > 0)
            std::cout << __LINE__ << " : At TransferVariablesToNodes, Variable " << rThisVariable.Name() << std::endl;

        int number_of_threads = 1;
        std::vector<unsigned int> node_partition;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
        double start_transfer = omp_get_wtime();
#endif
        OpenMPUtils::CreatePartition(number_of_threads, rTargetNodes.size(), node_partition);
        KRATOS_WATCH( number_of_threads );
        std::cout << "node_partition:";
        for (std::size_t i = 0; i < node_partition.size(); ++i)
            std::cout << " " << node_partition[i];
        std::cout << std::endl;
Kratos::progress_display show_progress( rTargetNodes.size() );
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for(int k = 0; k < number_of_threads; ++k)
        {
            NodesContainerType::ptr_iterator it_begin =
                rTargetNodes.ptr_begin() + node_partition[k];
            NodesContainerType::ptr_iterator it_end =
                rTargetNodes.ptr_begin() + node_partition[k+1];

            typename TVariableInitializer::DataType Tmp;
            for (NodesContainerType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                //Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
                PointType sourceLocalPoint;
                Element::Pointer sourceElement;

                ElementsContainerType pMasterElements;
                this->FindPotentialPartners(*(*it), pMasterElements);
                bool found = this->SearchPartner(*(*it), pMasterElements, sourceElement, sourceLocalPoint);

                TVariableInitializer::Initialize(Tmp);

                if (found)
                {
                    ValueVectorInOldMesh(Tmp, *sourceElement, sourceLocalPoint, rThisVariable);
                }
                else
                {
                    std::cout << "###### NO PARTNER FOUND IN OLD MESH : TransferVariablesToNodes(..."
                              << rThisVariable.Name() << "...) at node " << (*it)->Id() << ", " << (*it)->GetInitialPosition() << "#####"<< std::endl;
                    continue;
                }

                TVariableInitializer::Initialize((*it)->GetSolutionStepValue(rThisVariable), Tmp);

                ++show_progress;
            }
        }
    }
}; // Class VariableInterpolationUtility

} // namespace Kratos.

#endif /* KRATOS_VARIABLE_INTERPOLATION_UTILITY_INCLUDED  defined */
