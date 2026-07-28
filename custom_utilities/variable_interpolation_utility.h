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
#include "includes/variables.h"
#include "utilities/openmp_utils.h"
#include "utilities/progress.h"
#include "utilities/timer.h"
#include "custom_utilities/variable_utility.h"

namespace Kratos
{

/**
 * Utility to transfer the variables with the efficient search functionality
 */
template<class TEntitiesContainerType>
class VariableInterpolationUtility : public VariableUtility<TEntitiesContainerType>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(VariableInterpolationUtility);

    typedef VariableUtility<TEntitiesContainerType> BaseType;
    typedef typename BaseType::EntityType EntityType;
    typedef typename EntityType::GeometryType GeometryType;
    typedef typename GeometryType::PointType NodeType;
    typedef typename NodeType::PointType PointType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef ModelPart::NodesContainerType NodesContainerType;

    using DoubleVariableInitializer = typename BaseType::DoubleVariableInitializer;
    using Array1DVariableInitializer = typename BaseType::Array1DVariableInitializer;

    template<int TSize>
    using VectorVariableInitializer = typename BaseType::template VectorVariableInitializer<TSize>;

    VariableInterpolationUtility(const TEntitiesContainerType& pElements)
        : BaseType(pElements), mSearchTolerance(1e-8)
    {
    }

    VariableInterpolationUtility(const TEntitiesContainerType& pElements, const int EchoLevel)
        : BaseType(pElements, EchoLevel), mSearchTolerance(1e-8)
    {
    }

    ~VariableInterpolationUtility() override
    {
    }

    /// Set the search tolerance
    void SetSearchTolerance(double value)
    {
        mSearchTolerance = value;
    }

    /// Get the elements of which the BV contains the point
    TEntitiesContainerType FindPotentialPartners( const PointType& rSourcePoint ) const
    {
        TEntitiesContainerType pMasterElements;
        this->FindPotentialPartners( rSourcePoint, pMasterElements );
        return pMasterElements;
    }

    /// Get the element containing the point
    typename EntityType::Pointer SearchPartner( const PointType& rSourcePoint, TEntitiesContainerType& pMasterElements ) const
    {
        typename EntityType::Pointer pElement;

        PointType localPoint;
        this->SearchPartner( rSourcePoint, pMasterElements, pElement, localPoint);

        return pElement;
    }

    /// Transfer the double variable to node of the target model_part
    void TransferVariablesToNodes(ModelPart& rTarget, const Variable<double>& rThisVariable) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << rTarget.Name() << ", Variable<double> " << rThisVariable.Name() << std::endl;
        }

        TransferVariablesToNodes(rTarget.Nodes(), rThisVariable);
    }

    /// Transfer the double variable to node of the target node mesh
    void TransferVariablesToNodes(NodesContainerType& rTargetNodes, const Variable<double>& rThisVariable) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << " Variable<double> " << rThisVariable.Name() << std::endl;
        }

        TransferVariablesToNodesImpl<DoubleVariableInitializer>(rTargetNodes, rThisVariable);
    }

    /// Transfer the double variable to node of the target model_part
    void TransferVariablesToNodes(ModelPart& rTarget, const Variable<array_1d<double, 3> >& rThisVariable) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << rTarget.Name() << ", Variable<array_1d<double, 3> > " << rThisVariable.Name() << std::endl;
        }

        TransferVariablesToNodes(rTarget.Nodes(), rThisVariable);
    }

    /// Transfer the double variable to node of the target node mesh
    void TransferVariablesToNodes(NodesContainerType& rTargetNodes, const Variable<array_1d<double, 3> >& rThisVariable) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << " Variable<array_1d<double, 3> > " << rThisVariable.Name() << std::endl;
        }

        TransferVariablesToNodesImpl<Array1DVariableInitializer>(rTargetNodes, rThisVariable);
    }

    /// Transfer the double variable to node of the target model_part
    void TransferVariablesToNodes(ModelPart& rTarget, const Variable<Vector>& rThisVariable, const std::size_t& ncomponents) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << rTarget.Name() << ", Variable<Vector> " << rThisVariable.Name() << std::endl;
        }

        TransferVariablesToNodes(rTarget.Nodes(), rThisVariable, ncomponents);
    }

    /// Transfer the double variable to node of the target node mesh
    void TransferVariablesToNodes(NodesContainerType& rTargetNodes, const Variable<Vector>& rThisVariable, const std::size_t& ncomponents) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToNodes(" << " Variable<Vector> " << rThisVariable.Name() << std::endl;
        }

        if (ncomponents == 3)
        {
            TransferVariablesToNodesImpl<VectorVariableInitializer<3> >(rTargetNodes, rThisVariable);
        }
        else if (ncomponents == 6)
        {
            TransferVariablesToNodesImpl<VectorVariableInitializer<6> >(rTargetNodes, rThisVariable);
        }
    }

    /// Transfer the double variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(ModelPart& rTarget, const Variable<double>& rThisVariable) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << rTarget.Name() << ", Variable<double> " << rThisVariable.Name() << std::endl;
        }

        if constexpr (std::is_same<TEntitiesContainerType, ModelPart::ElementsContainerType>::value)
        {
            TransferVariablesToGaussPoints(rTarget.Elements(), rThisVariable, rTarget.GetProcessInfo());
        }
        else if constexpr (std::is_same<TEntitiesContainerType, ModelPart::ConditionsContainerType>::value)
        {
            TransferVariablesToGaussPoints(rTarget.Conditions(), rThisVariable, rTarget.GetProcessInfo());
        }
    }

    /// Transfer the double variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(TEntitiesContainerType& TargetMeshElementsArray, const Variable<double>& rThisVariable, const ProcessInfo& CurrentProcessInfo) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << " Variable<double> " << rThisVariable.Name() << std::endl;
        }

        TransferVariablesToGaussPointsImpl<DoubleVariableInitializer>(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
    }

    /// Transfer the array_1d variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(ModelPart& rTarget, const Variable<array_1d<double, 3> >& rThisVariable) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << rTarget.Name() << ", Variable<array_1d<double, 3> > " << rThisVariable.Name() << std::endl;
        }

        if constexpr (std::is_same<TEntitiesContainerType, ModelPart::ElementsContainerType>::value)
        {
            TransferVariablesToGaussPoints(rTarget.Elements(), rThisVariable, rTarget.GetProcessInfo());
        }
        else if constexpr (std::is_same<TEntitiesContainerType, ModelPart::ConditionsContainerType>::value)
        {
            TransferVariablesToGaussPoints(rTarget.Conditions(), rThisVariable, rTarget.GetProcessInfo());
        }
    }

    /// Transfer the double variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(TEntitiesContainerType& TargetMeshElementsArray, const Variable<array_1d<double, 3> >& rThisVariable, const ProcessInfo& CurrentProcessInfo) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << " Variable<array_1d<double, 3> > " << rThisVariable.Name() << std::endl;
        }

        TransferVariablesToGaussPointsImpl<Array1DVariableInitializer>(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
    }

    /// Transfer the vector variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(ModelPart& rTarget, const Variable<Vector>& rThisVariable, std::size_t ncomponents = 6) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << rTarget.Name() << ", Variable<Vector> " << rThisVariable.Name() << std::endl;
        }

        if constexpr (std::is_same<TEntitiesContainerType, ModelPart::ElementsContainerType>::value)
        {
            TransferVariablesToGaussPoints(rTarget.Elements(), rThisVariable, rTarget.GetProcessInfo(), ncomponents);
        }
        else if constexpr (std::is_same<TEntitiesContainerType, ModelPart::ConditionsContainerType>::value)
        {
            TransferVariablesToGaussPoints(rTarget.Conditions(), rThisVariable, rTarget.GetProcessInfo(), ncomponents);
        }
    }

    /// Transfer the vector variable to Gauss points of the target model_part
    void TransferVariablesToGaussPoints(TEntitiesContainerType& TargetMeshElementsArray, const Variable<Vector>& rThisVariable,
            const ProcessInfo& CurrentProcessInfo, std::size_t ncomponents = 6) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints(" << " Variable<Vector> " << rThisVariable.Name() << std::endl;
        }

        if (ncomponents == 3)
        {
            TransferVariablesToGaussPointsImpl<VectorVariableInitializer<3> >(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
        }
        else if (ncomponents == 6)
        {
            TransferVariablesToGaussPointsImpl<VectorVariableInitializer<6> >(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
        }
        else
            KRATOS_ERROR << "Number of component = " << ncomponents << " is not supported";
    }

    /// Compute the L-2 norm of the difference between two meshes for a double variable.
    /// It can be useful to compute the error norm when the master mesh is very fine and contains very accurate solution.
    double ComputeRelativeDifference(ModelPart& rTarget, const Variable<double>& rThisVariable,
            const ProcessInfo& CurrentProcessInfo) const
    {
        return ComputeRelativeDifference(GetEntities(rTarget), rThisVariable, CurrentProcessInfo);
    }

    /// Compute the L-2 norm of the difference between two meshes for a double variable.
    /// It can be useful to compute the error norm when the master mesh is very fine and contains very accurate solution.
    double ComputeRelativeDifference(TEntitiesContainerType& TargetMeshElementsArray,
            const Variable<double>& rThisVariable, const ProcessInfo& CurrentProcessInfo) const
    {
        return ComputeRelativeDifferenceImpl<DoubleVariableInitializer>(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
    }

    /// Compute the L-2 norm of the difference between two meshes for an array_1d variable.
    /// It can be useful to compute the error norm when the master mesh is very fine and contains very accurate solution.
    double ComputeRelativeDifference(ModelPart& rTarget, const Variable<array_1d<double, 3> >& rThisVariable,
            const ProcessInfo& CurrentProcessInfo) const
    {
        return ComputeRelativeDifference(GetEntities(rTarget), rThisVariable, CurrentProcessInfo);
    }

    /// Compute the L-2 norm of the difference between two meshes for an array_1d variable.
    /// It can be useful to compute the error norm when the master mesh is very fine and contains very accurate solution.
    double ComputeRelativeDifference(TEntitiesContainerType& TargetMeshElementsArray, const Variable<array_1d<double, 3> >& rThisVariable,
            const ProcessInfo& CurrentProcessInfo) const
    {
        return ComputeRelativeDifferenceImpl<Array1DVariableInitializer>(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
    }

    /// Compute the L-2 norm of the difference between two meshes for an array_1d variable.
    /// It can be useful to compute the error norm when the master mesh is very fine and contains very accurate solution.
    double ComputeRelativeDifference(ModelPart& rTarget, const Variable<Vector>& rThisVariable,
            const ProcessInfo& CurrentProcessInfo, std::size_t ncomponents = 6) const
    {
        return ComputeRelativeDifference(GetEntities(rTarget), rThisVariable, CurrentProcessInfo, ncomponents);
    }

    /// Compute the L-2 norm of the difference between two meshes for an array_1d variable.
    /// It can be useful to compute the error norm when the master mesh is very fine and contains very accurate solution.
    double ComputeRelativeDifference(TEntitiesContainerType& TargetMeshElementsArray,
            const Variable<Vector>& rThisVariable,
            const ProcessInfo& CurrentProcessInfo,
            std::size_t ncomponents = 6) const
    {
        if (ncomponents == 3)
        {
            return ComputeRelativeDifferenceImpl<VectorVariableInitializer<3> >(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
        }
        else if (ncomponents == 6)
        {
            return ComputeRelativeDifferenceImpl<VectorVariableInitializer<6> >(TargetMeshElementsArray, CurrentProcessInfo, rThisVariable);
        }
        else
            KRATOS_ERROR << "Number of component = " << ncomponents << " is not supported";

        return 0.0;
    }

    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return "VariableInterpolationUtility";
    }

    ///@}

protected:

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************

    /// Get the corresponding elements or conditions from the model_part
    static TEntitiesContainerType& GetEntities(ModelPart& rTarget)
    {
        if constexpr (std::is_same<TEntitiesContainerType, ModelPart::ElementsContainerType>::value)
        {
            return rTarget.Elements();
        }
        else if constexpr (std::is_same<TEntitiesContainerType, ModelPart::ConditionsContainerType>::value)
        {
            return rTarget.Conditions();
        }
        else
            KRATOS_ERROR << "Invalid operation";
    }

    /// Find the master element candidates that contains the point.
    /// REMARK: we should disable the move mesh flag if we want to search in the reference configuration
    virtual void FindPotentialPartners( const PointType& rSourcePoint, TEntitiesContainerType& pMasterElements ) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Find an element in pMasterElements contains rSourcePoint and assign it to pTargetElement.
    /// The rLocalTargetPoint is the local point in pTargetElement of rSourcePoint
    /// REMARK: we should disable the move mesh flag if we want to search in the reference configuration
    bool SearchPartner( const PointType& rSourcePoint, TEntitiesContainerType& pMasterElements,
                        typename EntityType::Pointer& pTargetElement, PointType& rLocalTargetPoint ) const
    {
        for ( auto it = pMasterElements.ptr_begin(); it != pMasterElements.ptr_end(); ++it )
        {
            const GeometryType& r_geom = (*it)->GetGeometry();

            r_geom.PointLocalCoordinates( rLocalTargetPoint, rSourcePoint, true, mSearchTolerance );
            bool is_inside = r_geom.IsInside( rLocalTargetPoint, mSearchTolerance );
            if ( is_inside )
            {
                pTargetElement = *it;
                return true;
            }
        }

        if (this->GetEchoLevel() > 4)
        {
            std::cout << " !!!! WARNING: NO ELEMENT FOUND TO CONTAIN " << rSourcePoint << " !!!! " << std::endl;
        }

        return false;
    }

    /// Interpolate the double value in the element
    void ValueVectorInOldMesh( double& newValue, const EntityType& oldElement, const PointType& localPoint,
                               const Variable<double>& rThisVariable ) const
    {
        Vector shape_functions_values;
        shape_functions_values = oldElement.GetGeometry().ShapeFunctionsValues(shape_functions_values, localPoint);

        newValue = 0.0;
        for (unsigned int i = 0; i < oldElement.GetGeometry().size(); ++i)
        {
            const double temp = oldElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);
            newValue += shape_functions_values[i] * temp;
        }
    }

    /// Interpolate the array_1d value in the element
    void ValueVectorInOldMesh( array_1d<double, 3>& newValue, const EntityType& oldElement, const PointType& localPoint,
                               const Variable<array_1d<double, 3> >& rThisVariable ) const
    {
        Vector shape_functions_values;
        shape_functions_values = oldElement.GetGeometry().ShapeFunctionsValues(shape_functions_values, localPoint);

        array_1d<double, 3> temp;
        noalias(newValue) = ZeroVector(3);
        for (unsigned int i = 0; i < oldElement.GetGeometry().size(); ++i)
        {
            noalias(temp) = oldElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);
            noalias(newValue) += shape_functions_values[i] * temp;
        }
    }

    /// Interpolate the vector value in the element
    void ValueVectorInOldMesh( Vector& newValue, const EntityType& oldElement, const PointType& localPoint,
                               const Variable<Vector>& rThisVariable ) const
    {
        Vector shape_functions_values;
        shape_functions_values = oldElement.GetGeometry().ShapeFunctionsValues(shape_functions_values, localPoint);

        const std::size_t ncomponents = newValue.size();
        noalias(newValue) = ZeroVector(ncomponents);
        Vector temp(ncomponents);
        for (unsigned int i = 0; i < oldElement.GetGeometry().size(); ++i)
        {
            noalias(temp) = oldElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);
            noalias(newValue) += shape_functions_values[i] * temp;
        }
    }

    /// Transfer the variable to Gauss points of the target mesh
    template<class TVariableInitializer>
    void TransferVariablesToGaussPointsImpl( TEntitiesContainerType& TargetMeshElementsArray,
            const ProcessInfo& CurrentProcessInfo,
            const typename TVariableInitializer::VariableType& rThisVariable) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToGaussPoints, Variable " << rThisVariable.Name() << std::endl;
        }

        int number_of_threads = 1;
        std::vector<unsigned int> element_partition;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
        double start_transfer = omp_get_wtime();
#endif
        OpenMPUtils::CreatePartition(number_of_threads, TargetMeshElementsArray.size(), element_partition);
        KRATOS_WATCH( number_of_threads );
        KRATOS_WATCH_STD_CON( element_partition )
        Kratos::progress_display show_progress( TargetMeshElementsArray.size() );
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int k = 0; k < number_of_threads; ++k)
        {
            auto it_begin = TargetMeshElementsArray.ptr_begin() + element_partition[k];
            auto it_end = TargetMeshElementsArray.ptr_begin() + element_partition[k + 1];
            for (auto it = it_begin; it != it_end; ++it)
            {
                if ( ((*it)->GetValue(IS_INACTIVE) == true) && !(*it)->Is(ACTIVE) )
                {
                    continue;
                }

                // KRATOS_WATCH((*it)->Id())
                // KRATOS_WATCH(typeid((*it)->GetGeometry()).name())
                const IntegrationPointsArrayType& integration_points
                    = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());
                // KRATOS_WATCH(integration_points.size())
                std::vector<typename TVariableInitializer::DataType> ValuesOnIntPoint(integration_points.size());
                TVariableInitializer::Initialize(ValuesOnIntPoint);
                for (unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    PointType sourceLocalPoint;
                    PointType targetLocalPoint;
                    noalias(targetLocalPoint) = integration_points[point];
                    PointType targetGlobalPoint;
                    (*it)->GetGeometry().GlobalCoordinates(targetGlobalPoint, targetLocalPoint);
//                    KRATOS_WATCH(targetGlobalPoint)
                    TEntitiesContainerType pMasterElements;
                    this->FindPotentialPartners(targetGlobalPoint, pMasterElements);
                    typename EntityType::Pointer sourceElement;
                    //Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
                    bool found = this->SearchPartner( targetGlobalPoint, pMasterElements, sourceElement, sourceLocalPoint );
                    if (found)
                    {
                        ValueVectorInOldMesh( ValuesOnIntPoint[point], *sourceElement, sourceLocalPoint, rThisVariable );
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
    void TransferVariablesToNodesImpl(NodesContainerType& rTargetNodes, const typename TVariableInitializer::VariableType& rThisVariable) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At TransferVariablesToNodes, Variable " << rThisVariable.Name() << std::endl;
        }

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
        {
            std::cout << " " << node_partition[i];
        }
        std::cout << std::endl;
        Kratos::progress_display show_progress( rTargetNodes.size() );
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int k = 0; k < number_of_threads; ++k)
        {
            NodesContainerType::ptr_iterator it_begin =
                rTargetNodes.ptr_begin() + node_partition[k];
            NodesContainerType::ptr_iterator it_end =
                rTargetNodes.ptr_begin() + node_partition[k + 1];

            typename TVariableInitializer::DataType Tmp;
            for (NodesContainerType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                //Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
                PointType sourceLocalPoint;
                typename EntityType::Pointer sourceElement;

                TEntitiesContainerType pMasterElements;
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
                              << rThisVariable.Name() << "...) at node " << (*it)->Id() << ", " << (*it)->GetInitialPosition() << "#####" << std::endl;
                    continue;
                }

                TVariableInitializer::Initialize((*it)->GetSolutionStepValue(rThisVariable), Tmp);

                ++show_progress;
            }
        }
    }

    /// Compute the relative difference: sqrt(int_{\Omega^s} ||u_s - u_m||^2 dV) / sqrt(int__{\Omega^s} ||u_m||^2 dV)
    template<class TVariableInitializer>
    double ComputeRelativeDifferenceImpl( TEntitiesContainerType& TargetMeshElementsArray,
            const ProcessInfo& CurrentProcessInfo,
            const typename TVariableInitializer::VariableType& rThisVariable) const
    {
        if (this->GetEchoLevel() > 0)
        {
            std::cout << __LINE__ << " : At ComputeRelativeDifference, Variable " << rThisVariable.Name() << std::endl;
        }

        int number_of_threads = 1;
        std::vector<unsigned int> element_partition;
        std::vector<double> error_partition;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
        double start_compute = omp_get_wtime();
#endif
        OpenMPUtils::CreatePartition(number_of_threads, TargetMeshElementsArray.size(), element_partition);
        error_partition.resize(number_of_threads);
        std::fill(error_partition.begin(), error_partition.end(), 0.0);
        KRATOS_WATCH( number_of_threads );
        KRATOS_WATCH_STD_CON( element_partition )
        Kratos::progress_display show_progress( TargetMeshElementsArray.size() );
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int k = 0; k < number_of_threads; ++k)
        {
            auto it_begin = TargetMeshElementsArray.begin() + element_partition[k];
            auto it_end = TargetMeshElementsArray.begin() + element_partition[k + 1];
            double nom = 0.0, denom = 0.0;
            for (auto it = it_begin; it != it_end; ++it)
            {
                if ( (it->GetValue(IS_INACTIVE) == true) && !it->Is(ACTIVE) )
                {
                    continue;
                }

                const IntegrationPointsArrayType& integration_points
                    = it->GetGeometry().IntegrationPoints(it->GetIntegrationMethod());

                typename GeometryType::JacobiansType J0;
                it->GetGeometry().Jacobian0(J0, it->GetIntegrationMethod());

                // Extract value of rVariable in target mesh
                std::vector<typename TVariableInitializer::DataType> ValuesOnIntPoint(integration_points.size());
                it->CalculateOnIntegrationPoints( rThisVariable, ValuesOnIntPoint, CurrentProcessInfo );

                for (unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    // Calculate value of rVariable in source mesh
                    PointType targetLocalPoint;
                    noalias(targetLocalPoint) = integration_points[point];
                    PointType targetGlobalPoint;
                    it->GetGeometry().GlobalCoordinates(targetGlobalPoint, targetLocalPoint);

                    TEntitiesContainerType pMasterElements;
                    this->FindPotentialPartners(targetGlobalPoint, pMasterElements);

                    PointType sourceLocalPoint;
                    typename EntityType::Pointer sourceElement;
                    typename TVariableInitializer::DataType sourceValue;
                    TVariableInitializer::Initialize( sourceValue );
                    bool found = this->SearchPartner( targetGlobalPoint, pMasterElements, sourceElement, sourceLocalPoint );
                    if (found)
                    {
                        ValueVectorInOldMesh( sourceValue, *sourceElement, sourceLocalPoint, rThisVariable );
                    }
                    else
                    {
                        std::cout << "###### NO PARTNER FOUND IN OLD MESH : TransferVariablesToGaussPoints(..."
                                  << rThisVariable.Name() << "...) at point " << targetGlobalPoint << "#####" << std::endl;
                        continue;
                    }

                    // Calculate the difference
                    double IntegrationWeight = integration_points[point].Weight();
                    IntegrationWeight *= std::sqrt(MathUtils<double>::Det(Matrix(prod(trans(J0[point]), J0[point]))));

                    const auto& targetValue = ValuesOnIntPoint[point];
// std::cout << "  targetValue: " << targetValue << ", sourceValue: " << sourceValue << std::endl;
                    double tmp1 = TVariableInitializer::Norm(targetValue - sourceValue);
                    double tmp2 = TVariableInitializer::Norm(sourceValue);
                    nom += tmp1*tmp1*IntegrationWeight;
                    denom += tmp2*tmp2*IntegrationWeight;
                }

                ++show_progress;
            }

            if (denom == 0.0)
                KRATOS_ERROR << "The value of Variable " << rThisVariable << " in source mesh is zero";

            error_partition[k] = std::sqrt(nom / denom);
        }

#ifdef _OPENMP
        double stop_compute = omp_get_wtime();
        std::cout << "ComputeRelativeDifference time: " << stop_compute - start_compute << std::endl;
#endif

        return std::accumulate(error_partition.begin(), error_partition.end(), 0.0);
    }

private:

    double mSearchTolerance;

}; // Class VariableInterpolationUtility

} // namespace Kratos.

#endif /* KRATOS_VARIABLE_INTERPOLATION_UTILITY_INCLUDED  defined */
