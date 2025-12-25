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
*   Date:                $Date: 24 May 2018 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_VARIABLE_PROJECTION_UTILITY_INCLUDED )
#define  KRATOS_VARIABLE_PROJECTION_UTILITY_INCLUDED

//System includes
#ifdef _OPENMP
#include <omp.h>
#endif

//External includes

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "utilities/openmp_utils.h"
#include "utilities/progress.h"
#include "custom_utilities/variable_utility.h"

namespace Kratos
{

/**
 * Utility to project the variables from Gauss points to nodes using L2-projection.
 */
// The transfer of Gaussian variables to nodal Variables is via L_2-Minimization
// see Jiao + Heath "Common-refinement-based data tranfer ..."
// International Journal for numerical methods in engineering 61 (2004) 2402--2427
// for general description of L_2-Minimization
template<class TEntitiesContainerType>
class VariableProjectionUtility : public VariableUtility<TEntitiesContainerType>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( VariableProjectionUtility );

    typedef VariableUtility<TEntitiesContainerType> BaseType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::EntityType EntityType;
    typedef typename BaseType::SparseSpaceType SparseSpaceType;
    typedef typename BaseType::DenseSpaceType DenseSpaceType;
    typedef typename BaseType::LinearSolverType LinearSolverType;
    typedef typename EntityType::GeometryType GeometryType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename GeometryType::JacobiansType JacobiansType;
    typedef ModelPart::NodesContainerType NodesContainerType;

    /**
     * Constructor.
     */
    VariableProjectionUtility(const TEntitiesContainerType& pElements, LinearSolverType::Pointer pLinearSolver)
    : BaseType(pElements), mpLinearSolver(pLinearSolver)
    {
        this->Initialize(pElements);
        std::cout << "VariableProjectionUtility created" << std::endl;
    }

    VariableProjectionUtility(const TEntitiesContainerType& pElements, LinearSolverType::Pointer pLinearSolver, const int EchoLevel)
    : BaseType(pElements, EchoLevel), mpLinearSolver(pLinearSolver)
    {
        this->Initialize(pElements);
        if (this->GetEchoLevel() > 0)
            std::cout << "VariableProjectionUtility created" << std::endl;
    }

    /**
     * Destructor.
     */
    ~VariableProjectionUtility() override
    {}

    /**
     * Operations
     */

    /**
     * Compute the projected nodal values, based on the integration point values
     */
    std::map<IndexType, double> ComputeNodalValues(const std::map<IndexType, std::vector<double> >& rElementalValues)
    {
        NodesContainerType pActiveNodes;
        std::map<std::size_t, std::size_t> NodeRowId;
        this->ExtractActiveNodes(BaseType::mpElements, pActiveNodes, NodeRowId);
        // KRATOS_WATCH(pActiveNodes.size())

        if (pActiveNodes.size() == 0)
        {
            std::cout << "Number of active nodes is null. There is nothing to transfer" << std::endl;
        }

        typename SparseSpaceType::VectorType g(pActiveNodes.size());
        typename SparseSpaceType::VectorType b(pActiveNodes.size());

        int number_of_threads = 1;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
#endif

#ifdef _OPENMP
        //create the array of lock
        unsigned int M_size = pActiveNodes.size();
        std::vector<omp_lock_t> lock_array(M_size);
        for (unsigned int i = 0; i < M_size; ++i)
            omp_init_lock(&lock_array[i]);
#endif

        noalias(g) = ZeroVector(M_size);
        noalias(b) = ZeroVector(M_size);

        std::vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, BaseType::mpElements.size(), element_partition);

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int k = 0; k < number_of_threads; ++k)
        {
            auto it_begin = BaseType::mpElements.begin() + element_partition[k];
            auto it_end = BaseType::mpElements.begin() + element_partition[k+1];
            std::map<std::size_t, std::size_t>::const_iterator it_id;

            for( auto it = it_begin; it != it_end; ++it )
            {
                if ( ! ( ( it->GetValue(IS_INACTIVE) == false ) || it->Is(ACTIVE) ) )
                    continue;

                auto& rGeometry = it->GetGeometry();

                #ifdef ENABLE_BEZIER_GEOMETRY
                rGeometry.Initialize(it->GetIntegrationMethod());
                #endif

                const IntegrationPointsArrayType& integration_points
                    = rGeometry.IntegrationPoints( it->GetIntegrationMethod());

                JacobiansType J(integration_points.size());
                J = rGeometry.Jacobian(J, it->GetIntegrationMethod());

                auto itv = rElementalValues.find(it->Id());
                if (itv == rElementalValues.end())
                    KRATOS_ERROR << "Elemental values for element " << it->Id() << " is not provided";

                const std::vector<double>& ValuesOnIntPoint = itv->second;
                // KRATOS_WATCH_STD_CON(ValuesOnIntPoint)

                if (ValuesOnIntPoint.size() != integration_points.size())
                    KRATOS_ERROR << "The size of provided input "
                                 << "(" << itv->second.size() << ")"
                                 << " for element " << it->Id() << " is not sufficient";

                const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues(it->GetIntegrationMethod());

                double DetJ;
                for (unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    if (rGeometry.WorkingSpaceDimension() == rGeometry.LocalSpaceDimension())
                        DetJ = MathUtils<double>::Det(J[point]);
                    else
                        DetJ = std::sqrt(MathUtils<double>::Det(Matrix(prod(trans(J[point]), J[point]))));

                    double dV = DetJ*integration_points[point].Weight();

                    for (unsigned int prim = 0; prim < rGeometry.size(); ++prim)
                    {
                        unsigned int row = NodeRowId[rGeometry[prim].Id()];
#ifdef _OPENMP
                        omp_set_lock(&lock_array[row]);
#endif
                        b(row) += ValuesOnIntPoint[point] * Ncontainer(point, prim) * dV;
#ifdef _OPENMP
                        omp_unset_lock(&lock_array[row]);
#endif
                    }
                }

                #ifdef ENABLE_BEZIER_GEOMETRY
                rGeometry.Clean();
                #endif
            }
        }

#ifdef _OPENMP
        for(unsigned int i = 0; i < M_size; ++i)
            omp_destroy_lock(&lock_array[i]);
#endif

        mpLinearSolver->Solve(mProjectionMatrix, g, b);

        std::map<IndexType, double> results;
        for(auto it = pActiveNodes.begin() ; it != pActiveNodes.end() ; it++)
        {
            results[it->Id()] = g(NodeRowId[it->Id()]);
        }

        return results;
    }

    /**
     * Transfer the double variable from Gauss point to nodes
     * Here a different element set is allowed to compute the right hand side. This element set shall have the same connectivities as the element set used in Initialize.
     */
    void TransferVariablesToNodes( const Variable<double>& rThisVariable, const ProcessInfo& rProcessInfo )
    {
        this->TransferVariablesToNodes( rThisVariable, rThisVariable, rProcessInfo );
    }

    /**
     * Transfer the double variable from Gauss point to nodes
     * Here a different element set is allowed to compute the right hand side. This element set shall have the same connectivities as the element set used in Initialize.
     */
    void TransferVariablesToNodes( const Variable<double>& rIntegrationPointVariable, const Variable<double>& rNodalVariable, const ProcessInfo& rProcessInfo )
    {
        NodesContainerType pActiveNodes;
        std::map<std::size_t, std::size_t> NodeRowId;
        this->ExtractActiveNodes(BaseType::mpElements, pActiveNodes, NodeRowId);
        // KRATOS_WATCH(pActiveNodes.size())

        if (pActiveNodes.size() == 0)
        {
            std::cout << "Number of active nodes is null. There is nothing to transfer" << std::endl;
        }

        typename SparseSpaceType::VectorType g(pActiveNodes.size());
        typename SparseSpaceType::VectorType b(pActiveNodes.size());

        int number_of_threads = 1;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
#endif

#ifdef _OPENMP
        //create the array of lock
        unsigned int M_size = pActiveNodes.size();
        std::vector<omp_lock_t> lock_array(M_size);
        for (unsigned int i = 0; i < M_size; ++i)
            omp_init_lock(&lock_array[i]);
#endif

        noalias(g) = ZeroVector(M_size);
        noalias(b) = ZeroVector(M_size);

        std::vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, BaseType::mpElements.size(), element_partition);

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int k = 0; k < number_of_threads; ++k)
        {
            auto it_begin = BaseType::mpElements.begin() + element_partition[k];
            auto it_end = BaseType::mpElements.begin() + element_partition[k+1];
            std::map<std::size_t, std::size_t>::const_iterator it_id;

            for( auto it = it_begin; it != it_end; ++it )
            {
                if ( ! ( ( it->GetValue(IS_INACTIVE) == false ) || it->Is(ACTIVE) ) )
                    continue;

                auto& rGeometry = it->GetGeometry();

                #ifdef ENABLE_BEZIER_GEOMETRY
                rGeometry.Initialize(it->GetIntegrationMethod());
                #endif

                const IntegrationPointsArrayType& integration_points
                    = rGeometry.IntegrationPoints( it->GetIntegrationMethod());

                JacobiansType J(integration_points.size());
                J = rGeometry.Jacobian(J, it->GetIntegrationMethod());

                std::vector<double> ValuesOnIntPoint(integration_points.size());
                it->CalculateOnIntegrationPoints(rIntegrationPointVariable, ValuesOnIntPoint, rProcessInfo);

                const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues(it->GetIntegrationMethod());

                double DetJ;
                for (unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    if (rGeometry.WorkingSpaceDimension() == rGeometry.LocalSpaceDimension())
                        DetJ = MathUtils<double>::Det(J[point]);
                    else
                        DetJ = sqrt(MathUtils<double>::Det(Matrix(prod(trans(J[point]), J[point]))));

                    double dV = DetJ*integration_points[point].Weight();

                    for (unsigned int prim = 0; prim < rGeometry.size(); ++prim)
                    {
                        unsigned int row = NodeRowId[rGeometry[prim].Id()];
#ifdef _OPENMP
                        omp_set_lock(&lock_array[row]);
#endif
                        b(row) += ValuesOnIntPoint[point] * Ncontainer(point, prim) * dV;
#ifdef _OPENMP
                        omp_unset_lock(&lock_array[row]);
#endif
                    }
                }

                #ifdef ENABLE_BEZIER_GEOMETRY
                rGeometry.Clean();
                #endif
            }
        }

#ifdef _OPENMP
        for(unsigned int i = 0; i < M_size; ++i)
            omp_destroy_lock(&lock_array[i]);
#endif

        mpLinearSolver->Solve(mProjectionMatrix, g, b);

        for(auto it = pActiveNodes.begin() ; it != pActiveNodes.end() ; it++)
        {
            it->GetSolutionStepValue(rNodalVariable) = g(NodeRowId[it->Id()]);
        }

        std::cout << "TransferVariablesToNodes from " << rIntegrationPointVariable.Name()
                  << " to " << rNodalVariable.Name()
                  << " completed" << std::endl;
    }

protected:

    /// Initialize the projection matrix
    void Initialize( const TEntitiesContainerType& pElements ) final
    {
        NodesContainerType pActiveNodes;
        std::map<std::size_t, std::size_t> NodeRowId;
        this->ExtractActiveNodes(pElements, pActiveNodes, NodeRowId);
        this->ConstructLHSMatrix(mProjectionMatrix, pElements, NodeRowId);
    }

private:

    LinearSolverType::Pointer mpLinearSolver;
    SparseSpaceType::MatrixType mProjectionMatrix;

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    /// Extract the active nodes and construct the nodal id row map
    void ExtractActiveNodes(const TEntitiesContainerType& pElements,
        NodesContainerType& pActiveNodes, std::map<std::size_t, std::size_t>& NodeRowId) const
    {
        // extract the active nodes
        for(auto it = pElements.begin(); it != pElements.end(); ++it)
        {
            if( it->GetValue(IS_INACTIVE) == false || it->Is(ACTIVE) )
            {
                for(std::size_t i = 0; i < it->GetGeometry().size(); ++i)
                {
                    pActiveNodes.push_back( it->GetGeometry().pGetPoint(i) );
                }
            }
        }

        pActiveNodes.Unique();
        // KRATOS_WATCH(pActiveNodes.size())

        // assign each node an id. That id is the row of this node in the global L2 projection matrix*/
        std::size_t cnt = 0;
        for(auto it = pActiveNodes.begin(); it != pActiveNodes.end(); ++it )
        {
            NodeRowId[it->Id()] = cnt++;
        }
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    /// Construct the L2 projection matrix
    void ConstructLHSMatrix(SparseSpaceType::MatrixType& rA, const TEntitiesContainerType& pElements,
            const std::map<std::size_t, std::size_t>& NodeRowId) const
    {
        // set up system
        if(rA.size1() != NodeRowId.size() || rA.size2() != NodeRowId.size())
            rA.resize(NodeRowId.size(), NodeRowId.size(), false);
        noalias(rA) = ZeroMatrix(NodeRowId.size(), NodeRowId.size());

        int number_of_threads = 1;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
#endif
        std::vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);
//        Kratos::progress_display show_progress( pElements.size() );

        // create the structure for M a priori
        this->ConstructMatrixStructure(rA, pElements, NodeRowId);

#ifdef _OPENMP
        //create the array of lock
        std::vector< omp_lock_t > lock_array(rA.size1());
        unsigned int system_size = rA.size1();
        for(unsigned int i = 0; i < system_size; ++i)
            omp_init_lock(&lock_array[i]);
#endif

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for(int k = 0; k < number_of_threads; ++k)
        {
            auto it_begin = pElements.begin() + element_partition[k];
            auto it_end = pElements.begin() + element_partition[k+1];
            std::map<std::size_t, std::size_t>::const_iterator it_id;

            for( auto it = it_begin; it != it_end; ++it )
            {
                if( it->GetValue(IS_INACTIVE) == true && !it->Is(ACTIVE) )
                    continue;

                auto& rGeometry = it->GetGeometry();

                #ifdef ENABLE_BEZIER_GEOMETRY
                rGeometry.Initialize(it->GetIntegrationMethod());
                #endif

                unsigned int dim = rGeometry.WorkingSpaceDimension();

                const IntegrationPointsArrayType& integration_points
                    = rGeometry.IntegrationPoints(it->GetIntegrationMethod());

                JacobiansType J(integration_points.size());
                J = rGeometry.Jacobian(J, it->GetIntegrationMethod());

                const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues(it->GetIntegrationMethod());

                double DetJ;
                for(unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    if (rGeometry.WorkingSpaceDimension() == rGeometry.LocalSpaceDimension())
                        DetJ = MathUtils<double>::Det(J[point]);
                    else
                        DetJ = std::sqrt(MathUtils<double>::Det(Matrix(prod(trans(J[point]), J[point]))));

                    double dV = DetJ*integration_points[point].Weight();

                    for(unsigned int prim = 0 ; prim < rGeometry.size(); ++prim)
                    {
                        it_id = NodeRowId.find(rGeometry[prim].Id());
                        std::size_t row = it_id->second;
#ifdef _OPENMP
                        omp_set_lock(&lock_array[row]);
#endif
                        for(unsigned int sec = 0 ; sec < rGeometry.size(); ++sec)
                        {
                            it_id = NodeRowId.find(rGeometry[sec].Id());
                            std::size_t col = it_id->second;
                            rA(row, col) += Ncontainer(point, prim)*Ncontainer(point, sec) * dV;
                        }
#ifdef _OPENMP
                        omp_unset_lock(&lock_array[row]);
#endif
                    }
                }

                #ifdef ENABLE_BEZIER_GEOMETRY
                rGeometry.Clean();
                #endif
//                ++show_progress;
            }
        }
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    void ConstructMatrixStructure(SparseSpaceType::MatrixType& A,
            const TEntitiesContainerType& pElements, const std::map<std::size_t, std::size_t>& NodeRowId) const
    {
        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);

        typename EntityType::EquationIdVectorType ids;
        std::map<std::size_t, std::size_t>::const_iterator it;
        for(auto i_element = pElements.begin(); i_element != pElements.end() ; ++i_element)
        {
            if( !(i_element)->GetValue( IS_INACTIVE ) || (i_element)->Is(ACTIVE) )
            {
                ids.resize((i_element)->GetGeometry().size());
                for(unsigned int i = 0; i < (i_element)->GetGeometry().size();  ++i)
                {
                    it = NodeRowId.find((i_element)->GetGeometry()[i].Id());
                    if (it == NodeRowId.end())
                        KRATOS_ERROR << "Error: not existing entry " << (i_element)->GetGeometry()[i].Id() << " in NodeRowId";
                    ids[i] = it->second;
                }

                for(std::size_t i = 0 ; i < ids.size() ; ++i)
                {
                    if(ids[i] < equation_size)
                    {
                        std::vector<std::size_t>& row_indices = indices[ids[i]];
                        for(std::size_t j = 0 ; j < ids.size() ; j++)
                        {
                            if(ids[j] < equation_size)
                            {
                                AddUnique(row_indices,ids[j]);
                            }
                        }
                    }
                }
            }
        }

        //allocating the memory needed
        int data_size = 0;
        for(std::size_t i = 0 ; i < indices.size() ; i++)
        {
            data_size += indices[i].size();
        }
        A.reserve(data_size, false);

        //filling with zero the matrix (creating the structure)
#ifndef _OPENMP
        for(std::size_t i = 0 ; i < indices.size() ; i++)
        {
            std::vector<std::size_t>& row_indices = indices[i];
            std::sort(row_indices.begin(), row_indices.end());

            for(std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end() ; it++)
            {
                A.push_back(i,*it,0.00);
            }
            row_indices.clear();
        }
#else
        int number_of_threads = omp_get_max_threads();
        std::vector<unsigned int> matrix_partition;
        OpenMPUtils::CreatePartition(number_of_threads, indices.size(), matrix_partition);
        for( int k=0; k<number_of_threads; ++k )
        {
            #pragma omp parallel
            if( omp_get_thread_num() == k )
            {
                for( std::size_t i = matrix_partition[k]; i < matrix_partition[k+1]; ++i )
                {
                    std::vector<std::size_t>& row_indices = indices[i];
                    std::sort(row_indices.begin(), row_indices.end());

                    for(std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end() ; ++it)
                    {
                        A.push_back(i, *it, 0.00);
                    }

                    row_indices.clear();
                }
            }
        }
#endif
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate) const
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while ( i != endit && (*i) != candidate)
        {
            ++i;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }
    }

}; // Class VariableProjectionUtility

} // namespace Kratos.

#endif /* KRATOS_VARIABLE_PROJECTION_UTILITY_INCLUDED  defined */
