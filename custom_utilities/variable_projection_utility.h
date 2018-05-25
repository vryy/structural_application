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
 * Utility to project the variables from Gauss points to nodes using L2-projection.
 */
class VariableProjectionUtility
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
    VariableProjectionUtility(LinearSolverType::Pointer pLinearSolver)
    : mpLinearSolver(pLinearSolver), mEchoLevel(-1)
    {
        std::cout << "VariableProjectionUtility created" << std::endl;
    }

    /**
     * Destructor.
     */
    virtual ~VariableProjectionUtility()
    {}

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

    /// Initialize the projection matrix
    void BeginProjection(ElementsContainerType& pElements)
    {
        NodesContainerType pActiveNodes;
        std::map<std::size_t, std::size_t> NodeRowId;        
        this->ExtractActiveNodes(pElements, pActiveNodes, NodeRowId);
        this->ConstructLHSMatrix(mProjectionMatrix, pElements, NodeRowId);
    }

    /**
     * Transfer the double variable from Gauss point to nodes
     * Here a different element set is allowed to compute the right hand side. This element set shall have the same connectivities as the element set used in InitializeProjection.
     */
    void TransferVariablesToNodes(Variable<double>& rThisVariable,
            ElementsContainerType& pElements, const ProcessInfo& rProcessInfo)
    {
        NodesContainerType pActiveNodes;
        std::map<std::size_t, std::size_t> NodeRowId;        
        this->ExtractActiveNodes(pElements, pActiveNodes, NodeRowId);

        //SetUpEquationSystem
        SparseSpaceType::VectorType g(pActiveNodes.size());
        SparseSpaceType::VectorType b(pActiveNodes.size());

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

        //Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization
        // see Jiao + Heath "Common-refinement-based data tranfer ..."
        // International Journal for numerical methods in engineering 61 (2004) 2402--2427
        // for general description of L_2-Minimization

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);

#ifdef _OPENMP
            #pragma omp parallel for
#endif
        for (int k = 0; k < number_of_threads; ++k)
        {
            ElementsContainerType::ptr_iterator it_begin = pElements.ptr_begin() + element_partition[k];
            ElementsContainerType::ptr_iterator it_end = pElements.ptr_begin() + element_partition[k+1];

            for ( ElementsContainerType::ptr_iterator it = it_begin; it != it_end; ++it )
            {
                if ( ! ( ( (*it)->GetValue(IS_INACTIVE) == false ) || (*it)->Is(ACTIVE) ) )
                    continue;

                const IntegrationPointsArrayType& integration_points
                    = (*it)->GetGeometry().IntegrationPoints( (*it)->GetIntegrationMethod());

                GeometryType::JacobiansType J(integration_points.size());
                J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());

                std::vector<double> ValuesOnIntPoint(integration_points.size());
                (*it)->GetValueOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, rProcessInfo);

                const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                double DetJ;
                for (unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    DetJ = MathUtils<double>::Det(J[point]);

                    double dV = DetJ*integration_points[point].Weight();

                    for (unsigned int prim = 0; prim<(*it)->GetGeometry().size(); ++prim)
                    {
                        unsigned int row = NodeRowId[(*it)->GetGeometry()[prim].Id()];
#ifdef _OPENMP
                        omp_set_lock(&lock_array[row]);
#endif
                        b(row) += ValuesOnIntPoint[point] * Ncontainer(point, prim) * dV;
#ifdef _OPENMP
                        omp_unset_lock(&lock_array[row]);
#endif
                    }
                }
            }
        }
/*KRATOS_WATCH(b)*/
        mpLinearSolver->Solve(mProjectionMatrix, g, b);
/*KRATOS_WATCH(g)*/
        for(NodesContainerType::iterator it = pActiveNodes.begin() ; it != pActiveNodes.end() ; it++)
        {
            it->GetSolutionStepValue(rThisVariable) = g(NodeRowId[it->Id()]);
        }

#ifdef _OPENMP
        for(unsigned int i = 0; i < M_size; ++i)
            omp_destroy_lock(&lock_array[i]);
#endif
        std::cout << "TransferVariablesToNodes for " << rThisVariable.Name() << " completed" << std::endl;
    }

    void EndProjection()
    {
        mProjectionMatrix.resize(0, 0, false);
    }

protected:

    LinearSolverType::Pointer mpLinearSolver;

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    /// Extract the active nodes and construct the nodal id row map
    void ExtractActiveNodes(ElementsContainerType& pElements,
        NodesContainerType& pActiveNodes, std::map<std::size_t, std::size_t>& NodeRowId) const
    {
        // extract the active nodes
        for(ModelPart::ElementsContainerType::ptr_iterator it = pElements.ptr_begin();
                it != pElements.ptr_end(); ++it)
        {
            if( (*it)->GetValue(IS_INACTIVE) == false || (*it)->Is(ACTIVE) )
            {
                for(std::size_t i = 0; i < (*it)->GetGeometry().size(); ++i)
                {
                    pActiveNodes.push_back( (*it)->GetGeometry().pGetPoint(i) );
                }
            }
        }

        pActiveNodes.Unique();
        KRATOS_WATCH(pActiveNodes.size())

        // assign each node an id. That id is the row of this node in the global L2 projection matrix*/
        std::size_t cnt = 0;
        for(NodesContainerType::iterator it = pActiveNodes.begin(); it != pActiveNodes.end(); ++it )
        {
            NodeRowId[it->Id()] = cnt++;
        }
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    /// Construct the L2 projection matrix
    void ConstructLHSMatrix(SparseSpaceType::MatrixType& rA, ElementsContainerType& pElements,
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
        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
//        boost::progress_display show_progress( pElements.size() );

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
            ModelPart::ElementsContainerType::ptr_iterator it_begin = pElements.ptr_begin() + element_partition[k];
            ModelPart::ElementsContainerType::ptr_iterator it_end = pElements.ptr_begin() + element_partition[k+1];
            std::map<std::size_t, std::size_t>::const_iterator it_id;

            for( ModelPart::ElementsContainerType::ptr_iterator it = it_begin; it != it_end; ++it )
            {
                if( (*it)->GetValue(IS_INACTIVE) == true && !(*it)->Is(ACTIVE) )
                    continue;

                unsigned int dim = (*it)->GetGeometry().WorkingSpaceDimension();

                const GeometryType::IntegrationPointsArrayType& integration_points
                    = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                GeometryType::JacobiansType J(integration_points.size());
                J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());

                const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                double DetJ;
                for(unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    DetJ = MathUtils<double>::Det(J[point]);

                    double dV = DetJ*integration_points[point].Weight();

                    for(unsigned int prim = 0 ; prim < (*it)->GetGeometry().size(); ++prim)
                    {
                        it_id = NodeRowId.find((*it)->GetGeometry()[prim].Id());
                        std::size_t row = it_id->second;
#ifdef _OPENMP
                        omp_set_lock(&lock_array[row]);
#endif
                        for(unsigned int sec = 0 ; sec < (*it)->GetGeometry().size(); ++sec)
                        {
                            it_id = NodeRowId.find((*it)->GetGeometry()[sec].Id());
                            std::size_t col = it_id->second;
                            rA(row, col) += Ncontainer(point, prim)*Ncontainer(point, sec) * dV;
                        }
#ifdef _OPENMP
                        omp_unset_lock(&lock_array[row]);
#endif
                    }
                }
//                ++show_progress;
            }
        }
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    void ConstructMatrixStructure(SparseSpaceType::MatrixType& A,
            ElementsContainerType& pElements, const std::map<std::size_t, std::size_t>& NodeRowId) const
    {
        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);

        Element::EquationIdVectorType ids;
        std::map<std::size_t, std::size_t>::const_iterator it;
        for(ModelPart::ElementsContainerType::iterator i_element = pElements.begin();
                i_element != pElements.end() ; ++i_element)
        {
            if( !(i_element)->GetValue( IS_INACTIVE ) || (i_element)->Is(ACTIVE) )
            {
                ids.resize((i_element)->GetGeometry().size());
                for(unsigned int i = 0; i < (i_element)->GetGeometry().size();  ++i)
                {
                    it = NodeRowId.find((i_element)->GetGeometry()[i].Id());
                    if (it == NodeRowId.end())
                        KRATOS_THROW_ERROR(std::logic_error, "Error: not existing entry in NodeRowId", "")
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
        vector<unsigned int> matrix_partition;
        CreatePartition(number_of_threads, indices.size(), matrix_partition);
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

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    inline void CreatePartition(unsigned int number_of_threads,const int number_of_rows, vector<unsigned int>& partitions) const
    {
        partitions.resize(number_of_threads+1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for(unsigned int i = 1; i<number_of_threads; i++)
            partitions[i] = partitions[i-1] + partition_size ;
    }

private:

    SparseSpaceType::MatrixType mProjectionMatrix;
    int mEchoLevel;

};//Class VariableProjectionUtility

}//namespace Kratos.

#endif /* KRATOS_VARIABLE_PROJECTION_UTILITY_INCLUDED  defined */
