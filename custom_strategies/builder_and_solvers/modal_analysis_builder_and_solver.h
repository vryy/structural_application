/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2008-04-29 12:23:09 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_MODAL_ANALYSIS_BUILDER_AND_SOLVER )
#define  KRATOS_MODAL_ANALYSIS_BUILDER_AND_SOLVER


/* System includes */
#include <set>
#include <omp.h>

/* External includes */


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "linear_solvers/power_iteration_eigenvalue_solver.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.

Detail class definition.

Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


*/
template<class TSparseSpace,
         class TDenseSpace , //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ModalAnalysisBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( ModalAnalysisBuilderAndSolver );


    typedef BuilderAndSolver<TSparseSpace,TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    ModalAnalysisBuilderAndSolver(typename TLinearSolver::Pointer pNewLinearSystemSolver)
    : BaseType(pNewLinearSystemSolver), mMaxEigenSolutions(1), mTolerance(1.0e-8), mMaxIterations(1000)
    {}


    /** Destructor.
    */
    virtual ~ModalAnalysisBuilderAndSolver() {}


    /*@} */
    /**@name Operators
    */
    /*@{ */

    void SetMaxEigenSolutions(std::size_t n) {mMaxEigenSolutions = n;}
    void SetMaxIterations(std::size_t n) {mMaxIterations = n;}
    void SetTolerance(double tol) {mTolerance = tol;}

    std::size_t GetNumberOfEigenSolutions() const {return mEigenvalues.size();}

    TSystemVectorType GetEigenVector(std::size_t i) const
    {
        TSystemVectorType V;
        V.resize(mEigenvectors.size2(), false);

        noalias(V) = row(mEigenvectors, i);

        return V;
    }

    //**************************************************************************
    //**************************************************************************
    virtual void SetUpDofSet( typename TSchemeType::Pointer pScheme, ModelPart& r_model_part )
    {
        KRATOS_TRY

        KRATOS_WATCH("setting up the dofs");
        //Gets the array of elements from the modeler
        ElementsArrayType& pElements = r_model_part.Elements();

        Element::DofsVectorType ElementalDofList;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();
        //mDofSet.clear();

        //double StartTime = GetTickCount();
        for (typename ElementsArrayType::iterator it=pElements.begin();
                it!=pElements.end(); ++it)
        {
            // gets list of Dof involved on every element
            pScheme->GetDofList(*it,ElementalDofList,CurrentProcessInfo);
            for(typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ;
                    i != ElementalDofList.end() ; ++i)
            {
                Doftemp.push_back(*i);
                //mDofSet.push_back(*i);
            }
        }

        //taking in account conditions
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::iterator it=pConditions.begin();
                it!=pConditions.end(); ++it)
        {
            // gets list of Dof involved on every element
            pScheme->GetDofList(*it,ElementalDofList,CurrentProcessInfo);
            for(typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ;
                    i != ElementalDofList.end() ; ++i)
            {
                //mDofSet.push_back(*i);
                Doftemp.push_back(*i);
            }
        }
        Doftemp.Unique();
        BaseType::mDofSet = Doftemp;

        //throws an execption if there are no Degrees of freedom involved in the analysis
        if (BaseType::mDofSet.size()==0)
            KRATOS_THROW_ERROR(std::logic_error, "No degrees of freedom!", "");
        BaseType::mDofSetIsInitialized = true;
        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    virtual void SetUpSystem(
        ModelPart& r_model_part
    )
    {
        // Set equation id for degrees of freedom
        // the free degrees of freedom are positioned at the beginning of the system,
        // while the fixed one are at the end (in opposite order).
        //
        // that means that if the EquationId is greater than "mEquationSystemSize"
        // the pointed degree of freedom is restrained
        //
        int free_id = 0;
        int fix_id = BaseType::mDofSet.size();

        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            if (dof_iterator->IsFixed())
                dof_iterator->SetEquationId(--fix_id);
            else
                dof_iterator->SetEquationId(free_id++);

        BaseType::mEquationSystemSize = fix_id;
    }

    //**************************************************************************
    //**************************************************************************
    virtual void ResizeAndInitializeEigenSystem(
        TSystemMatrixPointerType& pK,
        TSystemMatrixPointerType& pM,
        ElementsArrayType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo
    )
    {
        KRATOS_TRY
        if(pK == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0,0) );
            pK.swap(pNewA);
        }
        if(pM == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewM = TSystemMatrixPointerType(new TSystemMatrixType(0,0) );
            pM.swap(pNewM);
        }
        TSystemMatrixType& K = *pK;
        TSystemMatrixType& M = *pM;

        //resizing the system vectors and matrix
        if (K.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            K.resize(BaseType::mEquationSystemSize,BaseType::mEquationSystemSize,false);
            ConstructMatrixStructure(K,rElements,rConditions,CurrentProcessInfo);
        }
        else
        {
            if(K.size1() != BaseType::mEquationSystemSize || K.size2() != BaseType::mEquationSystemSize)
            {
                KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
                K.resize(BaseType::mEquationSystemSize,BaseType::mEquationSystemSize,true);
                ConstructMatrixStructure(K,rElements,rConditions,CurrentProcessInfo);
            }
        }

        // here we enforce K and M has the same matrix structure
        M = K;

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    virtual void BuildEigenSystem( typename TSchemeType::Pointer pScheme,
                        ModelPart& r_model_part,
                        TSystemMatrixType& K,
                        TSystemMatrixType& M)
    {
        KRATOS_TRY

        boost::timer building_time;

        //build matrices
        BuildSystemMatrices( pScheme, r_model_part, K, M );

        //elapsed time
        if(BaseType::GetEchoLevel() > 0)
        {
            std::cout << "Building Time : " << building_time.elapsed() << std::endl;
        }

        if (BaseType::GetEchoLevel() == 3)
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "stiffness Matrix = " << K << std::endl;
            std::cout << "mass Matrix = " << M << std::endl;
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    virtual void BuildAndSolve( typename TSchemeType::Pointer pScheme,
                        ModelPart& r_model_part,
                        TSystemMatrixType& A,
                        TSystemVectorType& Dx,
                        TSystemVectorType& b )
    {
        KRATOS_TRY

        boost::timer building_time;

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& pConditions = r_model_part.Conditions();

        //construct mass matrix structure
        TSystemMatrixType M = TSystemMatrixType( A.size1(), A.size2() );
        ConstructMatrixStructure( M, pElements, pConditions, r_model_part.GetProcessInfo() );

        //build matrices
        BuildSystemMatrices( pScheme, r_model_part, A, M );

        //elapsed time
        if(BaseType::GetEchoLevel() > 0)
        {
            std::cout << "Building Time : " << building_time.elapsed() << std::endl;
        }

        if (BaseType::GetEchoLevel() == 3)
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "stiffness Matrix = " << A << std::endl;
            std::cout << "mass Matrix = " << M << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        boost::timer solve_time;

        typedef PowerIterationEigenvalueSolver<TSparseSpace, TDenseSpace, TLinearSolver> EigenvalueSolverType;
        EigenvalueSolverType eigenvalue_solver( mTolerance, mMaxIterations, mMaxEigenSolutions, BaseType::mpLinearSystemSolver );

        eigenvalue_solver.Solve( A, M, mEigenvalues, mEigenvectors);

        if(BaseType::GetEchoLevel() > 0)
        {
            std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
        }
        if (BaseType::GetEchoLevel() == 3)
        {
            std::cout << "after the solution of the system" << std::endl;
            std::cout << "Eigenvalues = " << mEigenvalues << std::endl;
            std::cout << "Eigenvectors = " << mEigenvectors << std::endl;
        }

        KRATOS_CATCH("")
    }

    /**
    this function is intended to be called at the end of the solution step to clean up memory
    storage not needed
    */
    void Clear()
    {
        this->mDofSet = DofsArrayType();
        
        if(this->mpReactionsVector != NULL)
        {
            TSparseSpace::Clear( (this->mpReactionsVector) );
        }

        if (this->GetEchoLevel()>0)
        {
            std::cout << "ModalAnalysisBuilderAndSolver Clear Function called" << std::endl;
        }
    }

    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */
    //**************************************************************************
    virtual void ConstructMatrixStructure(
        TSystemMatrixType& A,
        ElementsContainerType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo)
    {

        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);
        //              std::vector<std::vector<std::size_t> > dirichlet_indices(TSystemSpaceType::Size1(mDirichletMatrix));

        Element::EquationIdVectorType ids(3,0);
        for(typename ElementsContainerType::iterator i_element = rElements.begin() ; i_element != rElements.end() ; i_element++)
        {
            (i_element)->EquationIdVector(ids, CurrentProcessInfo);

            for(std::size_t i = 0 ; i < ids.size() ; i++)
                if(ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for(std::size_t j = 0 ; j < ids.size() ; j++)
                        if(ids[j] < equation_size)
                        {
                            AddUnique(row_indices,ids[j]);
                            //indices[ids[i]].push_back(ids[j]);
                        }
                }
        }

        for(typename ConditionsArrayType::iterator i_condition = rConditions.begin() ; i_condition != rConditions.end() ; i_condition++)
        {
            (i_condition)->EquationIdVector(ids, CurrentProcessInfo);
            for(std::size_t i = 0 ; i < ids.size() ; i++)
                if(ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for(std::size_t j = 0 ; j < ids.size() ; j++)
                        if(ids[j] < equation_size)
                        {
                            AddUnique(row_indices,ids[j]);
                            //  indices[ids[i]].push_back(ids[j]);
                        }
                }
        }

        //allocating the memory needed
        int data_size = 0;
        for(std::size_t i = 0 ; i < indices.size() ; i++)
        {
            data_size += indices[i].size();
        }
        A.reserve(data_size,false);

        //filling with zero the matrix (creating the structure)
        for(std::size_t i = 0 ; i < indices.size() ; i++)
        {
            std::vector<std::size_t>& row_indices = indices[i];
            std::sort(row_indices.begin(), row_indices.end());

            for(std::vector<std::size_t>::iterator it= row_indices.begin(); it != row_indices.end() ; it++)
            {
                A.push_back(i,*it,0.00);
//                  A()(i,*it) = 0.00;
            }
            //row_indices = std::vector<std::size_t>();
            row_indices.clear();
        }
    }

    //**************************************************************************
    void AssembleLHS(
        TSystemMatrixType& A,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local=0; i_local<local_size; i_local++)
        {
            unsigned int i_global=EquationId[i_local];
            if ( i_global < BaseType::mEquationSystemSize )
            {
                for (unsigned int j_local=0; j_local<local_size; j_local++)
                {
                    unsigned int j_global=EquationId[j_local];
                    if ( j_global < BaseType::mEquationSystemSize )
                    {
                        A(i_global,j_global) += LHS_Contribution(i_local,j_local);
                    }
                }
            }
        }
    }



    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    std::size_t mMaxEigenSolutions;
    std::size_t mMaxIterations;
    double mTolerance;

    LocalSystemVectorType mEigenvalues;
    LocalSystemMatrixType mEigenvectors;

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */

    //**************************************************************************
    //**************************************************************************
    void BuildSystemMatrices(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& K,
        TSystemMatrixType& M )
    {
        KRATOS_TRY
        if(!pScheme)
            KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //create a partition of the element array
        int number_of_threads = omp_get_max_threads();
        std::vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);

        double start_prod = omp_get_wtime();

        #pragma omp parallel for
        for(int k=0; k<number_of_threads; k++)
        {
            //contributions to the system
            LocalSystemMatrixType K_Contribution = LocalSystemMatrixType(0,0);
            LocalSystemMatrixType M_Contribution = LocalSystemMatrixType(0,0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            //vector containing the localization in the system of the different
            //terms
            Element::EquationIdVectorType EquationId;
            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
            typename ElementsArrayType::iterator
            it_begin=pElements.begin()+element_partition[k];
            typename ElementsArrayType::iterator
            it_end=pElements.begin()+element_partition[k+1];

            // assemble all elements
            for (typename ElementsArrayType::iterator it=it_begin; it!=it_end; ++it)
            {
                //calculate elemental contribution
                pScheme->CalculateSystemContributions(*it, K_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                it->CalculateMassMatrix( M_Contribution, CurrentProcessInfo );
                //(*it)->CalculateLocalSystem( K_Contribution,RHS_Contribution,CurrentProcessInfo );

                #pragma omp critical
                {
                    //assemble the elemental contribution
                    AssembleLHS(K,K_Contribution,EquationId);
                    AssembleLHS(M,M_Contribution,EquationId);
                    // clean local elemental memory
                    pScheme->CleanMemory(*it);
                }
            }
        }

        std::vector<unsigned int> condition_partition;
        CreatePartition(number_of_threads, ConditionsArray.size(), condition_partition);

        #pragma omp parallel for
        for(int k=0; k<number_of_threads; k++)
        {
            //contributions to the system
            LocalSystemMatrixType K_Contribution = LocalSystemMatrixType(0,0);
            LocalSystemMatrixType M_Contribution = LocalSystemMatrixType(0,0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            Condition::EquationIdVectorType EquationId;

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            typename ConditionsArrayType::iterator
            it_begin=ConditionsArray.begin()+condition_partition[k];
            typename ConditionsArrayType::iterator
            it_end=ConditionsArray.begin()+condition_partition[k+1];

            // assemble all conditions
            for (typename ConditionsArrayType::iterator it=it_begin; it!=it_end; ++it)
            {
                //calculate elemental contribution
                pScheme->CalculateSystemContributions(*it, K_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                it->CalculateMassMatrix( M_Contribution, CurrentProcessInfo );
                //(*it)->CalculateLocalSystem( K_Contribution,RHS_Contribution,CurrentProcessInfo );

                #pragma omp critical
                {
                    //assemble the elemental contribution
                    AssembleLHS(K,K_Contribution,EquationId);
                    AssembleLHS(M,M_Contribution,EquationId);
                }
            }
        }
        double stop_prod = omp_get_wtime();
        std::cout << "time: " << stop_prod - start_prod << std::endl;
        KRATOS_WATCH("finished parallel building");

        KRATOS_CATCH("")
    }

    //******************************************************************************************
    //******************************************************************************************
    inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while ( i != endit && (*i) != candidate)
        {
            i++;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }

    }

    //******************************************************************************************
    //******************************************************************************************
    inline void CreatePartition(unsigned int number_of_threads,const int number_of_rows, std::vector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads+1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for(int i = 1; i<number_of_threads; i++)
            partitions[i] = partitions[i-1] + partition_size ;
    }

    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class ModalAnalysisBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_MODAL_ANALYSIS_BUILDER_AND_SOLVER  defined */

