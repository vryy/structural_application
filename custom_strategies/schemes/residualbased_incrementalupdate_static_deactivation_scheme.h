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
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 20 Apr 2020 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_INCREMENTAL_STATIC_DEACTIVATION_SCHEME )
#define  KRATOS_RESIDUALBASED_INCREMENTAL_STATIC_DEACTIVATION_SCHEME


/* System includes */


/* External includes */


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "custom_elements/prescribed_object.h"

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

This class provides the implementation of the basic tasks that are needed by the solution strategy.
It is intended to be the place for tailoring the solution strategies to problem specific tasks.

Detail class definition.

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
         class TDenseSpace //= DenseSpace<double>
         >
class ResidualBasedIncrementalUpdateStaticDeactivationScheme : public Scheme<TSparseSpace,TDenseSpace>
{

public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedIncrementalUpdateStaticDeactivationScheme);

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef TSparseSpace SparseSpaceType;
    typedef TDenseSpace DenseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    ResidualBasedIncrementalUpdateStaticDeactivationScheme()
        : Scheme<TSparseSpace,TDenseSpace>()
    {}

    /** Destructor.
    */
    ~ResidualBasedIncrementalUpdateStaticDeactivationScheme() override {}

    /*@} */
    /**@name Operators
    */
    /*@{ */

    /// Copy assignment
    ResidualBasedIncrementalUpdateStaticDeactivationScheme& operator=(const ResidualBasedIncrementalUpdateStaticDeactivationScheme& rOther)
    {
        return *this;
    }

    /*@} */
    /**@name Operations */
    /*@{ */

    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY

        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
            }
        }

        KRATOS_CATCH("")
    }

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY

        //initialize solution step for all of the elements
        ElementsArrayType& pElements = r_model_part.Elements();
        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        bool element_is_active;
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            element_is_active = true;
            if(it->IsDefined(ACTIVE))
                element_is_active = it->Is(ACTIVE);

            if (element_is_active)
                it->InitializeSolutionStep(CurrentProcessInfo);
        }

        bool condition_is_active;
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            condition_is_active = true;
            if(it->IsDefined(ACTIVE))
                condition_is_active = it->Is(ACTIVE);

            if (condition_is_active)
                it->InitializeSolutionStep(CurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

    void InitializeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY

        ElementsArrayType& pElements = r_model_part.Elements();
        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        bool element_is_active;
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            element_is_active = true;
            if(it->IsDefined(ACTIVE))
                element_is_active = it->Is(ACTIVE);

            if (element_is_active)
                (it)->InitializeNonLinearIteration(CurrentProcessInfo);
        }

        bool condition_is_active;
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            condition_is_active = true;
            if(it->IsDefined(ACTIVE))
                condition_is_active = it->Is(ACTIVE);

            if (condition_is_active)
                (it)->InitializeNonLinearIteration(CurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

    void FinalizeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY

        // to account for prescribed displacement, the displacement at prescribed nodes need to be updated
        double curr_disp, delta_disp;
        for (ModelPart::NodesContainerType::iterator it_node = r_model_part.Nodes().begin(); it_node != r_model_part.Nodes().end(); ++it_node)
        {
            if (it_node->IsFixed(DISPLACEMENT_X))
            {
                curr_disp = it_node->GetSolutionStepValue(DISPLACEMENT_X);
                delta_disp = it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
                it_node->GetSolutionStepValue(DISPLACEMENT_X) = curr_disp + delta_disp;
                it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X) = 0.0; // set the prescribed displacement to zero to avoid update in the second step
            }

            if (it_node->IsFixed(DISPLACEMENT_Y))
            {
                curr_disp = it_node->GetSolutionStepValue(DISPLACEMENT_Y);
                delta_disp = it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
                it_node->GetSolutionStepValue(DISPLACEMENT_Y) = curr_disp + delta_disp;
                it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y) = 0.0; // set the prescribed displacement to zero to avoid update in the second step
            }

            if (it_node->IsFixed(DISPLACEMENT_Z))
            {
                curr_disp = it_node->GetSolutionStepValue(DISPLACEMENT_Z);
                delta_disp = it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                it_node->GetSolutionStepValue(DISPLACEMENT_Z) = curr_disp + delta_disp;
                it_node->GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z) = 0.0; // set the prescribed displacement to zero to avoid update in the second step
            }
        }

        // invoking the element and condition finalization after an iteration

        ElementsArrayType& pElements = r_model_part.Elements();
        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        bool element_is_active;
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            element_is_active = true;
            if(it->IsDefined(ACTIVE))
                element_is_active = it->Is(ACTIVE);

            if (element_is_active)
                (it)->FinalizeNonLinearIteration(CurrentProcessInfo);
        }

        bool condition_is_active;
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            condition_is_active = true;
            if(it->IsDefined(ACTIVE))
                condition_is_active = it->Is(ACTIVE);

            if (condition_is_active)
                (it)->FinalizeNonLinearIteration(CurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY
        //finalizes solution step for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            typename ElementsArrayType::iterator ElementsBegin = rElements.begin() + ElementPartition[k];
            typename ElementsArrayType::iterator ElementsEnd = rElements.begin() + ElementPartition[k + 1];

            bool element_is_active;
            for (typename ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                element_is_active = true;
                if(itElem->IsDefined(ACTIVE))
                    element_is_active = itElem->Is(ACTIVE);

                if (element_is_active)
                    itElem->FinalizeSolutionStep(CurrentProcessInfo);
            }
        }

        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            typename ConditionsArrayType::iterator ConditionsBegin = rConditions.begin() + ConditionPartition[k];
            typename ConditionsArrayType::iterator ConditionsEnd = rConditions.begin() + ConditionPartition[k + 1];

            bool condition_is_active;
            for (typename ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; itCond++)
            {
                condition_is_active = true;
                if(itCond->IsDefined(ACTIVE))
                    condition_is_active = itCond->Is(ACTIVE);

                if (condition_is_active)
                    itCond->FinalizeSolutionStep(CurrentProcessInfo);
            }
        }

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        try
        {
            const PrescribedObject& rObject = dynamic_cast<PrescribedObject&>(rCurrentElement);
            rObject.ApplyPrescribedDofs(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        }
        catch (std::bad_cast& bc)
        {
            // DO NOTHING
        }

        KRATOS_CATCH("")
    }


    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        rCurrentCondition.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
        rCurrentCondition.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedIncrementalUpdateStaticDeactivationScheme";
    }

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

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


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

}; /* Class Scheme */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_INCREMENTAL_STATIC_DEACTIVATION_SCHEME  defined */
