/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 7 Mar 2022 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_ARC_LENGTH_STATIC_SCHEME_H_INCLUDED )
#define  KRATOS_RESIDUALBASED_ARC_LENGTH_STATIC_SCHEME_H_INCLUDED


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
    virtual ~ResidualBasedIncrementalUpdateStaticDeactivationScheme() {}


    /*@} */
    /**@name Operators
    */
    /*@{ */

    /**
    Performing the update of the solution.
    */
    //***************************************************************************
    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) final
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

    /**
    Function called once at the beginning of each solution step.
    The basic operations to be carried in there are the following:
    - managing variables to be kept constant over the time step
    (for example time-Scheme constants depending on the actual time step)
     */
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) final
    {
        KRATOS_TRY
        //initialize solution step for all of the elements
        ElementsArrayType& pElements = r_model_part.Elements();
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        bool element_is_active;
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            element_is_active = true;
            if(it->IsDefined(ACTIVE))
                element_is_active = it->Is(ACTIVE);

            if (element_is_active)
                (it)->InitializeSolutionStep(CurrentProcessInfo);
        }

        bool condition_is_active;
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            condition_is_active = true;
            if(it->IsDefined(ACTIVE))
                condition_is_active = it->Is(ACTIVE);

            if (condition_is_active)
                (it)->InitializeSolutionStep(CurrentProcessInfo);
        }
        KRATOS_CATCH("")
    }

    /**
    function to be called when it is needed to initialize an iteration.
    it is designed to be called at the beginning of each non linear iteration

      take care: the elemental function with the same name is NOT called here.
      The function is called in the builder for memory efficiency
     */
    void InitializeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) final
    {
        KRATOS_TRY
        ElementsArrayType& pElements = r_model_part.Elements();
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

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

    /**
    function to be called when it is needed to finalize an iteration.
    it is designed to be called at the end of each non linear iteration
     */
    void FinalizeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) final
    {
        KRATOS_TRY

        ElementsArrayType& pElements = r_model_part.Elements();
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

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

        KRATOS_CATCH("")
    }

    /**
    function called once at the end of a solution step, after convergence is reached if
    an iterative process is needed
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) final
    {
        KRATOS_TRY
        //finalizes solution step for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

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

    /** this function is designed to be called in the builder and solver to introduce
    the selected time integration scheme. It "asks" the matrix needed to the element and
    performs the operations needed to introduce the seected time integration scheme.

    this function calculates at the same time the contribution to the LHS and to the RHS
    of the system
    */
    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) final
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


    //***************************************************************************
    //***************************************************************************
    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) final
    {
        KRATOS_TRY

        rCurrentElement.CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************
    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) final
    {
        KRATOS_TRY

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    /** functions totally analogous to the precedent but applied to
    the "condition" objects
    */
    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) final
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
    ) final
    {
        KRATOS_TRY

        rCurrentCondition.CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
        rCurrentCondition.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
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

#endif /* KRATOS_RESIDUALBASED_ARC_LENGTH_STATIC_SCHEME_H_INCLUDED  defined */

