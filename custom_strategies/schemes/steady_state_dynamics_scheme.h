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
*   Date:                $Date: 18/3/2025 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/


#if !defined(KRATOS_STRUCTURAL_APPLICATION_STEADY_STATE_DYNAMICS_SCHEME_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APPLICATION_STEADY_STATE_DYNAMICS_SCHEME_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "includes/fnv_1a_hash.h"
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include "custom_python3/legacy_structural_app_vars.h"
#else
#include "includes/legacy_structural_app_vars.h"
#endif
#include "solving_strategies/schemes/scheme.h"
#include "structural_application_variables.h"

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
/**
 * Steady state dynamics scheme for analysis in frequency domain
 * Reference: tn6.pdf
 */
template<class TSparseSpace, class TDenseSpace, class TModelPartType = ModelPart>
class SteadyStateDynamicsScheme : public Scheme<TSparseSpace, TDenseSpace, TModelPartType>
{
public:
    /**@name Type Definitions */
    KRATOS_CLASS_POINTER_DEFINITION( SteadyStateDynamicsScheme );

    typedef Scheme<TSparseSpace, TDenseSpace, TModelPartType> BaseType;

    typedef typename BaseType::ModelPartType ModelPartType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename DataTypeToValueType<TDataType>::value_type ValueType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::ElementType ElementType;

    typedef typename BaseType::ConditionType ConditionType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    /**
     * Constructor
     */
    SteadyStateDynamicsScheme(const ValueType omega) : BaseType(), mOmega(omega)
    {
        std::cout << "Steady State Dynamics Scheme!!!!!!!!!!!!!!!!!!!!!Omega = " << mOmega << std::endl;
    }

    /** Destructor.*/
    ~SteadyStateDynamicsScheme() override
    {}

    /**@name Operators */

    /// Copy assignment
    SteadyStateDynamicsScheme& operator=(const SteadyStateDynamicsScheme& rOther)
    {
        mOmega = rOther.mOmega;
        return *this;
    }

    /*@} */
    /**@name Operations */
    /*@{ */

    ValueType GetOmega() const
    {
        return mOmega;
    }

    void SetOmega(ValueType omega)
    {
        mOmega = omega;
    }

    /// Initialize the scheme
    void Initialize(ModelPartType& r_model_part) override
    {
        BaseType::Initialize(r_model_part);

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        CurrentProcessInfo[TIME_INTEGRATION_SCHEME] = FNV1a32Hash::CalculateHash(Info().c_str());

        KRATOS_WATCH(CurrentProcessInfo[TIME_INTEGRATION_SCHEME])
    }

    void Update(
        ModelPartType& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY


        for (typename DofsArrayType::iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            typename ModelPartType::NodeType& rNode = r_model_part.GetNode(dof_iterator->Id());

            if (dof_iterator->GetVariable() == VARSELC(TDataType, DISPLACEMENT, X))
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(VARSELC(TDataType, DISPLACEMENT, X))
                    += Dx[dof_iterator->EquationId()];
                    rNode.GetSolutionStepValue(VARSELC(TDataType, ACCELERATION, X))
                    -= mOmega*mOmega*Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == VARSELC(TDataType, DISPLACEMENT, Y))
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(VARSELC(TDataType, DISPLACEMENT, Y))
                    += Dx[dof_iterator->EquationId()];
                    rNode.GetSolutionStepValue(VARSELC(TDataType, ACCELERATION, Y))
                    -= mOmega*mOmega*Dx[dof_iterator->EquationId()];
                }
            }
            else if (dof_iterator->GetVariable() == VARSELC(TDataType, DISPLACEMENT, Z))
            {
                if (dof_iterator->IsFree())
                {
                    rNode.GetSolutionStepValue(VARSELC(TDataType, DISPLACEMENT, Z))
                    += Dx[dof_iterator->EquationId()];
                    rNode.GetSolutionStepValue(VARSELC(TDataType, ACCELERATION, Z))
                    -= mOmega*mOmega*Dx[dof_iterator->EquationId()];
                }
            }
        }

        KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************

    void CalculateSystemContributions(
        ElementType& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename ElementType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        this->CalculateSystemContributionsImpl(rCurrentElement,
                LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(
        ConditionType& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename ConditionType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        this->CalculateSystemContributionsImpl(rCurrentCondition,
                LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

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
        return "SteadyStateDynamicsScheme";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Omega: " << mOmega;
    }

    /*@} */
    /**@name Friends */
    /*@{ */

protected:
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

private:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */
    ValueType mOmega;
    /*@} */
    /**@name Private Operators*/
    /*@{ */
    /*@} */
    /**@name Private Operations*/
    /*@{ */

    template<typename TEntityType>
    inline void CalculateSystemContributionsImpl(
        TEntityType& rCurrentElement,
        typename TEntityType::MatrixType& LHS_Contribution,
        typename TEntityType::VectorType& RHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) const
    {
        KRATOS_TRY

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        typename TEntityType::VectorType acceleration;
        rCurrentElement.GetSecondDerivativesVector(acceleration, 0);

        if (acceleration.size() > 0)
        {
            typename TEntityType::MatrixType MassMatrix;

            rCurrentElement.CalculateMassMatrix(MassMatrix, CurrentProcessInfo);

            noalias(RHS_Contribution) -= prod(MassMatrix, acceleration);

            noalias(LHS_Contribution) -= mOmega*mOmega*MassMatrix;
        }

        KRATOS_CATCH("")
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
}; /* class SteadyStateDynamicsScheme */

}  /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APPLICATION_STEADY_STATE_DYNAMICS_SCHEME_H_INCLUDED  defined */
