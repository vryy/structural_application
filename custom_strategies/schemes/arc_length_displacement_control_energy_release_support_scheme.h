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
*   Date:                $Date: 23/3/2022 $
*   Revision:            $Revision: 1.7 $
*
* ***********************************************************/


#if !defined(KRATOS_STRUCTURAL_APPLICATION_ARC_LENGTH_DISPLACEMENT_CONTROL_ENERGY_RELEASE_SUPPORT_SCHEME )
#define  KRATOS_STRUCTURAL_APPLICATION_ARC_LENGTH_DISPLACEMENT_CONTROL_ENERGY_RELEASE_SUPPORT_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include "custom_python3/legacy_structural_app_vars.h"
#else
#include "includes/legacy_structural_app_vars.h"
#endif
#include "custom_strategies/schemes/arc_length_displacement_control_support_scheme.h"
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
 * This scheme is a host scheme. This is used mainly with arc-length energy release control strategy for displacement control.
 * This scheme shall behave identically
 * as the underlying scheme when another solving strategy is used
 * Reference:
 * +    Verhoosel et al, A dissipation-based arc-length method for robust simulation of brittle and ductile failure, 2009
*/
template<class TSchemeType>
class ArcLengthDisplacementControlEnergyReleaseSupportScheme : public ArcLengthDisplacementControlSupportScheme<TSchemeType>
{
public:
    /**@name Type Definitions */
    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthDisplacementControlEnergyReleaseSupportScheme );

    typedef ArcLengthDisplacementControlSupportScheme<TSchemeType> BaseType;

    typedef typename BaseType::SparseSpaceType SparseSpaceType;
    typedef typename BaseType::DenseSpaceType DenseSpaceType;

    typedef typename SparseSpaceType::SizeType SizeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::ElementsArrayType ElementsArrayType;

    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    /**
     * Default constructor
     */
    ArcLengthDisplacementControlEnergyReleaseSupportScheme(typename TSchemeType::Pointer pScheme)
    : BaseType(pScheme)
    {
        std::cout << "ArcLengthDisplacementControlEnergyReleaseSupportScheme is invoked on " << pScheme->Info() << std::endl;
        BaseType::mName = "ArcLengthDisplacementControlEnergyReleaseSupportScheme<" + pScheme->Info() + ">";
    }

    /** Destructor.*/
    virtual ~ArcLengthDisplacementControlEnergyReleaseSupportScheme()
    {}

    /**@name Operators */

    /*@} */
    /**@name Operations */
    /*@{ */

    const TSystemVectorType& GetForceVector2() const
    {
        return mLastForceVector2;
    }

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        BaseType::InitializeSolutionStep(r_model_part, A, Dx, b);

        if (SparseSpaceType::Size(mCurrentForceVector2) != BaseType::mEquationSystemSize)
            SparseSpaceType::Resize(mCurrentForceVector2, BaseType::mEquationSystemSize);
        SparseSpaceType::SetToZero(mCurrentForceVector2);
    }

    void InitializeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        BaseType::InitializeNonLinIteration(r_model_part, A, Dx, b);

        SparseSpaceType::SetToZero(mCurrentForceVector2);
    }

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        BaseType::CalculateSystemContributions(rCurrentElement, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

        // KRATOS_WATCH(LHS_Contribution)
        // KRATOS_WATCH(RHS_Contribution)

        // compute the forces induced by displacement
        LocalSystemVectorType Force;
        DenseSpaceType::Resize(Force, DenseSpaceType::Size(RHS_Contribution));
        DenseSpaceType::SetToZero(Force);

        // assemble the force vector
        Element::DofsVectorType ElementalDofList;
        this->GetDofList(rCurrentElement, ElementalDofList, CurrentProcessInfo);

        std::vector<std::size_t> local_ids;
        std::vector<std::size_t> global_ids;
        std::size_t cnt = 0;
        for (Element::DofsVectorType::iterator it = ElementalDofList.begin(); it != ElementalDofList.end(); ++it)
        {
            if ((*it)->GetVariable() == DISPLACEMENT_X
             || (*it)->GetVariable() == DISPLACEMENT_Y
             || (*it)->GetVariable() == DISPLACEMENT_Z )
            {
                local_ids.push_back(cnt);
                global_ids.push_back((*it)->EquationId());
            }

            ++cnt;
        }

        // std::cout << "local_ids:";
        // for (std::size_t i = 0; i < local_ids.size(); ++i)
        //     std::cout << " " << local_ids[i];
        // std::cout << std::endl;

        // std::cout << "global_ids:";
        // for (std::size_t i = 0; i < global_ids.size(); ++i)
        //     std::cout << " " << global_ids[i];
        // std::cout << std::endl;

        std::size_t disp_size = local_ids.size();
        LocalSystemVectorType LocalValues;
        rCurrentElement.GetValuesVector(LocalValues);
        // KRATOS_WATCH(LocalValues)
        // LocalSystemVectorType Tmp1 = -RHS_Contribution - prod(LHS_Contribution, LocalValues);
        // KRATOS_WATCH(Tmp1)
        // LocalSystemVectorType Tmp2(18, 0.0);
        // for (int i = 0; i < 18; ++i)
        // {
        //     Tmp2[i] = -RHS_Contribution[i];
        //     for (int j = 0; j < 18; ++j)
        //         Tmp2[i] -= LHS_Contribution(i, j)*LocalValues[j];
        // }
        // KRATOS_WATCH(Tmp2)

        for (std::size_t i = 0; i < disp_size; ++i)
        {
            if (global_ids[i] < BaseType::mEquationSystemSize)
            {
                double v = -RHS_Contribution[local_ids[i]]; // assuming RHS only contain internal forces

                for (std::size_t j = 0; j < disp_size; ++j)
                    v -= LHS_Contribution(local_ids[i], local_ids[j]) * LocalValues[local_ids[j]];

                // KRATOS_WATCH(v)

                mCurrentForceVector2[global_ids[i]] += v;
            }
        }
    }

    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        BaseType::FinalizeSolutionStep(r_model_part, A, Dx, b);

        SparseSpaceType::Copy(mCurrentForceVector2, mLastForceVector2);
        double tmp = SparseSpaceType::TwoNorm(mLastForceVector2);
        std::cout << "Norm of the force vector 2: " << tmp << std::endl;
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


    /*@} */
    /**@name Friends */
    /*@{ */

private:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */

    TSystemVectorType mCurrentForceVector2;
    TSystemVectorType mLastForceVector2;

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
}; /* Class Scheme */

}  /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APPLICATION_ARC_LENGTH_DISPLACEMENT_CONTROL_ENERGY_RELEASE_SUPPORT_SCHEME  defined */
