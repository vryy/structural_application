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
*   Date:                $Date: 17/3/2022 $
*   Revision:            $Revision: 1.7 $
*
* ***********************************************************/


#if !defined(KRATOS_STRUCTURAL_APPLICATION_ARC_LENGTH_DISPLACEMENT_CONTROL_SUPPORT_SCHEME )
#define  KRATOS_STRUCTURAL_APPLICATION_ARC_LENGTH_DISPLACEMENT_CONTROL_SUPPORT_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
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
 * This scheme is a host scheme. This is used mainly with arc-length displacement control strategy
 * to compute the external forces induced by prescribed displacement. This scheme shall behave identically
 * as the underlying scheme when another solving strategy is used
*/
template<class TSchemeType>
class ArcLengthDisplacementControlSupportScheme : public TSchemeType
{
public:
    /**@name Type Definitions */
    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthDisplacementControlSupportScheme );

    typedef TSchemeType BaseType;

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
    ArcLengthDisplacementControlSupportScheme(typename BaseType::Pointer pScheme)
    : BaseType(), mEquationSystemSize(0)
    {
        BaseType::operator=(*pScheme);
        std::cout << "ArcLengthDisplacementControlSupportScheme is invoked on " << pScheme->Info() << std::endl;
        mName = "ArcLengthDisplacementControlSupportScheme<" + pScheme->Info() + ">";
    }

    /** Destructor.*/
    virtual ~ArcLengthDisplacementControlSupportScheme()
    {}

    /**@name Operators */

    /*@} */
    /**@name Operations */
    /*@{ */

    const TSystemVectorType& GetForceVector() const
    {
        return mForceVector;
    }

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        BaseType::InitializeSolutionStep(r_model_part, A, Dx, b);

        mEquationSystemSize = SparseSpaceType::Size(b);
        if (SparseSpaceType::Size(mForceVector) != mEquationSystemSize)
            SparseSpaceType::Resize(mForceVector, mEquationSystemSize);
        SparseSpaceType::SetToZero(mForceVector);
    }

    void InitializeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        BaseType::InitializeNonLinIteration(r_model_part, A, Dx, b);

        SparseSpaceType::SetToZero(mForceVector);
    }

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        BaseType::CalculateSystemContributions(rCurrentElement, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

        // compute the forces induced by displacement
        LocalSystemVectorType Force;
        DenseSpaceType::Resize(Force, DenseSpaceType::Size(RHS_Contribution));
        DenseSpaceType::SetToZero(Force);
        try
        {
            const PrescribedObject& rObject = dynamic_cast<PrescribedObject&>(rCurrentElement);
            rObject.ComputePrescribedForces(LHS_Contribution, Force, CurrentProcessInfo);
        }
        catch (std::bad_cast& bc)
        {
            KRATOS_THROW_ERROR(std::logic_error, "The element does not derive from PrescribedObject", "")
        }

        // assemble the force vector
        Element::DofsVectorType ElementalDofList;
        this->GetDofList(rCurrentElement, ElementalDofList, CurrentProcessInfo);

        std::size_t cnt = 0;
        for (Element::DofsVectorType::iterator it = ElementalDofList.begin(); it != ElementalDofList.end(); ++it)
        {
            if ((*it)->GetVariable() == DISPLACEMENT_X
             || (*it)->GetVariable() == DISPLACEMENT_Y
             || (*it)->GetVariable() == DISPLACEMENT_Z )
            {
                const auto eq_id = (*it)->EquationId();
                if (eq_id < mEquationSystemSize)
                {
                    mForceVector[eq_id] += Force[cnt];
                }
            }

            ++cnt;
        }
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
        return mName;
    }

    /*@} */
    /**@name Friends */
    /*@{ */

private:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */

    std::string mName;
    SizeType mEquationSystemSize;
    TSystemVectorType mForceVector;

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

#endif /* KRATOS_STRUCTURAL_APPLICATION_ARC_LENGTH_DISPLACEMENT_CONTROL_SUPPORT_SCHEME  defined */
