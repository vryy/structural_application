/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dxdvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.eDx
rrossi@cimne.upc.eDx
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
CLAIM, DxMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */
/* *********************************************************
 *
 *   Last Modified by:    $Author: janosch $
 *   Dxte:                $Dxte: 2007-04-13 15:59:32 $
 *   Revision:            $Revision: 1.2 $
 *
 * ***********************************************************/


#if !defined(KRATOS_MULTIPHASEFLOW_CRITERIA )
#define  KRATOS_MULTIPHASEFLOW_CRITERIA


/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"
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

/** Short class definition.
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
         class TDenseSpace
         >
class MultiPhaseFlowCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( MultiPhaseFlowCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    MultiPhaseFlowCriteria(
        TDataType RelativeTolerance,
        TDataType AbsoluteTolerance)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mRelativeTolerance = RelativeTolerance;
        mAbsoluteTolerance = AbsoluteTolerance;

        // mCheckType = 1; // only check the absolute criteria
        // mCheckType = 2; // only check the relative criteria
        // mCheckType = 3; // check both
        mCheckType = 4; // one of them
    }

    /** Destructor.
     */
    virtual ~MultiPhaseFlowCriteria()
    {
    }


    void SetType(const int& Type)
    {
        mCheckType = Type;
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    /*Criterias that need to be called before getting the solution */
    virtual bool PreCriteria(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        return true;
    }

    /*Criterias that need to be called after getting the solution */
    bool PostCriteria(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        if (Dx.size() != 0) //if we are solving for something
        {
            double norm_Dx = 0.0;
            double norm_b = 0.0;
            double norm_x = 0.0;
            double norm_Dx_WATER = 0.0;
            double norm_b_WATER = 0.0;
            double norm_x_WATER = 0.0;
            double norm_Dx_AIR = 0.0;
            double norm_b_AIR = 0.0;
            double norm_x_AIR = 0.0;

            bool HasWaterPres = false;
            bool HasAirPres = false;

            for (typename DofsArrayType::iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
            {
                if (i_dof->IsFree())
                {
                    if (i_dof->GetVariable() == DISPLACEMENT_X)
                    {
                        norm_Dx += Dx[i_dof->EquationId()] * Dx[i_dof->EquationId()];
                        norm_b += b[i_dof->EquationId()] * b[i_dof->EquationId()];
                        norm_x += i_dof->GetSolutionStepValue(DISPLACEMENT_X) * i_dof->GetSolutionStepValue(DISPLACEMENT_X);
                    }
                    if (i_dof->GetVariable() == DISPLACEMENT_Y)
                    {
                        norm_Dx += Dx[i_dof->EquationId()] * Dx[i_dof->EquationId()];
                        norm_b += b[i_dof->EquationId()] * b[i_dof->EquationId()];
                        norm_x += i_dof->GetSolutionStepValue(DISPLACEMENT_Y) * i_dof->GetSolutionStepValue(DISPLACEMENT_Y);
                    }
                    if (i_dof->GetVariable() == DISPLACEMENT_Z)
                    {
                        norm_Dx += Dx[i_dof->EquationId()] * Dx[i_dof->EquationId()];
                        norm_b += b[i_dof->EquationId()] * b[i_dof->EquationId()];
                        norm_x += i_dof->GetSolutionStepValue(DISPLACEMENT_Z) * i_dof->GetSolutionStepValue(DISPLACEMENT_Z);
                    }
                    if (i_dof->GetVariable() == WATER_PRESSURE)
                    {
                        HasWaterPres = true;

                        norm_Dx_WATER += Dx[i_dof->EquationId()] * Dx[i_dof->EquationId()];
                        norm_b_WATER += b[i_dof->EquationId()] * b[i_dof->EquationId()];
                        norm_x_WATER += i_dof->GetSolutionStepValue(WATER_PRESSURE) * i_dof->GetSolutionStepValue(WATER_PRESSURE);
                    }
                    if (i_dof->GetVariable() == AIR_PRESSURE)
                    {
                        HasAirPres = true;

                        norm_Dx_AIR += Dx[i_dof->EquationId()] * Dx[i_dof->EquationId()];
                        norm_b_AIR += b[i_dof->EquationId()] * b[i_dof->EquationId()];
                        norm_x_AIR += i_dof->GetSolutionStepValue(AIR_PRESSURE) * i_dof->GetSolutionStepValue(AIR_PRESSURE);
                    }
                }
            }

            norm_x = sqrt(norm_x);
            norm_Dx = sqrt(norm_Dx);
            norm_b = sqrt(norm_b);

            if (HasWaterPres)
            {
                norm_x_WATER = sqrt(norm_x_WATER);
                norm_Dx_WATER = sqrt(norm_Dx_WATER);
                norm_b_WATER = sqrt(norm_b_WATER);
            }

            if (HasAirPres)
            {
                norm_x_AIR = sqrt(norm_x_AIR);
                norm_Dx_AIR = sqrt(norm_Dx_AIR);
                norm_b_AIR = sqrt(norm_b_AIR);
            }

            double ratioDisp = 1.0;

            double ratioWater = 1.0;

            double ratioAir = 1.0;

            if (norm_x > 0)
                ratioDisp = norm_Dx / norm_x;
            if (norm_Dx == 0.0)
                ratioDisp = 0.0;

            if (norm_x_WATER > 0)
                ratioWater = norm_Dx_WATER / norm_x_WATER;
            if (norm_Dx_WATER == 0.0)
                ratioWater = 0.0;

            if (norm_x_AIR > 0)
                ratioAir = norm_Dx_AIR / norm_x_AIR;
            if (norm_Dx_AIR == 0.0)
                ratioAir = 0.0;

            std::cout << "********************************************CONVERGENCE CRITERIA FOR MULTIPHASE PROBLEMS********************************************" << std::endl;
            std::cout.precision(6);
            std::cout.setf(std::ios::scientific);
            std::cout << "** expected values: \t\t\t\t\t\tabs_tol = " << mAbsoluteTolerance << "\t\t\t\trel_tol = " << mRelativeTolerance << " **" << std::endl;
            std::cout << "** obtained values displacement:\tratio = " << ratioDisp << "\t||Dx|| = " << norm_Dx << "\t||x|| = " << norm_x << "\t||b|| = " << norm_b << " **" << std::endl;
            if (HasWaterPres)
            {
                std::cout << "** obtained values water pressure:\tratio = " << ratioWater << "\t||Dx|| = " << norm_Dx_WATER << "\t||x|| = " << norm_x_WATER << " \t||b|| = " << norm_b_WATER << " **" << std::endl;
                if (HasAirPres)
                {
                    std::cout << "** obtained values air pressure:\tratio = " << ratioAir << "\t||Dx|| = " << norm_Dx_AIR << "\t||x|| = " << norm_x_AIR << "\t||b|| = " << norm_b_AIR << " **" << std::endl;

                    std::cout << "** obtained values total:\t\tratio = " << ratioAir + ratioWater + ratioDisp << "\tchange = " << norm_Dx + norm_Dx_WATER + norm_Dx_AIR << "\tabsolute = " << norm_x_AIR + norm_x_WATER + norm_x << "\tenergy = " << norm_b_AIR + norm_b_WATER + norm_b << " **" << std::endl;
                }
                else
                {
                    std::cout << "** obtained values total:\t\tratio = " << ratioWater + ratioDisp << "\t change = " << norm_Dx + norm_Dx_WATER << "\tabsolute = " << norm_x_WATER + norm_x << "\tenergy = " << norm_b_WATER + norm_b << " **" << std::endl;
                }
            }
            std::cout << "************************************************************************************************************************************" << std::endl;

            bool disp_reason_1 = (norm_b <= mAbsoluteTolerance);
            bool disp_reason_2 = false;
            int disp_reason_2_case = 0;
            if (ratioDisp <= mRelativeTolerance)
            {
                disp_reason_2 = true;
                disp_reason_2_case = 1;
            }
            else
            {
                if (norm_x < mAbsoluteTolerance) // case that zero state is solved
                {
                    if (norm_Dx < mAbsoluteTolerance)
                    {
                        disp_reason_2 = true;
                        disp_reason_2_case = 2;
                    }
                }
            }
            bool disp_converged;

            if (mCheckType == 1)
                disp_converged = disp_reason_1;
            else if (mCheckType == 2)
                disp_converged = disp_reason_2;
            else if (mCheckType == 3)
                disp_converged = disp_reason_1 && disp_reason_2;
            else if (mCheckType == 4)
                disp_converged = disp_reason_1 || disp_reason_2;

            bool water_converged;
            bool water_reason_1;
            bool water_reason_2;
            if(HasWaterPres)
            {
                water_reason_1 = (norm_b_WATER <= mAbsoluteTolerance);
                water_reason_2 = (ratioWater <= mRelativeTolerance);
                if (mCheckType == 1)
                    water_converged = water_reason_1;
                else if (mCheckType == 2)
                    water_converged = water_reason_2;
                else if (mCheckType == 3)
                    water_converged = water_reason_1 && water_reason_2;
                else if (mCheckType == 4)
                    water_converged = water_reason_1 || water_reason_2;
            }
            else
                water_converged = true;

            bool air_converged;
            bool air_reason_1;
            bool air_reason_2;
            if(HasAirPres)
            {
                air_reason_1 = (norm_b_AIR <= mAbsoluteTolerance);
                air_reason_2 = (ratioAir <= mRelativeTolerance);
                if (mCheckType == 1)
                    air_converged = air_reason_1;
                else if (mCheckType == 2)
                    air_converged = air_reason_2;
                else if (mCheckType == 3)
                    air_converged = air_reason_1 && air_reason_2;
                else if (mCheckType == 4)
                    air_converged = air_reason_1 || air_reason_2;
            }
            else
                air_converged = true;

            if(disp_converged && water_converged && air_converged)
            {
                std::cout << "Congratulations the solution strategy is converged." << std::endl;
                std::cout << "Reason for converged displacement:";
                if (mCheckType == 1)
                {
                    std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                }
                else if (mCheckType == 2)
                {
                    if (disp_reason_2_case == 1)
                        std::cout << " {(||Dx||/||x|| = " << ratioDisp << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    else if (disp_reason_2_case == 2)
                        std::cout << " {(||Dx|| = " << norm_Dx << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                  << " and {(||x|| = " << norm_x << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                  << std::endl;
                }
                else if (mCheckType == 3)
                {
                    std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                    if (disp_reason_2 == 1)
                        std::cout << "and {(||Dx||/||x|| = " << ratioDisp << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    else if (disp_reason_2_case == 2)
                        std::cout << "and {(||Dx|| = " << norm_Dx << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                  << " and {(||x|| = " << norm_x << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                  << std::endl;
                }
                else if (mCheckType == 4)
                {
                    if(disp_reason_1) std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                    if (disp_reason_2)
                    {
                        if (disp_reason_2_case == 1)
                            std::cout << " {(||Dx||/||x|| = " << ratioDisp << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                        else if (disp_reason_2_case == 2)
                            std::cout << " {(||Dx|| = " << norm_Dx << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                      << " and {(||x|| = " << norm_x << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                      << std::endl;
                    }
                }

                if(HasWaterPres)
                {
                    std::cout << "Reason for converged water pressure:";
                    if (mCheckType == 1)
                    {
                        std::cout << " {(||b|| (water) = " << norm_b_WATER << ") <= (expected ||b|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 2)
                    {
                        std::cout << " {(||Dx||/||x|| (water) = " << ratioWater << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 3)
                    {
                        std::cout << " {(||b|| (water) = " << norm_b_WATER << ") <= (expected ||b|| = " << mRelativeTolerance << ")}" << std::endl;
                        std::cout << "and {(||Dx||/||x|| (water) = " << ratioWater << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 4)
                    {
                        if(water_reason_1) std::cout << " {(||b|| (water) = " << norm_b_WATER << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                        if(water_reason_2) std::cout << " {(||Dx||/||x|| (water) = " << ratioWater << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                }

                if(HasAirPres)
                {
                    std::cout << "Reason for converged air pressure:";
                    if (mCheckType == 4)
                    {
                        std::cout << " {(||b|| (air) = " << norm_b_AIR << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 2)
                    {
                        std::cout << " {(||Dx||/||x|| (air) = " << ratioAir << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 3)
                    {
                        std::cout << " {(||b|| (air) = " << norm_b_AIR << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                        std::cout << "and {(||Dx||/||x|| (air) = " << ratioAir << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 4)
                    {
                        if(air_reason_1) std::cout << " {(||b|| (air) = " << norm_b_AIR << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                        if(air_reason_2) std::cout << " {(||Dx||/||x|| (air) = " << ratioAir << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                }

                return true;
            }
            else
                return false;
        }
        else   //in this case all the displacements are imposed!
        {
            return true;
        }

    }

    void Initialize(
        ModelPart& r_model_part
    )
    {
    }

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
    }

    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
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

    TDataType mRelativeTolerance;
    TDataType mAbsoluteTolerance;
    int mCheckType;

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

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_MULTIPHASEFLOW_CRITERIA  defined */
