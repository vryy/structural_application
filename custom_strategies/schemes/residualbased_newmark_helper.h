/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 23/07/2023 $
*   Revision:            $Revision: 1.9 $
*
* ***********************************************************/


#if !defined(KRATOS_STRUCTURAL_APPLICATION_RESIDUALBASED_NEWMARK_HELPER )
#define  KRATOS_STRUCTURAL_APPLICATION_RESIDUALBASED_NEWMARK_HELPER

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "structural_application_variables.h"

namespace Kratos
{

/**
 * Helper for Time Integration scheme based on Generalized-Alpha Method
 * TMassDampingType = 0: inertial forces and damping forces are computed linearly, i.e.
 *          finert = Mass * acceleration
 *      and
 *          fdamp = Damp * velocity
 *
 * TMassDampingType = 1: inertial forces and damping forces are assumed nonlinear, the
 * contribution has to be taken into account at the element level
 */
template<int TMassDampingType>
struct ResidualBasedNewmarkHelper;

/**
 * Local assembly assuming linear inertial and damping contribution
 */
template<>
struct ResidualBasedNewmarkHelper<0>
{

    template<class TEntityType, typename LocalSystemMatrixType, typename LocalSystemVectorType>
    static void CalculateSystemContributions(
        TEntityType& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        const double alpha_f = CurrentProcessInfo[NEWMARK_ALPHAF];
        const double alpha_m = CurrentProcessInfo[NEWMARK_ALPHAM];
        const double beta = CurrentProcessInfo[NEWMARK_BETA];
        const double gamma = CurrentProcessInfo[NEWMARK_GAMMA];

        ///

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        if (CurrentProcessInfo[QUASI_STATIC_ANALYSIS])
        {
            Matrix DampingMatrix;

            rCurrentElement.CalculateDampingMatrix(DampingMatrix, CurrentProcessInfo);

            if (norm_frobenius(DampingMatrix) > 0.0)
            {
                AddDampingToRHS(rCurrentElement, RHS_Contribution, DampingMatrix, CurrentProcessInfo);
                AssembleTimeSpaceLHS_QuasiStatic(LHS_Contribution, DampingMatrix, CurrentProcessInfo, alpha_f, beta, gamma);
            }
            else
                // the element that has no damping matrix shall not contribute the viscosity part to the system
                AssembleTimeSpaceLHS_QuasiStatic_NoDamping(LHS_Contribution, CurrentProcessInfo, alpha_f);
        }
        else
        {
            Matrix DampingMatrix;

            Matrix MassMatrix;

            rCurrentElement.CalculateDampingMatrix(DampingMatrix, CurrentProcessInfo);

            rCurrentElement.CalculateMassMatrix(MassMatrix, CurrentProcessInfo);

            if ((norm_frobenius(DampingMatrix) == 0.0) && (norm_frobenius(MassMatrix) > 0.0))
            {
                AddInertiaToRHS(rCurrentElement, RHS_Contribution, MassMatrix, CurrentProcessInfo);
                AssembleTimeSpaceLHS_Dynamics_NoDamping(LHS_Contribution, MassMatrix, CurrentProcessInfo, alpha_f, alpha_m, beta);
            }
            else if ((norm_frobenius(DampingMatrix) > 0.0) && (norm_frobenius(MassMatrix) == 0.0))
            {
                AddDampingToRHS(rCurrentElement, RHS_Contribution, DampingMatrix, CurrentProcessInfo);
                AssembleTimeSpaceLHS_QuasiStatic(LHS_Contribution, DampingMatrix, CurrentProcessInfo, alpha_f, beta, gamma);
            }
            else if ((norm_frobenius(DampingMatrix) > 0.0) && (norm_frobenius(MassMatrix) > 0.0))
            {
                AddInertiaToRHS(rCurrentElement, RHS_Contribution, MassMatrix, CurrentProcessInfo);
                AddDampingToRHS(rCurrentElement, RHS_Contribution, DampingMatrix, CurrentProcessInfo);
                AssembleTimeSpaceLHS_Dynamics(LHS_Contribution, DampingMatrix, MassMatrix, CurrentProcessInfo, alpha_f, alpha_m, beta, gamma);
            }
            else
            {
                AssembleTimeSpaceLHS_Dynamics_NoMass_NoDamping(LHS_Contribution, CurrentProcessInfo, alpha_f);
            }
        }
    }

    template<class TEntityType, typename LocalSystemVectorType>
    static void CalculateRHSContribution(
        TEntityType& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        if (CurrentProcessInfo[QUASI_STATIC_ANALYSIS])
        {
            Matrix DampingMatrix;

            rCurrentElement.CalculateDampingMatrix(DampingMatrix, CurrentProcessInfo);

            if (norm_frobenius(DampingMatrix) > 0.0)
            {
                AddDampingToRHS(rCurrentElement, RHS_Contribution, DampingMatrix, CurrentProcessInfo);
            }
        }
        else
        {
            Matrix DampingMatrix;

            Matrix MassMatrix;

            rCurrentElement.CalculateDampingMatrix(DampingMatrix, CurrentProcessInfo);

            rCurrentElement.CalculateMassMatrix(MassMatrix, CurrentProcessInfo);

            if ((norm_frobenius(DampingMatrix) == 0.0) && (norm_frobenius(MassMatrix) > 0.0))
            {
                AddInertiaToRHS(rCurrentElement, RHS_Contribution, MassMatrix, CurrentProcessInfo);
            }
            else if ((norm_frobenius(DampingMatrix) > 0.0) && (norm_frobenius(MassMatrix) == 0.0))
            {
                AddDampingToRHS(rCurrentElement, RHS_Contribution, DampingMatrix, CurrentProcessInfo);
            }
            else if ((norm_frobenius(DampingMatrix) > 0.0) && (norm_frobenius(MassMatrix) > 0.0))
            {
                AddInertiaToRHS(rCurrentElement, RHS_Contribution, MassMatrix, CurrentProcessInfo);
                AddDampingToRHS(rCurrentElement, RHS_Contribution, DampingMatrix, CurrentProcessInfo);
            }
            else
            {
                // DO NOTHING
            }
        }
    }

    template<class TEntityType, typename LocalSystemMatrixType>
    static void CalculateLHSContribution(
        TEntityType& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        const double alpha_f = CurrentProcessInfo[NEWMARK_ALPHAF];
        const double alpha_m = CurrentProcessInfo[NEWMARK_ALPHAM];
        const double beta = CurrentProcessInfo[NEWMARK_BETA];
        const double gamma = CurrentProcessInfo[NEWMARK_GAMMA];

        ///

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        if (CurrentProcessInfo[QUASI_STATIC_ANALYSIS])
        {
            Matrix DampingMatrix;

            rCurrentElement.CalculateDampingMatrix(DampingMatrix, CurrentProcessInfo);

            if (norm_frobenius(DampingMatrix) > 0.0)
                AssembleTimeSpaceLHS_QuasiStatic(LHS_Contribution, DampingMatrix, CurrentProcessInfo, alpha_f, beta, gamma);
            else
                // the element that has no damping matrix shall not contribute the viscosity part to the system
                AssembleTimeSpaceLHS_QuasiStatic_NoDamping(LHS_Contribution, CurrentProcessInfo, alpha_f);
        }
        else
        {
            Matrix DampingMatrix;

            Matrix MassMatrix;

            rCurrentElement.CalculateDampingMatrix(DampingMatrix, CurrentProcessInfo);

            rCurrentElement.CalculateMassMatrix(MassMatrix, CurrentProcessInfo);

            if ((norm_frobenius(DampingMatrix) == 0.0) && (norm_frobenius(MassMatrix) > 0.0))
            {
                AssembleTimeSpaceLHS_Dynamics_NoDamping(LHS_Contribution, MassMatrix, CurrentProcessInfo, alpha_f, alpha_m, beta);
            }
            else if ((norm_frobenius(DampingMatrix) > 0.0) && (norm_frobenius(MassMatrix) == 0.0))
            {
                AssembleTimeSpaceLHS_QuasiStatic(LHS_Contribution, DampingMatrix, CurrentProcessInfo, alpha_f, beta, gamma);
            }
            else if ((norm_frobenius(DampingMatrix) > 0.0) && (norm_frobenius(MassMatrix) > 0.0))
            {
                AssembleTimeSpaceLHS_Dynamics(LHS_Contribution, DampingMatrix, MassMatrix, CurrentProcessInfo, alpha_f, alpha_m, beta, gamma);
            }
            else
            {
                AssembleTimeSpaceLHS_Dynamics_NoMass_NoDamping(LHS_Contribution, CurrentProcessInfo, alpha_f);
            }
        }
    }

    template<class TEntityType, typename LocalSystemMatrixType, typename LocalSystemVectorType>
    static void AddInertiaToRHS(
        const TEntityType& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        const LocalSystemMatrixType& MassMatrix,
        const ProcessInfo& CurrentProcessInfo)
    {
        // adding inertia contribution
        Vector acceleration;
        rCurrentElement.GetSecondDerivativesVector(acceleration, 0);
        noalias(RHS_Contribution) -= prod(MassMatrix, acceleration);
    }

    template<class TEntityType, typename LocalSystemMatrixType, typename LocalSystemVectorType>
    static void AddDampingToRHS(
        const TEntityType& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        const LocalSystemMatrixType& DampingMatrix,
        const ProcessInfo& CurrentProcessInfo)
    {
        // adding damping contribution
        Vector velocity;
        rCurrentElement.GetFirstDerivativesVector(velocity, 0);
        noalias(RHS_Contribution) -= prod(DampingMatrix, velocity);
    }

    template<typename LocalSystemMatrixType>
    static void AssembleTimeSpaceLHS_QuasiStatic(
        LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemMatrixType& DampingMatrix,
        const ProcessInfo& CurrentProcessInfo,
        const double alpha_f, const double beta, const double gamma)
    {
        const double Dt = CurrentProcessInfo[DELTA_TIME];
        double aux;

        // adding stiffness contribution to the dynamic stiffness
        aux = (1 - alpha_f);
        LHS_Contribution *= aux;

        // adding damping contribution to the dynamic stiffness
        aux = (1 - alpha_f)*gamma/(beta*Dt);
        noalias(LHS_Contribution) += aux * DampingMatrix;
    }

    template<typename LocalSystemMatrixType>
    static void AssembleTimeSpaceLHS_QuasiStatic_NoDamping(
        LocalSystemMatrixType& LHS_Contribution,
        const ProcessInfo& CurrentProcessInfo,
        const double alpha_f)
    {
        const double Dt = CurrentProcessInfo[DELTA_TIME];
        double aux;

        // adding stiffness contribution to the dynamic stiffness
        aux = (1 - alpha_f);
        LHS_Contribution *= aux;
    }

    template<typename LocalSystemMatrixType>
    static void AssembleTimeSpaceLHS_Dynamics(
        LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemMatrixType& DampingMatrix,
        const LocalSystemMatrixType& MassMatrix,
        const ProcessInfo& CurrentProcessInfo,
        const double alpha_f, const double alpha_m, const double beta, const double gamma)
    {
        const double Dt = CurrentProcessInfo[DELTA_TIME];
        double aux;

        // adding stiffness contribution to the dynamic stiffness
        aux = (1 - alpha_f);
        LHS_Contribution *= aux;

        // adding damping contribution to the dynamic stiffness
        aux = (1 - alpha_f) * gamma/(beta*Dt);
        noalias(LHS_Contribution) += aux * DampingMatrix;

        // adding mass contribution to the dynamic stiffness
        aux = (1 - alpha_m) / (beta*pow(Dt, 2));
        noalias(LHS_Contribution) += aux * MassMatrix;
    }

    template<typename LocalSystemMatrixType>
    static void AssembleTimeSpaceLHS_Dynamics_NoDamping(
        LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemMatrixType& MassMatrix,
        const ProcessInfo& CurrentProcessInfo,
        const double alpha_f, const double alpha_m, const double beta)
    {
        const double Dt = CurrentProcessInfo[DELTA_TIME];
        double aux;

        // adding stiffness contribution to the dynamic stiffness
        aux = (1 - alpha_f);
        LHS_Contribution *= aux;

        // adding mass contribution to the dynamic stiffness
        aux = (1 - alpha_m) / (beta*pow(Dt, 2));
        noalias(LHS_Contribution) += aux * MassMatrix;
    }

    template<typename LocalSystemMatrixType>
    static void AssembleTimeSpaceLHS_Dynamics_NoMass_NoDamping(
        LocalSystemMatrixType& LHS_Contribution,
        const ProcessInfo& CurrentProcessInfo,
        const double alpha_f)
    {
        // adding stiffness contribution to the dynamic stiffness
        LHS_Contribution *= (1 - alpha_f);
    }

}; /* struct ResidualBasedNewmarkHelper<0> */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/**
 * Local assembly assuming nonlinear inertial and damping contribution.
 * It is noted that, the element and condition shall add the contribution
 * to the inertia and damping terms with respect to the time integration scheme
 * manually.
 */
template<>
struct ResidualBasedNewmarkHelper<1>
{

    template<class TEntityType, typename LocalSystemMatrixType, typename LocalSystemVectorType>
    static void CalculateSystemContributions(
        TEntityType& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        const double alpha_f = CurrentProcessInfo[NEWMARK_ALPHAF];

        ///

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        LHS_Contribution *= (1 - alpha_f);

        if (CurrentProcessInfo[QUASI_STATIC_ANALYSIS])
        {
            rCurrentElement.CalculateLocalVelocityContribution(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        }
        else
        {
            rCurrentElement.CalculateLocalVelocityContribution(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

            rCurrentElement.CalculateLocalAccelerationContribution(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        }
    }

    template<class TEntityType, typename LocalSystemVectorType>
    static void CalculateRHSContribution(
        TEntityType& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        if (CurrentProcessInfo[QUASI_STATIC_ANALYSIS])
        {
            rCurrentElement.AddDampingForces(RHS_Contribution, 1.0, CurrentProcessInfo);
        }
        else
        {
            rCurrentElement.AddDampingForces(RHS_Contribution, 1.0, CurrentProcessInfo);

            rCurrentElement.AddInertiaForces(RHS_Contribution, 1.0, CurrentProcessInfo);
        }
    }

    template<class TEntityType, typename LocalSystemMatrixType>
    static void CalculateLHSContribution(
        TEntityType& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_ERROR << __FUNCTION__ << " is not implemented. Call CalculateSystemContributions instead.";
    }

}; /* struct ResidualBasedNewmarkHelper<1> */

}  /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APPLICATION_RESIDUALBASED_NEWMARK_HELPER  defined */
