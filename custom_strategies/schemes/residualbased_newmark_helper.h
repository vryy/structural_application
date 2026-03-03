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
 * Helper for Time Integration scheme based on Generalized-Alpha Method.
 *
 * In the typical dynamics system, the element equilibrium is represented by
 *  rinert(u,udd) + rvis(u,ud) + rint(u) = rext
 * where
 *  -   rinert(u,udd) = M(u)*udd
 *  -   rvist(u,ud) = D(u,ud)*ud
 *  -   u,ud,udd are time-interpolated primal variables, e.g., displacement, velocity, acceleration.
 *      And U is the primal unknown, e.g., u_{n+1}
 * The linearization of the time integration scheme follows by:
 *      LHS = d(rint)/d(u)*d(u)/d(U)
 *          + d(rvis)/(du)*d(u)/d(U) + d(rvis)/(dud)*d(ud)/d(U)
 *          + d(rinert)/(du)*d(u)/d(U) + d(rinert)/(dudd)*d(udd)/d(U)
 *          = K*d(u)/d(U)
 *          + [d(D)/(du)*ud]*d(u)/d(U) + [d(D)/d(ud)*ud + D]*d(ud)/d(U)
 *          + [d(M)/(du)*udd]*d(u)/d(U) + M*d(udd)/d(U)
 *
 * The term d(D)/(du)*ud is called damping-induced stiffness matrix and
 * d(M)/(du)*udd as mass-induced stiffness matrix. The handling of these two matrices
 * defines the type of element which deals with nonlinear mass damping.
 */
template<int TMassDampingType>
struct ResidualBasedNewmarkHelper;

/**
 * In this instance of helper (TMassDampingType=0), following assumptions are made:
 * -    CalculateLocalSystem: K
 * -    CalculateDampingMatrix: D -> unused
 * -    CalculateMassMatrix: M -> unused
 * -    CalculateLocalVelocityContribution: -> unused
 * -    CalculateLocalAccelerationContribution: -> unused
 * The element is assumed as static only and computes the stiffness matrix using the
 * time-interpolated primal variables.
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

        ///

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        LHS_Contribution *= (1 - alpha_f);
    }

    template<class TEntityType, typename LocalSystemVectorType>
    static void CalculateRHSContribution(
        TEntityType& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    template<class TEntityType, typename LocalSystemMatrixType>
    static void CalculateLHSContribution(
        TEntityType& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        const double alpha_f = CurrentProcessInfo[NEWMARK_ALPHAF];

        ///

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        LHS_Contribution *= (1 - alpha_f);
    }
}; /* struct ResidualBasedNewmarkHelper<0> */

/**
 * In this instance of helper (TMassDampingType=1), corresponding to linear mass damping case,
 * following assumptions are made:
 * -    CalculateLocalSystem: K + d(D)/(du)*ud + d(M)/(du)*udd
 * -    CalculateDampingMatrix: D
 * -    CalculateMassMatrix: M
 * -    d(D)/d(ud)*ud must be zero. Hence D only depends on u: D = D(u)
 * Since D=D(u), the element with velocity-dependent damping may not work and quadratic
 * convergent will not be obtained.
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
        const double alpha_m = CurrentProcessInfo[NEWMARK_ALPHAM];
        const double beta = CurrentProcessInfo[NEWMARK_BETA];
        const double gamma = CurrentProcessInfo[NEWMARK_GAMMA];

        ///

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

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

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

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

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

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
}; /* struct ResidualBasedNewmarkHelper<1> */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/**
 * In this instance of helper (TMassDampingType=2), following assumptions are made:
 * -    CalculateLocalSystem: K
 * -    CalculateDampingMatrix: D -> unused
 * -    CalculateMassMatrix: M -> unused
 * -    CalculateLocalVelocityContribution: [d(D)/d(ud)*ud + D] as rDampingMatrix
 *          and [d(D)/(du)*ud] as damping-induced stiffness matrix
 * -    CalculateLocalAccelerationContribution: M as rMassMatrix
 *          and [d(M)/(du)*udd] as mass-induced stiffness matrix
 * It is noted that the damping matrix in CalculateLocalVelocityContribution is different
 * to the damping matrix returned by CalculateDampingMatrix
 */
template<>
struct ResidualBasedNewmarkHelper<2>
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

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        const double aux = (1 - alpha_f);
        LHS_Contribution *= aux;

        if (CurrentProcessInfo[QUASI_STATIC_ANALYSIS])
        {
            LocalSystemMatrixType DampingMatrix;
            LocalSystemMatrixType DampingInducedStiffnessMatrix;

            rCurrentElement.CalculateLocalVelocityContribution(DampingMatrix, DampingInducedStiffnessMatrix, RHS_Contribution, CurrentProcessInfo);

            if (norm_frobenius(DampingMatrix) > 0.0)
            {
                const double gamma = CurrentProcessInfo[NEWMARK_GAMMA];
                const double beta = CurrentProcessInfo[NEWMARK_BETA];
                const double Dt = CurrentProcessInfo[DELTA_TIME];

                const double aux1 = aux * gamma/(beta*Dt);
                noalias(LHS_Contribution) += aux1 * DampingMatrix;
            }

            if (norm_frobenius(DampingInducedStiffnessMatrix) > 0.0)
                noalias(LHS_Contribution) += DampingInducedStiffnessMatrix;
        }
        else
        {
            LocalSystemMatrixType DampingMatrix;
            LocalSystemMatrixType DampingInducedStiffnessMatrix;

            rCurrentElement.CalculateLocalVelocityContribution(DampingMatrix, DampingInducedStiffnessMatrix, RHS_Contribution, CurrentProcessInfo);

            if (norm_frobenius(DampingMatrix) > 0.0)
            {
                const double gamma = CurrentProcessInfo[NEWMARK_GAMMA];
                const double beta = CurrentProcessInfo[NEWMARK_BETA];
                const double Dt = CurrentProcessInfo[DELTA_TIME];

                const double aux1 = aux * gamma/(beta*Dt);
                noalias(LHS_Contribution) += aux1 * DampingMatrix;
            }

            if (norm_frobenius(DampingInducedStiffnessMatrix) > 0.0)
                noalias(LHS_Contribution) += DampingInducedStiffnessMatrix;

            //

            LocalSystemMatrixType MassMatrix;
            LocalSystemMatrixType MassInducedStiffnessMatrix;

            rCurrentElement.CalculateLocalAccelerationContribution(MassMatrix, MassInducedStiffnessMatrix, RHS_Contribution, CurrentProcessInfo);

            if (norm_frobenius(MassMatrix) > 0.0)
            {
                const double alpha_m = CurrentProcessInfo[NEWMARK_ALPHAM];
                const double beta = CurrentProcessInfo[NEWMARK_BETA];
                const double Dt = CurrentProcessInfo[DELTA_TIME];

                const double aux1 = (1 - alpha_m) / (beta*pow(Dt, 2));
                noalias(LHS_Contribution) += aux1 * MassMatrix;
            }

            if (norm_frobenius(MassInducedStiffnessMatrix) > 0.0)
                noalias(LHS_Contribution) += MassInducedStiffnessMatrix;
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

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

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
        KRATOS_ERROR << __FUNCTION__ << " is not supported. Call CalculateSystemContributions instead.";
    }
}; /* struct ResidualBasedNewmarkHelper<2> */

/**
 * In this instance of helper (TMassDampingType=3), following assumptions are made:
 * -    CalculateLocalSystem: K
 * -    CalculateDampingMatrix: D -> unused
 * -    CalculateMassMatrix: M -> unused
 * -    CalculateLocalVelocityContribution: [d(D)/(du)*ud]*d(u)/d(U) + [d(D)/d(ud)*ud + D]*d(ud)/d(U)
 *              as the left hand side term
 * -    CalculateLocalAccelerationContribution: [d(M)/(du)*udd]*d(u)/d(U) + M*d(udd)/d(U)
                as the left hand side term
 * The element has to deal with time integration, hence it has to access time integration parameters.
 * An example of this type of element is AD element, then one is refrained to perform the linearization
 * to obtain the damping- and mass-induced stiffness matrix.
 */
template<>
struct ResidualBasedNewmarkHelper<3>
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

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

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

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

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
        const double alpha_f = CurrentProcessInfo[NEWMARK_ALPHAF];

        ///

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        LHS_Contribution *= (1 - alpha_f);
    }
}; /* struct ResidualBasedNewmarkHelper<3> */

/**
 * In this instance of helper (TMassDampingType=3), following assumptions are made:
 * -    CalculateLocalSystem: everything
 * -    CalculateDampingMatrix: D -> unused
 * -    CalculateMassMatrix: M -> unused
 * -    CalculateLocalVelocityContribution: -> unused
 * -    CalculateLocalAccelerationContribution: -> unused
 * The element has to deal with time integration all at once in CalculateLocalSystem, hence it has to access time integration parameters.
 * An example of this type of element is AD element, then one is refrained to perform the linearization to obtain any kind of matrices.
 */
template<>
struct ResidualBasedNewmarkHelper<4>
{
    template<class TEntityType, typename LocalSystemMatrixType, typename LocalSystemVectorType>
    static void CalculateSystemContributions(
        TEntityType& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    template<class TEntityType, typename LocalSystemVectorType>
    static void CalculateRHSContribution(
        TEntityType& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    template<class TEntityType, typename LocalSystemMatrixType>
    static void CalculateLHSContribution(
        TEntityType& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        typename TEntityType::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);
        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
    }
}; /* struct ResidualBasedNewmarkHelper<4> */

}  /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APPLICATION_RESIDUALBASED_NEWMARK_HELPER  defined */
