//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 13/3/2022 $
//
//


#if !defined(KRATOS_ARC_LENGTH_ENERGY_RELEASE_CONTROL_PROCESS_H_INCLUDED )
#define  KRATOS_ARC_LENGTH_ENERGY_RELEASE_CONTROL_PROCESS_H_INCLUDED


#include "custom_processes/arc_length_control_process.h"


namespace Kratos
{


/**
 * Arc-length energy release control
 * In this constraint, all the free d.o.fs are accounted
 * Reference:
 * +    Gutierrez et al, Energy release control for numerical simulations of failure in quasi-brittle solids
 */
template<class TBuilderAndSolverType>
class ArcLengthEnergyReleaseControlProcess : public ArcLengthControlProcess<TBuilderAndSolverType>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthEnergyReleaseControlProcess );

    typedef ArcLengthControlProcess<TBuilderAndSolverType> BaseType;

    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    ArcLengthEnergyReleaseControlProcess(const double& Radius)
    : BaseType(), mRadius(Radius)
    {
        std::cout << "ArcLengthEnergyReleaseControlProcess is used, radius = " << mRadius << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    void SetRadius(const double& Radius)
    {
        mRadius = Radius;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Get the value of the constraint
    double GetValue() const override
    {
        const auto EquationSystemSize = this->GetBuilderAndSolver().GetEquationSystemSize();
        const DofsArrayType& rDofSet = this->GetBuilderAndSolver().GetDofSet();

        double f = 0.0;
        const TSystemVectorType& fhat = this->GetForceVector();

        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if ( dof_iterator->GetVariable() == DISPLACEMENT_X
              || dof_iterator->GetVariable() == DISPLACEMENT_Y
              || dof_iterator->GetVariable() == DISPLACEMENT_Z )
            {
                const auto row = dof_iterator->EquationId();
                if (row < EquationSystemSize)
                {
                    double fv = TSparseSpaceType::GetValue(fhat, row);
                    f += 0.5 * ( this->LambdaOld() * (dof_iterator->GetSolutionStepValue() - dof_iterator->GetSolutionStepValue(1))
                               - (this->Lambda() - this->LambdaOld()) * dof_iterator->GetSolutionStepValue(1) ) * fv;
                }
            }
        }

        return f - mRadius;
    }

    /// Get the derivatives of the constraint
    TSystemVectorType GetDerivativesDU() const override
    {
        const auto EquationSystemSize = this->GetBuilderAndSolver().GetEquationSystemSize();
        const DofsArrayType& rDofSet = this->GetBuilderAndSolver().GetDofSet();

        TSystemVectorType dfdu;
        TSparseSpaceType::Resize(dfdu, EquationSystemSize);
        TSparseSpaceType::SetToZero(dfdu);
        const TSystemVectorType& fhat = this->GetForceVector();

        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if ( dof_iterator->GetVariable() == DISPLACEMENT_X
              || dof_iterator->GetVariable() == DISPLACEMENT_Y
              || dof_iterator->GetVariable() == DISPLACEMENT_Z )
            {
                const auto row = dof_iterator->EquationId();
                if (row < EquationSystemSize)
                {
                    double fv = TSparseSpaceType::GetValue(fhat, row);
                    TSparseSpaceType::SetValue(dfdu, row, 0.5 * this->LambdaOld() * fv);
                }
            }
        }

        return dfdu;
    }

    /// Get the derivatives of the constraint
    double GetDerivativesDLambda() const override
    {
        const auto EquationSystemSize = this->GetBuilderAndSolver().GetEquationSystemSize();
        const DofsArrayType& rDofSet = this->GetBuilderAndSolver().GetDofSet();

        double df = 0.0;
        const TSystemVectorType& fhat = this->GetForceVector();

        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if ( dof_iterator->GetVariable() == DISPLACEMENT_X
              || dof_iterator->GetVariable() == DISPLACEMENT_Y
              || dof_iterator->GetVariable() == DISPLACEMENT_Z )
            {
                const auto row = dof_iterator->EquationId();
                if (row < EquationSystemSize)
                {
                    double fv = TSparseSpaceType::GetValue(fhat, row);
                    df -= 0.5 * dof_iterator->GetSolutionStepValue(1) * fv;
                }
            }
        }

        return df;
    }

    double Predict(const TSystemVectorType& rDeltaUl) const override
    {
        const auto EquationSystemSize = this->GetBuilderAndSolver().GetEquationSystemSize();
        const DofsArrayType& rDofSet = this->GetBuilderAndSolver().GetDofSet();

        /* compute the forward criteria */
        // assemble delta u and s0
        double norm_p2_dul = 0.0;
        TSystemVectorType Du;
        TSparseSpaceType::Resize(Du, EquationSystemSize);
        TSparseSpaceType::SetToZero(Du);
        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if ( dof_iterator->GetVariable() == DISPLACEMENT_X
              || dof_iterator->GetVariable() == DISPLACEMENT_Y
              || dof_iterator->GetVariable() == DISPLACEMENT_Z )
            {
                const auto row = dof_iterator->EquationId();
                // KRATOS_WATCH(dof_iterator->GetVariable().Name())
                if (row < EquationSystemSize)
                {
                    TSparseSpaceType::SetValue(Du, row, dof_iterator->GetSolutionStepValue(1) - dof_iterator->GetSolutionStepValue(2));

                    double tmp = TSparseSpaceType::GetValue(rDeltaUl, row);
                    norm_p2_dul += tmp*tmp;
                }
            }
        }

        double s0 = sqrt(norm_p2_dul);
        // KRATOS_WATCH(__LINE__)
        // KRATOS_WATCH(s0)

        if (this->IsForcedReverse() || this->IsForcedForward())
        {
            if (this->IsForcedReverse())
            {
                s0 *= -1.0;
                std::cout << "ArcLengthEnergyReleaseControlProcess: forward sign is forced to be reversed" << std::endl;
            }

            if (this->IsForcedForward())
            {
                std::cout << "ArcLengthEnergyReleaseControlProcess: forward sign is forced to remain" << std::endl;
            }
        }
        else
        {
            // compute the forward criteria
            double forward_criteria = TSparseSpaceType::Dot(Du, rDeltaUl);
            // KRATOS_WATCH(Du)
            // KRATOS_WATCH(*mp_delta_u_l)
            // KRATOS_WATCH(TSparseSpaceType::Dot(Du, *mp_delta_u_l))
            // KRATOS_WATCH(forward_criteria)

            if (forward_criteria < -1.0e-10)
            {
                s0 *= -1.0;
                std::cout << "ArcLengthEnergyReleaseControlProcess: forward sign is reversed" << std::endl;
            }
            else
            {
                std::cout << "ArcLengthEnergyReleaseControlProcess: forward sign remains" << std::endl;
            }
        }

        // compute delta lambda
        double delta_lambda = 1.0/s0 * mRadius;
        return delta_lambda;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ArcLengthEnergyReleaseControlProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ArcLengthEnergyReleaseControlProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Radius: " << mRadius << std::endl;
        BaseType::PrintData(rOStream);
    }

private:

    double mRadius;

}; // Class ArcLengthEnergyReleaseControlProcess

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_ENERGY_RELEASE_CONTROL_PROCESS_H_INCLUDED defined
