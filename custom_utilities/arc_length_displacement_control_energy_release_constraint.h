//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 23/3/2022 $
//
//


#if !defined(KRATOS_ARC_LENGTH_DISPLACEMENT_CONTROL_ENERGY_RELEASE_CONSTRAINT_H_INCLUDED )
#define  KRATOS_ARC_LENGTH_DISPLACEMENT_CONTROL_ENERGY_RELEASE_CONSTRAINT_H_INCLUDED


#include "custom_processes/arc_length_control_process.h"


namespace Kratos
{


/**
 * Arc-length energy release control
 * In this constraint, all the free d.o.fs are accounted
 * Reference:
 *  +   Verhoosel et al, A dissipation-based arc-length method for robust simulation of brittle and ductile failure, 2009
 */
template<class TBuilderAndSolverType>
class ArcLengthDisplacementControlEnergyReleaseConstraint : public ArcLengthConstraint<TBuilderAndSolverType>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthDisplacementControlEnergyReleaseConstraint );

    typedef ArcLengthConstraint<TBuilderAndSolverType> BaseType;

    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    ArcLengthDisplacementControlEnergyReleaseConstraint(const double Radius)
    : BaseType(Radius), mp_ext_f(NULL)
    {
        std::cout << "ArcLengthDisplacementControlEnergyReleaseConstraint is used, radius = " << BaseType::Radius() << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    /// Check if the force vector is required
    int NeedForceVector() const override
    {
        return 2;
    }

    /// Set the external force vector
    void SetForceVector(const TSystemVectorType& rValue) override
    {
        mp_ext_f = &rValue;
    }

    /// Get the external force vector
    const TSystemVectorType& GetForceVector() const override
    {
        if (mp_ext_f == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Force vector is not set", "")
        return *mp_ext_f;
    }

    ///@}
    ///@name Operators
    ///@{

    void Copy(const BaseType& rOther) override
    {
        BaseType::Copy(rOther);
        try
        {
            const ArcLengthDisplacementControlEnergyReleaseConstraint& thisValue = dynamic_cast<const ArcLengthDisplacementControlEnergyReleaseConstraint&>(rOther);
            this->mp_ext_f = &thisValue.GetForceVector();
        }
        catch(std::bad_cast& bc)
        {}
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
                    f += 0.5 * fv * (dof_iterator->GetSolutionStepValue() - dof_iterator->GetSolutionStepValue(1));
                }
            }
        }

        std::cout << "ArcLengthDisplacementControlEnergyReleaseConstraint: f = " << f << std::endl;

        return f - BaseType::Radius();
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
                    TSparseSpaceType::SetValue(dfdu, row, 0.5 * fv);
                }
            }
        }

        return dfdu;
    }

    /// Get the derivatives of the constraint
    double GetDerivativesDLambda() const override
    {
        return 0.0;
    }

    /// Compute a trial solution at the predictor stage
    double Predict(const TSystemVectorType& rDeltaUl, const int& rMode) const override
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

        if (rMode != 0)
        {
            if (rMode == -1)
            {
                s0 *= -1.0;
                std::cout << "ArcLengthDisplacementControlEnergyReleaseConstraint: forward sign is forced to be reversed" << std::endl;
            }

            if (rMode == 1)
            {
                std::cout << "ArcLengthDisplacementControlEnergyReleaseConstraint: forward sign is forced to remain" << std::endl;
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
                std::cout << "ArcLengthDisplacementControlEnergyReleaseConstraint: forward sign is reversed" << std::endl;
            }
            else
            {
                std::cout << "ArcLengthDisplacementControlEnergyReleaseConstraint: forward sign remains" << std::endl;
            }
        }

        // compute delta lambda
        double delta_lambda = 1.0/s0 * BaseType::Radius();
        return delta_lambda;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ArcLengthDisplacementControlEnergyReleaseConstraint";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ArcLengthDisplacementControlEnergyReleaseConstraint";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Radius: " << BaseType::Radius() << std::endl;
        BaseType::PrintData(rOStream);
    }

private:

    const TSystemVectorType* mp_ext_f;  // pointer to hold external force vector

}; // Class ArcLengthDisplacementControlEnergyReleaseConstraint

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_DISPLACEMENT_CONTROL_ENERGY_RELEASE_CONSTRAINT_H_INCLUDED defined
