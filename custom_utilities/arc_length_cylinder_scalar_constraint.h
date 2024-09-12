//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 11/3/2022 $
//
//


#if !defined(KRATOS_ARC_LENGTH_CYLINDER_SCALAR_CONSTRAINT_H_INCLUDED )
#define  KRATOS_ARC_LENGTH_CYLINDER_SCALAR_CONSTRAINT_H_INCLUDED


#include "containers/variable.h"
#include "custom_processes/arc_length_control_process.h"


namespace Kratos
{


/**
 * Arc-length sphere control
 * In this constraint, only the free d.o.f, defined by the scalar variable, are accounted
 * Reference:
 * +    Souza Neto, Computational Plasticity
 */
template<class TBuilderAndSolverType>
class ArcLengthCylinderScalarConstraint : public ArcLengthConstraint<TBuilderAndSolverType>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthCylinderScalarConstraint );

    typedef ArcLengthConstraint<TBuilderAndSolverType> BaseType;

    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    ArcLengthCylinderScalarConstraint(const Variable<double>& rVariable, const double Radius)
    : BaseType(Radius), mrVariable(rVariable)
    {
        std::cout << "ArcLengthCylinderScalarConstraint is used, variable = " << mrVariable.Name() << ", radius = " << BaseType::Radius() << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Get the value of the constraint
    double GetValue() const override
    {
        const auto EquationSystemSize = this->GetBuilderAndSolver().GetEquationSystemSize();
        const DofsArrayType& rDofSet = this->GetBuilderAndSolver().GetDofSet();

        double f = 0.0;

        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if ( dof_iterator->GetVariable() == mrVariable )
            {
                const auto row = dof_iterator->EquationId();
                if (row < EquationSystemSize)
                    f += pow(dof_iterator->GetSolutionStepValue() - dof_iterator->GetSolutionStepValue(1), 2);
            }
        }

        return sqrt(f) - BaseType::Radius();
    }

    /// Get the derivatives of the constraint
    TSystemVectorType GetDerivativesDU() const override
    {
        const auto EquationSystemSize = this->GetBuilderAndSolver().GetEquationSystemSize();
        const DofsArrayType& rDofSet = this->GetBuilderAndSolver().GetDofSet();

        TSystemVectorType dfdu;
        TSparseSpaceType::Resize(dfdu, EquationSystemSize);
        TSparseSpaceType::SetToZero(dfdu);

        double f = this->GetValue() + BaseType::Radius();

        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if ( dof_iterator->GetVariable() == mrVariable )
            {
                const auto row = dof_iterator->EquationId();
                if (row < EquationSystemSize)
                    TSparseSpaceType::SetValue(dfdu, row, dof_iterator->GetSolutionStepValue() - dof_iterator->GetSolutionStepValue(1));
            }
        }

        TSparseSpaceType::InplaceMult(dfdu, 1.0/f);
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
        // KRATOS_WATCH(rDeltaUl)
        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if ( dof_iterator->GetVariable() == mrVariable )
            {
                const auto row = dof_iterator->EquationId();
                // KRATOS_WATCH(dof_iterator->GetVariable().Name())
                if (row < EquationSystemSize)
                {
                    TSparseSpaceType::SetValue(Du, row, dof_iterator->GetSolutionStepValue(1) - dof_iterator->GetSolutionStepValue(2));

                    double tmp = TSparseSpaceType::GetValue(rDeltaUl, row);
                    // KRATOS_WATCH(tmp)
                    norm_p2_dul += tmp*tmp;
                }
            }
        }

        double s0 = sqrt(norm_p2_dul);

        if (rMode != 0)
        {
            if (rMode == -1)
            {
                s0 *= -1.0;
                std::cout << "ArcLengthCylinderScalarConstraint: forward sign is forced to be reversed" << std::endl;
            }

            if (rMode == 1)
            {
                std::cout << "ArcLengthCylinderScalarConstraint: forward sign is forced to remain" << std::endl;
            }
        }
        else
        {
            // compute the forward criteria
            double forward_criteria = TSparseSpaceType::Dot(Du, rDeltaUl);
            // KRATOS_WATCH(Du)
            // KRATOS_WATCH(rDeltaUl)
            // KRATOS_WATCH(TSparseSpaceType::Dot(Du, rDeltaUl))
            // KRATOS_WATCH(forward_criteria)

            if (forward_criteria < -1.0e-10)
            {
                s0 *= -1.0;
                std::cout << "ArcLengthCylinderScalarConstraint: forward sign is reversed" << std::endl;
            }
            else
            {
                std::cout << "ArcLengthCylinderScalarConstraint: forward sign remains" << std::endl;
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
        return "ArcLengthCylinderScalarConstraint";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ArcLengthCylinderScalarConstraint";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Variable: " << mrVariable.Name();
        rOStream << ", Radius: " << BaseType::Radius() << std::endl;
        BaseType::PrintData(rOStream);
    }

private:

    const Variable<double>& mrVariable;

}; // Class ArcLengthCylinderScalarConstraint

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_CYLINDER_SCALAR_CONSTRAINT_H_INCLUDED defined
