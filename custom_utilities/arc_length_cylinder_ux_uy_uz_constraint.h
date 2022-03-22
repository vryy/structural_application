//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 10/3/2022 $
//
//


#if !defined(KRATOS_ARC_LENGTH_CYLINDER_UX_UY_UZ_CONSTRAINT_H_INCLUDED )
#define  KRATOS_ARC_LENGTH_CYLINDER_UX_UY_UZ_CONSTRAINT_H_INCLUDED


#include "custom_processes/arc_length_control_process.h"


namespace Kratos
{


/**
 * Arc-length cylinder constraint
 * In this constraint, only the free displacement d.o.fs are accounted
 * Reference:
 * +    Souza Neto, Computational Plasticity
 */
template<class TBuilderAndSolverType>
class ArcLengthCylinderUxUyUzConstraint : public ArcLengthConstraint<TBuilderAndSolverType>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthCylinderUxUyUzConstraint );

    typedef ArcLengthConstraint<TBuilderAndSolverType> BaseType;

    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    ArcLengthCylinderUxUyUzConstraint(const double& Radius)
    : BaseType(), mRadius(Radius)
    {
        std::cout << "ArcLengthCylinderUxUyUzConstraint is used, radius = " << mRadius << std::endl;
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

        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if ( dof_iterator->GetVariable() == DISPLACEMENT_X
              || dof_iterator->GetVariable() == DISPLACEMENT_Y
              || dof_iterator->GetVariable() == DISPLACEMENT_Z )
            {
                const auto row = dof_iterator->EquationId();
                if (row < EquationSystemSize)
                    f += pow(dof_iterator->GetSolutionStepValue() - dof_iterator->GetSolutionStepValue(1), 2);
            }
        }

        return sqrt(f) - mRadius;
    }

    /// Get the derivatives of the constraint
    TSystemVectorType GetDerivativesDU() const override
    {
        const auto EquationSystemSize = this->GetBuilderAndSolver().GetEquationSystemSize();
        const DofsArrayType& rDofSet = this->GetBuilderAndSolver().GetDofSet();

        TSystemVectorType dfdu;
        TSparseSpaceType::Resize(dfdu, EquationSystemSize);
        TSparseSpaceType::SetToZero(dfdu);

        double f = this->GetValue() + mRadius;

        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if ( dof_iterator->GetVariable() == DISPLACEMENT_X
              || dof_iterator->GetVariable() == DISPLACEMENT_Y
              || dof_iterator->GetVariable() == DISPLACEMENT_Z )
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

        if (rMode != 0)
        {
            if (rMode == -1)
            {
                s0 *= -1.0;
                std::cout << "ArcLengthCylinderUxUyUzConstraint: forward sign is forced to be reversed" << std::endl;
            }

            if (rMode == 1)
            {
                std::cout << "ArcLengthCylinderUxUyUzConstraint: forward sign is forced to remain" << std::endl;
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
                std::cout << "ArcLengthCylinderUxUyUzConstraint: forward sign is reversed" << std::endl;
            }
            else
            {
                std::cout << "ArcLengthCylinderUxUyUzConstraint: forward sign remains" << std::endl;
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
        return "ArcLengthCylinderUxUyUzConstraint";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ArcLengthCylinderUxUyUzConstraint";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Radius: " << mRadius << std::endl;
        BaseType::PrintData(rOStream);
    }

private:

    double mRadius;

}; // Class ArcLengthCylinderUxUyUzConstraint

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_CYLINDER_UX_UY_UZ_CONSTRAINT_H_INCLUDED defined
