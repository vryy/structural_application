//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 10/3/2022 $
//
//


#if !defined(KRATOS_ARC_LENGTH_CYLINDER_UX_UY_UZ_CONSTRAINT_H_INCLUDED )
#define  KRATOS_ARC_LENGTH_CYLINDER_UX_UY_UZ_CONSTRAINT_H_INCLUDED


#include "custom_processes/arc_length_control_process.h"
#include "custom_utilities/sd_math_utils.h"


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

    ArcLengthCylinderUxUyUzConstraint(const double Radius)
    : BaseType(Radius)
    {
        std::cout << "ArcLengthCylinderUxUyUzConstraint is used, radius = " << BaseType::Radius() << std::endl;
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
            if ( dof_iterator->GetVariable() == DISPLACEMENT_X
              || dof_iterator->GetVariable() == DISPLACEMENT_Y
              || dof_iterator->GetVariable() == DISPLACEMENT_Z )
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
        double delta_lambda = 1.0/s0 * BaseType::Radius();
        return delta_lambda;
    }

    /// Solve for delta_lambda for the non-consistent scheme
    /// It is noted that the input taking in the incremental solution (u_(k+1)-u_k), not the delta solution (u_(k+1)-u_0)
    double SolveNonConsistent(const TSystemVectorType& r_d_ur, const TSystemVectorType& r_d_ul) const override
    {
        const auto EquationSystemSize = this->GetBuilderAndSolver().GetEquationSystemSize();
        const DofsArrayType& rDofSet = this->GetBuilderAndSolver().GetDofSet();

        // assemble current delta u
        TSystemVectorType Du0;
        TSystemVectorType Du1;
        TSparseSpaceType::Resize(Du0, EquationSystemSize);
        TSparseSpaceType::Resize(Du1, EquationSystemSize);
        TSparseSpaceType::SetToZero(Du0);
        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            if ( dof_iterator->GetVariable() == DISPLACEMENT_X
              || dof_iterator->GetVariable() == DISPLACEMENT_Y
              || dof_iterator->GetVariable() == DISPLACEMENT_Z )
            {
                const auto row = dof_iterator->EquationId();
                if (row < EquationSystemSize)
                    TSparseSpaceType::SetValue(Du0, row, dof_iterator->GetSolutionStepValue() - dof_iterator->GetSolutionStepValue(1));
            }
        }
        TSparseSpaceType::ScaleAndAdd(1.0, Du0, 1.0, r_d_ur, Du1);

        // compute current delta_lambda
        double delta_lambda = this->Lambda() - this->LambdaOld();

        double a = TSparseSpaceType::Dot(r_d_ul, r_d_ul);
        double b = 2.0 * TSparseSpaceType::Dot(Du1, r_d_ul);
        double c = TSparseSpaceType::Dot(Du1, Du1) - BaseType::Radius()*BaseType::Radius();

        double x[2];
        int solve_flag = SD_MathUtils<double>::SolveQuadratic(a, b, c, x);
        if (solve_flag == 0)
        {
            KRATOS_WATCH(a)
            KRATOS_WATCH(b)
            KRATOS_WATCH(c)
            KRATOS_ERROR << "The arc-length non-consistent scheme encounters complex solution";
        }
        else
        {
            if (solve_flag == 1)
            {
                return x[0];
            }
            else if (solve_flag == 2)
            {
                // select the appropriate root
                TSystemVectorType Du21, Du22;
                TSparseSpaceType::Resize(Du21, EquationSystemSize);
                TSparseSpaceType::Resize(Du22, EquationSystemSize);
                TSparseSpaceType::ScaleAndAdd(1.0, Du1, x[0], r_d_ul, Du21);
                TSparseSpaceType::ScaleAndAdd(1.0, Du1, x[1], r_d_ul, Du22);

                double v1 = TSparseSpaceType::Dot(Du21, Du0);
                double v2 = TSparseSpaceType::Dot(Du22, Du0);

                if (v1 > v2)
                    return x[0];
                else
                    return x[1];
            }
            else
                KRATOS_ERROR << "Invalid solve flag " << solve_flag;
        }

        // can't come here, to satisfy the compiler
        return 0.0;
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
        rOStream << " Radius: " << BaseType::Radius() << std::endl;
        BaseType::PrintData(rOStream);
    }

}; // Class ArcLengthCylinderUxUyUzConstraint

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_CYLINDER_UX_UY_UZ_CONSTRAINT_H_INCLUDED defined
