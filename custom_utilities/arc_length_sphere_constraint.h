//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 8/3/2022 $
//
//


#if !defined(KRATOS_ARC_LENGTH_SPHERE_CONSTRAINT_H_INCLUDED )
#define  KRATOS_ARC_LENGTH_SPHERE_CONSTRAINT_H_INCLUDED


#include "custom_processes/arc_length_control_process.h"
#include "custom_utilities/sd_math_utils.h"


namespace Kratos
{


/**
 * Arc-length sphere control
 * In this constraint, all the free d.o.fs are accounted
 * Reference:
 * +    Crisfield 1981
 */
template<class TBuilderAndSolverType>
class ArcLengthSphereConstraint : public ArcLengthConstraint<TBuilderAndSolverType>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthSphereConstraint );

    typedef ArcLengthConstraint<TBuilderAndSolverType> BaseType;

    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    ArcLengthSphereConstraint(const double& Psi, const double& Radius)
    : BaseType(Radius), mPsi(Psi)
    {
        std::cout << "ArcLengthSphereConstraint is used, radius = " << BaseType::Radius() << ", load scale factor = " << mPsi << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    void SetScale(const double& Psi)
    {
        mPsi = Psi;
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    ArcLengthSphereConstraint& operator=(const ArcLengthSphereConstraint& rOther)
    {
        BaseType::operator=(rOther);
        mPsi = rOther.mPsi;
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
            const auto row = dof_iterator->EquationId();
            if (row < EquationSystemSize)
                f += pow(dof_iterator->GetSolutionStepValue() - dof_iterator->GetSolutionStepValue(1), 2);
        }

        f += pow(mPsi * (this->Lambda() - this->LambdaOld()), 2);

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
            const auto row = dof_iterator->EquationId();
            if (row < EquationSystemSize)
                TSparseSpaceType::SetValue(dfdu, row, dof_iterator->GetSolutionStepValue() - dof_iterator->GetSolutionStepValue(1));
        }

        TSparseSpaceType::InplaceMult(dfdu, 1.0/f);
        return dfdu;
    }

    /// Get the derivatives of the constraint
    double GetDerivativesDLambda() const override
    {
        double f = this->GetValue() + BaseType::Radius();
        return mPsi*mPsi*(this->Lambda() - this->LambdaOld()) / f;
    }

    /// Compute a trial solution at the predictor stage
    double Predict(const TSystemVectorType& rDeltaUl, const int& rMode) const override
    {
        const auto EquationSystemSize = this->GetBuilderAndSolver().GetEquationSystemSize();
        const DofsArrayType& rDofSet = this->GetBuilderAndSolver().GetDofSet();

        // compute s0
        double norm_dul = TSparseSpaceType::TwoNorm(rDeltaUl);
        double s0 = sqrt(mPsi*mPsi + norm_dul*norm_dul);
        // KRATOS_WATCH(s0)

        /* compute the forward criteria */
        // assemble delta u
        TSystemVectorType Du;
        TSparseSpaceType::Resize(Du, EquationSystemSize);
        TSparseSpaceType::SetToZero(Du);
        for (typename DofsArrayType::const_iterator dof_iterator = rDofSet.begin(); dof_iterator != rDofSet.end(); ++dof_iterator)
        {
            const auto row = dof_iterator->EquationId();
            if (row < EquationSystemSize)
                TSparseSpaceType::SetValue(Du, row, dof_iterator->GetSolutionStepValue(1) - dof_iterator->GetSolutionStepValue(2));
        }

        if (rMode != 0)
        {
            if (rMode == -1)
            {
                s0 *= -1.0;
                std::cout << "ArcLengthSphereConstraint: forward sign is forced to be reversed" << std::endl;
            }

            if (rMode == 1)
            {
                std::cout << "ArcLengthSphereConstraint: forward sign is forced to remain" << std::endl;
            }

            if (rMode == 2)
            {
                // compute the forward criteria
                double forward_criteria = TSparseSpaceType::Dot(Du, rDeltaUl) + mPsi*mPsi*(this->LambdaOld() - this->LambdaOldOld());

                if (forward_criteria < -1.0e-10)
                {
                    std::cout << "ArcLengthSphereConstraint: forward sign is checked to be reversed. Proceed? (y/n)" << std::endl;

                    char c;
                    std::cin >> c;

                    if (c == 'y')
                    {
                        s0 *= -1.0;
                       std::cout << "ArcLengthSphereConstraint: forward sign is proceeded to be reversed" << std::endl;
                    }
                    else
                    {
                       std::cout << "ArcLengthSphereConstraint: forward sign is not proceeded to be reversed" << std::endl;
                    }
                }

            }
        }
        else
        {
            // compute the forward criteria
            double forward_criteria = TSparseSpaceType::Dot(Du, rDeltaUl) + mPsi*mPsi*(this->LambdaOld() - this->LambdaOldOld());
            // KRATOS_WATCH(Du)
            // KRATOS_WATCH(Dx)
            // KRATOS_WATCH(TSparseSpaceType::Dot(Du, Dx))
            // KRATOS_WATCH(this->LambdaOld())
            // KRATOS_WATCH(this->LambdaOldOld())
            // KRATOS_WATCH(forward_criteria)

            if (forward_criteria < -1.0e-10)
            {
                s0 *= -1.0;
                std::cout << "ArcLengthSphereConstraint: forward sign is reversed" << std::endl;
            }
            else
            {
                std::cout << "ArcLengthSphereConstraint: forward sign remains" << std::endl;
            }
        }

        // compute delta lambda
        double delta_lambda = mPsi/s0 * BaseType::Radius();
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
            const auto row = dof_iterator->EquationId();
            if (row < EquationSystemSize)
                TSparseSpaceType::SetValue(Du0, row, dof_iterator->GetSolutionStepValue() - dof_iterator->GetSolutionStepValue(1));
        }
        TSparseSpaceType::ScaleAndAdd(1.0, Du0, 1.0, r_d_ur, Du1);

        // compute current delta_lambda
        double delta_lambda = this->Lambda() - this->LambdaOld();

        double a = TSparseSpaceType::Dot(r_d_ul, r_d_ul) + mPsi*mPsi;
        double b = 2.0 * (TSparseSpaceType::Dot(Du1, r_d_ul) + mPsi*mPsi*delta_lambda);
        double c = TSparseSpaceType::Dot(Du1, Du1) + mPsi*mPsi*delta_lambda*delta_lambda - BaseType::Radius()*BaseType::Radius();

        double x[2];
        int solve_flag = SD_MathUtils<double>::SolveQuadratic(a, b, c, x);
        if (solve_flag == 0)
        {
            KRATOS_WATCH(a)
            KRATOS_WATCH(b)
            KRATOS_WATCH(c)
            KRATOS_THROW_ERROR(std::logic_error, "The arc-length non-consistent scheme encounters complex solution", "")
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

                double v1 = TSparseSpaceType::Dot(Du21, Du0) + mPsi*mPsi*delta_lambda*(delta_lambda + x[0]);
                double v2 = TSparseSpaceType::Dot(Du22, Du0) + mPsi*mPsi*delta_lambda*(delta_lambda + x[1]);

                // double norm_0 = TSparseSpaceType::TwoNorm(Du0);
                // double norm_1 = TSparseSpaceType::TwoNorm(Du21);
                // double norm_2 = TSparseSpaceType::TwoNorm(Du22);
                // double v1 = TSparseSpaceType::Dot(Du21, Du0) / (norm_0*norm_1);
                // double v2 = TSparseSpaceType::Dot(Du22, Du0) / (norm_0*norm_2);

                if (v1 > v2)
                    return x[0];
                else
                    return x[1];
            }
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ArcLengthSphereConstraint";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ArcLengthSphereConstraint";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Radius: " << BaseType::Radius()
                 << ", Psi: " << mPsi << std::endl;
        BaseType::PrintData(rOStream);
    }

private:

    double mPsi;

}; // Class ArcLengthSphereConstraint

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_SPHERE_CONSTRAINT_H_INCLUDED defined
