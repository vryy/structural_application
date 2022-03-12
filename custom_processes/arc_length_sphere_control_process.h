//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 8/3/2022 $
//
//


#if !defined(KRATOS_ARC_LENGTH_SPHERE_CONTROL_PROCESS_H_INCLUDED )
#define  KRATOS_ARC_LENGTH_SPHERE_CONTROL_PROCESS_H_INCLUDED


#include "custom_processes/arc_length_control_process.h"


namespace Kratos
{


/**
 * Arc-length sphere control
 * In this constraint, all the free d.o.fs are accounted
 * Reference:
 * +    Crisfield 1981
 */
template<class TBuilderAndSolverType>
class ArcLengthSphereControlProcess : public ArcLengthControlProcess<TBuilderAndSolverType>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthSphereControlProcess );

    typedef ArcLengthControlProcess<TBuilderAndSolverType> BaseType;

    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    ArcLengthSphereControlProcess(const double& Psi, const double& Radius)
    : BaseType(), mPsi(Psi), mRadius(Radius)
    {
        std::cout << "ArcLengthCylinderControlProcess is used, radius = " << mRadius << ", load scale factor = " << mPsi << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    void SetScale(const double& Psi)
    {
        mPsi = Psi;
    }

    void SetRadius(const double& Radius)
    {
        mRadius = Radius;
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    ArcLengthSphereControlProcess& operator=(const ArcLengthSphereControlProcess& rOther)
    {
        BaseType::operator=(rOther);
        mPsi = rOther.mPsi;
        mRadius = rOther.mRadius;
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
        double f = this->GetValue() + mRadius;
        return mPsi*mPsi*(this->Lambda() - this->LambdaOld()) / f;
    }

    double Predict(const TSystemVectorType& rDeltaUl) const override
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

        if (this->IsForcedReverse() || this->IsForcedForward())
        {
            if (this->IsForcedReverse())
            {
                s0 *= -1.0;
                std::cout << "ArcLengthSphereControlProcess: forward sign is forced to be reversed" << std::endl;
            }

            if (this->IsForcedForward())
            {
                std::cout << "ArcLengthSphereControlProcess: forward sign is forced to remain" << std::endl;
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
                std::cout << "ArcLengthSphereControlProcess: forward sign is reversed" << std::endl;
            }
            else
            {
                std::cout << "ArcLengthSphereControlProcess: forward sign remains" << std::endl;
            }
        }

        // compute delta lambda
        double delta_lambda = mPsi/s0 * mRadius;
        return delta_lambda;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ArcLengthSphereControlProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ArcLengthSphereControlProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Radius: " << mRadius
                 << ", Psi: " << mPsi << std::endl;
        BaseType::PrintData(rOStream);
    }

private:

    double mPsi;
    double mRadius;

}; // Class ArcLengthSphereControlProcess

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_SPHERE_CONTROL_PROCESS_H_INCLUDED defined
