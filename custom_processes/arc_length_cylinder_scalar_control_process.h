//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 11/3/2022 $
//
//


#if !defined(KRATOS_ARC_LENGTH_CYLINDER_SCALAR_CONTROL_PROCESS_H_INCLUDED )
#define  KRATOS_ARC_LENGTH_CYLINDER_SCALAR_CONTROL_PROCESS_H_INCLUDED


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
template<class TBuilderAndSolveType>
class ArcLengthCylinderScalarControlProcess : public ArcLengthControlProcess<TBuilderAndSolveType>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthCylinderScalarControlProcess );

    typedef ArcLengthControlProcess<TBuilderAndSolveType> BaseType;

    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    ArcLengthCylinderScalarControlProcess(const Variable<double>& rVariable, const double& Radius)
    : BaseType(), mrVariable(rVariable), mRadius(Radius)
    {
        std::cout << "ArcLengthCylinderScalarControlProcess is used, variable = " << mrVariable.Name() << ", radius = " << mRadius << std::endl;
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
            if ( dof_iterator->GetVariable() == mrVariable )
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

        if (this->IsForcedReverse() || this->IsForcedForward())
        {
            if (this->IsForcedReverse())
            {
                s0 *= -1.0;
                std::cout << "ArcLengthCylinderScalarControlProcess: forward sign is forced to be reversed" << std::endl;
            }

            if (this->IsForcedForward())
            {
                std::cout << "ArcLengthCylinderScalarControlProcess: forward sign is forced to remain" << std::endl;
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
                std::cout << "ArcLengthCylinderScalarControlProcess: forward sign is reversed" << std::endl;
            }
            else
            {
                std::cout << "ArcLengthCylinderScalarControlProcess: forward sign remains" << std::endl;
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
        return "ArcLengthCylinderScalarControlProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ArcLengthCylinderScalarControlProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Variable: " << mrVariable.Name();
        rOStream << ", Radius: " << mRadius << std::endl;
        BaseType::PrintData(rOStream);
    }

private:

    const Variable<double>& mrVariable;
    double mRadius;

}; // Class ArcLengthCylinderScalarControlProcess

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_CYLINDER_SCALAR_CONTROL_PROCESS_H_INCLUDED defined