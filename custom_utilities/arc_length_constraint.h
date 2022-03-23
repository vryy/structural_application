//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 15/3/2022 $
//
//


#if !defined(KRATOS_ARC_LENGTH_CONSTRAINT_H_INCLUDED )
#define  KRATOS_ARC_LENGTH_CONSTRAINT_H_INCLUDED


#include "includes/define.h"
#include "includes/model_part.h"


namespace Kratos
{


/**
 * Abstract class for arc-length constraint
 */
template<class TBuilderAndSolverType>
class ArcLengthConstraint
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthConstraint );

    typedef typename TBuilderAndSolverType::TSparseSpaceType TSparseSpaceType;
    typedef typename TBuilderAndSolverType::TSystemMatrixType TSystemMatrixType;
    typedef typename TBuilderAndSolverType::TSystemVectorType TSystemVectorType;

    typedef typename TBuilderAndSolverType::TSchemeType TSchemeType;
    typedef typename TSchemeType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ArcLengthConstraint(const double& Radius)
    : mRadius(Radius), mpLambda(NULL), mpLambdaOld(NULL), mpLambdaOldOld(NULL)
    , mp_model_part(NULL), mp_builder_and_solver(NULL)
    {}

    /// Copy constructor
    ArcLengthConstraint(ArcLengthConstraint& rOther)
    : mRadius(rOther.mRadius), mpLambda(rOther.mpLambda), mpLambdaOld(rOther.mpLambdaOld), mpLambdaOldOld(rOther.mpLambdaOldOld)
    , mp_model_part(rOther.mp_model_part), mp_builder_and_solver(rOther.mp_builder_and_solver)
    {}

    ///@}
    ///@name Access
    ///@{

    /// Set the arc length
    void SetRadius(const double& Radius)
    {
        mRadius = Radius;
    }

    /// Set the multiplier
    void SetLambda(const double& Lambda)
    {
        mpLambda = &Lambda;
    }

    /// Set the multiplier
    void SetLambdaOld(const double& Lambda)
    {
        mpLambdaOld = &Lambda;
    }

    /// Set the multiplier
    void SetLambdaOldOld(const double& Lambda)
    {
        mpLambdaOldOld = &Lambda;
    }

    /// Get the model_part
    const ModelPart& GetModelPart() const
    {
        if (mp_model_part == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "ModelPart is not set", "")
        return *mp_model_part;
    }

    /// Set the model_part
    void SetModelPart(const ModelPart& r_model_part)
    {
        mp_model_part = &r_model_part;
    }

    /// Get the builder_and_solver
    const TBuilderAndSolverType& GetBuilderAndSolver() const
    {
        if (mp_builder_and_solver == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "BuilderAndSolver is not set", "")
        return *mp_builder_and_solver;
    }

    /// Set the builder_and_solver
    void SetBuilderAndSolver(const TBuilderAndSolverType& r_builder_and_solver)
    {
        mp_builder_and_solver = &r_builder_and_solver;
    }

    /// Get the type of required force vector
    virtual int NeedForceVector() const
    {
        return 0;
    }

    /// Set the external force vector
    virtual void SetForceVector(const TSystemVectorType& rValue)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    /// Get the external force vector
    virtual const TSystemVectorType& GetForceVector() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    ///@}
    ///@name Operators
    ///@{

    /// Copy the internal data from other constraint
    virtual void Copy(const ArcLengthConstraint& rOther)
    {
        this->SetModelPart(rOther.GetModelPart());
        this->SetBuilderAndSolver(rOther.GetBuilderAndSolver());
    }

    ///@}
    ///@name Operations
    ///@{

    /// Solve for delta_lambda, given the incremental solutions
    double Solve(const TSystemVectorType& rDeltaUr, const TSystemVectorType& rDeltaUl, const int& rMode) const
    {
        if (rMode == 0)
        {
            return this->SolveConsistent(rDeltaUr, rDeltaUl);
        }
        else if (rMode == 1)
        {
            return this->SolveNonConsistent(rDeltaUr, rDeltaUl);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid arc-length solve option", rMode)
    }

    /// Solve for delta_lambda for the non-consistent scheme
    virtual double SolveNonConsistent(const TSystemVectorType& rDeltaUr, const TSystemVectorType& rDeltaUl) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling unsupported function", __FUNCTION__)
    }

    /// Get the value of the constraint
    virtual double GetValue() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    /// Get the derivatives of the constraint
    virtual TSystemVectorType GetDerivativesDU() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    /// Get the derivatives of the constraint
    virtual double GetDerivativesDLambda() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    /// Compute a trial solution at the predictor stage
    virtual double Predict(const TSystemVectorType& rDeltaUl, const int& rMode) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ArcLengthConstraint";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ArcLengthConstraint";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

protected:

    const double& Radius() const
    {
        return mRadius;
    }

    const double& Lambda() const
    {
        return *mpLambda;
    }

    const double& LambdaOld() const
    {
        return *mpLambdaOld;
    }

    const double& LambdaOldOld() const
    {
        return *mpLambdaOldOld;
    }

private:

    double mRadius;

    // pointer to the multiplier
    const double* mpLambda;
    const double* mpLambdaOld;
    const double* mpLambdaOldOld;

    // pointer to model_part
    const ModelPart* mp_model_part;

    // pointer to builder_and_solver
    const TBuilderAndSolverType* mp_builder_and_solver;

    /// Solve for delta_lambda for the consistent scheme
    double SolveConsistent(const TSystemVectorType& rDeltaUr, const TSystemVectorType& rDeltaUl) const
    {
        double f = this->GetValue();
        TSystemVectorType dfdu = this->GetDerivativesDU();
        double dfdl = this->GetDerivativesDLambda();
        double delta_lambda = -(f + TSparseSpaceType::Dot(dfdu, rDeltaUr)) / (dfdl + TSparseSpaceType::Dot(dfdu, rDeltaUl));
        return delta_lambda;
    }

}; // Class ArcLengthConstraint

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_CONSTRAINT_H_INCLUDED defined
