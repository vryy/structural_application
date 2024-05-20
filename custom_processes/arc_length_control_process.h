//   Project Name:        Kratos
//   Modified by:         $Author: hbui $
//   Date:                $Date: 8/3/2022 $
//
//


#if !defined(KRATOS_ARC_LENGTH_CONTROL_PROCESS_H_INCLUDED )
#define  KRATOS_ARC_LENGTH_CONTROL_PROCESS_H_INCLUDED


#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_utilities/arc_length_constraint.h"


namespace Kratos
{


/**
 * Arc-length control process using the consistent scheme to enforce the constraint
 * Reference:
 * +    [1] GM lecture note on solution method
 * +    [2] Souza de Neto, Computational Plasticity
 */
template<class TBuilderAndSolverType>
class ArcLengthControlProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ArcLengthControlProcess );

    typedef typename TBuilderAndSolverType::TSparseSpaceType TSparseSpaceType;
    typedef typename TBuilderAndSolverType::TSystemMatrixType TSystemMatrixType;
    typedef typename TBuilderAndSolverType::TSystemVectorType TSystemVectorType;
    typedef ArcLengthConstraint<TBuilderAndSolverType> ArcLengthConstraintType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ArcLengthControlProcess(typename ArcLengthConstraintType::Pointer pConstraint)
    : mLambda(0.0), mLambdaOld(0.0), mLambdaOldOld(0.0)
    , mIsPredictor(false), mForcedFlag(0), mSolveMode(0)
    , mp_delta_u_r(NULL), mp_delta_u_l(NULL)
    , mpConstraint(pConstraint)
    {
        mpConstraint->SetLambda(mLambda);
        mpConstraint->SetLambdaOld(mLambdaOld);
        mpConstraint->SetLambdaOldOld(mLambdaOldOld);
    }

    /// Copy constructor
    ArcLengthControlProcess(const ArcLengthControlProcess& rOther)
    : mLambda(rOther.mLambda), mLambdaOld(rOther.mLambdaOld), mLambdaOldOld(rOther.mLambdaOldOld)
    , mpConstraint(rOther.mpConstraint)
    {
        mpConstraint->SetLambda(mLambda);
        mpConstraint->SetLambdaOld(mLambdaOld);
        mpConstraint->SetLambdaOldOld(mLambdaOldOld);
    }

    ///@}
    ///@name Access
    ///@{

    /// Get the constraint
    typename ArcLengthConstraintType::Pointer pGetConstraint() const
    {
        return mpConstraint;
    }

    /// Get the multiplier
    double GetLambda() const
    {
        return mLambda;
    }

    /// Get the multiplier
    double GetLambdaOld() const
    {
        return mLambdaOld;
    }

    /// Get the multiplier increment
    double GetDeltaLambda() const
    {
        return mLambda - mLambdaOld;
    }

    /// Get the multiplier increment
    double GetDeltaLambdaOld() const
    {
        return mLambdaOld - mLambdaOldOld;
    }

    /// Set vector delta_u_r
    void SetDeltaUr(TSystemVectorType& Dx)
    {
        mp_delta_u_r = &Dx;
    }

    /// Set vector delta_u_l
    void SetDeltaUl(TSystemVectorType& Dx)
    {
        mp_delta_u_l = &Dx;
    }

    /// Set the predictor flag
    void SetPredictor(const bool rValue)
    {
        mIsPredictor = rValue;
    }

    /// Set the force flag
    void SetForcedMode(const int rValue)
    {
        mForcedFlag = rValue;
    }

    /// Set the solve option
    void SetSolveMode(const int rValue)
    {
        mSolveMode = rValue;
    }

    /// Get the solve option
    int GetSolveMode() const
    {
        return mSolveMode;
    }

    /// Set the model_part
    void SetModelPart(const ModelPart& r_model_part)
    {
        mpConstraint->SetModelPart(r_model_part);
    }

    /// Set the builder_and_solver
    void SetBuilderAndSolver(const TBuilderAndSolverType& r_builder_and_solver)
    {
        mpConstraint->SetBuilderAndSolver(r_builder_and_solver);
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    ArcLengthControlProcess& operator=(const ArcLengthControlProcess& rOther)
    {
        this->Copy(rOther);
    }

    // Copy data from other process
    void Copy(const ArcLengthControlProcess& rOther)
    {
        mLambda = rOther.mLambda;
        mLambdaOld = rOther.mLambdaOld;
        mLambdaOldOld = rOther.mLambdaOldOld;
        mpConstraint->Copy(*rOther.mpConstraint);
    }

    ///@}
    ///@name Operations
    ///@{

    /// Reset the multiplier
    void Reset()
    {
        mLambda = 0.0;
        mLambdaOld = 0.0;
        mLambdaOldOld = 0.0;
    }

    /// Update the constraint with new multiplier
    void Update(const double& DeltaLambda)
    {
        mLambda += DeltaLambda;
    }

    /// To be called at the beginning of the solution step
    void ExecuteInitializeSolutionStep() override
    {}

    /// Compute the new solution based on delta_u_r and delta_u_l
    /// On output delta_u_r will be overrided by delta_u
    virtual void Execute(TSystemVectorType& rDeltaUr, const TSystemVectorType& rDeltaUl)
    {
        double delta_lambda;

        if (this->IsPredictor())
        {
            // compute initial value of delta lambda
            delta_lambda = this->GetConstraint().Predict(rDeltaUl, this->ForcedMode());

            // reset delta u
            TSparseSpaceType::SetToZero(rDeltaUr);
        }
        else
        {
            // compute delta lambda
            delta_lambda = this->GetConstraint().Solve(rDeltaUr, rDeltaUl, mSolveMode);
        }

        // update delta u
        TSparseSpaceType::UnaliasedAdd(rDeltaUr, delta_lambda, rDeltaUl);

        // update lambda
        this->Update(delta_lambda);
    }

    /// Compute the new solution
    void Execute() override
    {
        if (mp_delta_u_r == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Delta_U_r is not set", "")

        if (mp_delta_u_l == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Delta_U_l is not set", "")

        this->Execute(*mp_delta_u_r, *mp_delta_u_l);
    }

    /// To be called at the end of the solution step
    void ExecuteFinalizeSolutionStep() override
    {
        // record current multiplier
        mLambdaOldOld = mLambdaOld;
        mLambdaOld = mLambda;

        // reset temporary vectors
        mp_delta_u_r = NULL;
        mp_delta_u_l = NULL;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ArcLengthControlProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info() << ", Constraint: ";
        mpConstraint->PrintInfo(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Lambda: " << mLambda
                 << ", Lambda_old: " << mLambdaOld
                 << ", Lambda_old_old: " << mLambdaOldOld;
        rOStream << " Constraint Data:";
        mpConstraint->PrintData(rOStream);
    }

protected:

    ArcLengthConstraintType& GetConstraint()
    {
        return *mpConstraint;
    }

    const ArcLengthConstraintType& GetConstraint() const
    {
        return *mpConstraint;
    }

    bool IsPredictor() const
    {
        return mIsPredictor;
    }

    int ForcedMode() const
    {
        return mForcedFlag;
    }

    double Lambda() const
    {
        return mLambda;
    }

    double LambdaOld() const
    {
        return mLambdaOld;
    }

    double LambdaOldOld() const
    {
        return mLambdaOldOld;
    }

private:

    double mLambda;
    double mLambdaOld;
    double mLambdaOldOld;

    bool mIsPredictor;
    int mForcedFlag; // this flag is to allow user to intervene the arc-length control by forcing the reverse of loading. Use it with care.
                     // 0: no forcing; -1: force reverse; 1: force forward
    int mSolveMode;  // solve option for delta_lambda; 0: consistent scheme; 1: non-consistent scheme (2)

    TSystemVectorType* mp_delta_u_r;
    TSystemVectorType* mp_delta_u_l;

    typename ArcLengthConstraintType::Pointer mpConstraint;

}; // Class ArcLengthControlProcess

///@}
///@name Input and output
///@{

/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                   TrussElement& rThis);
*/

/// output stream function
template<class TBuilderAndSolverType>
inline std::ostream& operator << (std::ostream& rOStream, const ArcLengthControlProcess<TBuilderAndSolverType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_CONTROL_PROCESS_H_INCLUDED defined
