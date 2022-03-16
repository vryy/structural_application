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


namespace Kratos
{


/**
 * Abstract class for arc-length control
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

    typedef typename TBuilderAndSolverType::TSchemeType TSchemeType;
    typedef typename TSchemeType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ArcLengthControlProcess()
    : mLambda(0.0), mLambdaOld(0.0), mLambdaOldOld(0.0)
    , mIsPredictor(false), mIsForcedReverse(false), mIsForcedForward(false)
    , mp_delta_u_r(NULL), mp_delta_u_l(NULL)
    , mp_model_part(NULL), mp_builder_and_solver(NULL)
    , mp_ext_f(NULL)
    {}

    /// Copy constructor
    ArcLengthControlProcess(const ArcLengthControlProcess& rOther)
    : mLambda(rOther.mLambda), mLambdaOld(rOther.mLambdaOld), mLambdaOldOld(rOther.mLambdaOldOld)
    , mp_model_part(rOther.mp_model_part), mp_builder_and_solver(rOther.mp_builder_and_solver)
    {}

    ///@}
    ///@name Access
    ///@{

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
    void SetPredictor(const bool& rValue)
    {
        mIsPredictor = rValue;
    }

    /// Set the force reverse flag
    void SetForcedReverse(const bool& rValue)
    {
        mIsForcedReverse = rValue;
    }

    /// Set the force forward flag
    void SetForcedForward(const bool& rValue)
    {
        mIsForcedForward = rValue;
    }

    /// Set the model_part
    void SetModelPart(const ModelPart& r_model_part)
    {
        mp_model_part = &r_model_part;
    }

    /// Set the builder_and_solver
    void SetBuilderAndSolver(const TBuilderAndSolverType& r_builder_and_solver)
    {
        mp_builder_and_solver = &r_builder_and_solver;
    }

    /// Set the external force vector
    void SetForceVector(const TSystemVectorType& rValue)
    {
        mp_ext_f = &rValue;
    }

    /// Set the stiffness matrix
    void SetTangentMatrix(const TSystemMatrixType& rValue)
    {
        mp_Kt = &rValue;
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
        mp_model_part = rOther.mp_model_part;
        mp_builder_and_solver = rOther.mp_builder_and_solver;
        mp_ext_f = rOther.mp_ext_f;
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

    /// Compute a trial solution at the predictor state
    virtual double Predict(const TSystemVectorType& rDeltaUl) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    /// To be called at the beginning of the solution step
    void ExecuteInitializeSolutionStep() override
    {}

    /// Compute the new solution based on delta_u_r and delta_u_l
    /// On output delta_u_r will be overrided by delta_u
    virtual void Execute(TSystemVectorType& rDeltaUr, const TSystemVectorType& rDeltaUl)
    {
        double delta_lambda;

        if (mIsPredictor)
        {
            // compute initial value of delta lambda
            delta_lambda = this->Predict(rDeltaUl);

            // reset delta u
            TSparseSpaceType::SetToZero(rDeltaUr);
        }
        else
        {
            // compute delta lambda
            double f = this->GetValue();
            TSystemVectorType dfdu = this->GetDerivativesDU();
            double dfdl = this->GetDerivativesDLambda();
            delta_lambda = -(f + TSparseSpaceType::Dot(dfdu, rDeltaUr)) / (dfdl + TSparseSpaceType::Dot(dfdu, rDeltaUl));
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
        rOStream << "ArcLengthControlProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Lambda: " << mLambda
                 << ", Lambda_old: " << mLambdaOld
                 << ", Lambda_old_old: " << mLambdaOldOld;
    }

protected:

    const ModelPart& GetModelPart() const
    {
        if (mp_model_part == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "ModelPart is not set", "")
        return *mp_model_part;
    }

    const TBuilderAndSolverType& GetBuilderAndSolver() const
    {
        if (mp_builder_and_solver == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "BuilderAndSolver is not set", "")
        return *mp_builder_and_solver;
    }

    const bool& IsForcedReverse() const
    {
        return mIsForcedReverse;
    }

    const bool& IsForcedForward() const
    {
        return mIsForcedForward;
    }

    const double& Lambda() const
    {
        return mLambda;
    }

    const double& LambdaOld() const
    {
        return mLambdaOld;
    }

    const double& LambdaOldOld() const
    {
        return mLambdaOldOld;
    }

    const TSystemVectorType& GetForceVector() const
    {
        if (mp_ext_f == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Force vector is not set", "")
        return *mp_ext_f;
    }

    const TSystemMatrixType& GetTangentMatrix() const
    {
        if (mp_Kt == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Force vector is not set", "")
        return *mp_Kt;
    }

private:

    const ModelPart* mp_model_part;
    const TBuilderAndSolverType* mp_builder_and_solver;

    double mLambda;
    double mLambdaOld;
    double mLambdaOldOld;

    bool mIsPredictor;
    bool mIsForcedReverse; // this flag is to allow user to intervene the arc-length control by forcing the reverse of loading. Use it with care.
    bool mIsForcedForward;

    TSystemVectorType* mp_delta_u_r;
    TSystemVectorType* mp_delta_u_l;

    const TSystemVectorType* mp_ext_f;  // pointer to hold external force vector if necessary
    const TSystemMatrixType* mp_Kt;     // pointer to hold the system stiffness/tangent matrix

}; // Class ArcLengthControlProcess

}  // namespace Kratos.

#endif // KRATOS_ARC_LENGTH_CONTROL_PROCESS_H_INCLUDED defined
