/*
LICENSE: see material_point_application/LICENSE.txt
*/
/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 5 Oct 2021 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/eigen_utility.h"
#include "constitutive_laws/multiplicative_finite_strain_bridging_constitutive_law.h"
#include "structural_application_variables.h"

namespace Kratos
{

//**********************************************************************
MultiplicativeFiniteStrainBridgingConstitutiveLaw::MultiplicativeFiniteStrainBridgingConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//**********************************************************************
MultiplicativeFiniteStrainBridgingConstitutiveLaw::MultiplicativeFiniteStrainBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : ConstitutiveLaw(), mpConstitutiveLaw(pConstitutiveLaw)
{
}

//**********************************************************************
MultiplicativeFiniteStrainBridgingConstitutiveLaw::~MultiplicativeFiniteStrainBridgingConstitutiveLaw()
{
}

//**********************************************************************
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw::Has( const Variable<int>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw::Has( const Variable<double>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw::Has( const Variable<Vector>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool MultiplicativeFiniteStrainBridgingConstitutiveLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
int& MultiplicativeFiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
double& MultiplicativeFiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Vector& MultiplicativeFiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Matrix& MultiplicativeFiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<ConstitutiveLaw::Pointer>& rThisVariable, ConstitutiveLaw::Pointer rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == CONSTITUTIVE_LAW)
        mpConstitutiveLaw = rValue;
    else
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->InitializeMaterial(props, geom, ShapeFunctionsValues);

    const unsigned int dim = geom.WorkingSpaceDimension();

    mLastF = ZeroMatrix(dim, dim);
    mCurrentF = mLastF;
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->ResetMaterial(props, geom, ShapeFunctionsValues);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);

    noalias(mLastF) = mCurrentF;
}

//**********************************************************************
void MultiplicativeFiniteStrainBridgingConstitutiveLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    const GeometryType& rGeometry = rValues.GetElementGeometry();
    const Properties& rProperties = rValues.GetMaterialProperties();
    const unsigned int dim = rGeometry.WorkingSpaceDimension();
    const unsigned int strain_size = dim*(dim+1) / 2;

    const Matrix& F = rValues.GetDeformationGradientF();
    const double& J = rValues.GetDeterminantF();

    // create the infinitesimal strain vector
    Vector StrainVector(strain_size);
    rValues.SetStrainVector(StrainVector);

    // compute elastic left Cauchy-Green tensor from previous step
    // It is noted that the constitutive law must be able to provide ELASTIC_STRAIN_VECTOR without resizing elastic_strain
    Vector elastic_strain(strain_size);
    mpConstitutiveLaw->GetValue(ELASTIC_STRAIN_VECTOR, elastic_strain);

    Matrix elastic_strain_tensor(3, 3), left_cauchy_green_tensor_n(3, 3);
    SD_MathUtils<double>::StrainVectorToTensor(elastic_strain, elastic_strain_tensor);
    EigenUtility::ComputeIsotropicTensorFunction(exp2, left_cauchy_green_tensor_n, elastic_strain_tensor);

    // compute trial elastic left Cauchy-Green tensor
    Matrix Fincr(3, 3), invFn(3, 3);
    double detFn;
    MathUtils<double>::InvertMatrix(mLastF, invFn, detFn);
    noalias(Fincr) = prod(mCurrentF, invFn);

    Matrix left_cauchy_green_tensor_trial(3, 3);
    noalias(left_cauchy_green_tensor_trial) = prod(Fincr, Matrix(prod(left_cauchy_green_tensor_n, trans(Fincr))));

    // compute trial elastic logarithmic strain tensor
    std::vector<double> left_cauchy_green_tensor_trial_pri(3);
    std::vector<Matrix> left_cauchy_green_tensor_trial_eigprj(3);
    EigenUtility::calculate_principle_stresses(
            left_cauchy_green_tensor_trial(0, 0),
            left_cauchy_green_tensor_trial(1, 1),
            left_cauchy_green_tensor_trial(2, 2),
            left_cauchy_green_tensor_trial(0, 1),
            left_cauchy_green_tensor_trial(1, 2),
            left_cauchy_green_tensor_trial(0, 2),
            left_cauchy_green_tensor_trial_pri[0],
            left_cauchy_green_tensor_trial_pri[1],
            left_cauchy_green_tensor_trial_pri[2],
            left_cauchy_green_tensor_trial_eigprj[0],
            left_cauchy_green_tensor_trial_eigprj[1],
            left_cauchy_green_tensor_trial_eigprj[2]);
    ComputeIsotropicTensorFunction(log2, elastic_strain_tensor, left_cauchy_green_tensor_trial_pri, left_cauchy_green_tensor_trial_eigprj);
    SD_MathUtils<double>::StrainTensorToVector(elastic_strain_tensor, StrainVector);

    // integrate the (small strain) constitutive law
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);

    // transform the stress to Cauchy stress
    if (rValues.IsSetStressVector())
    {
        Vector& StressVector = rValues.GetStressVector();
        StressVector /= J;
    }

    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();

        Fourth_Order_Tensor D;
        SD_MathUtils<double>::MatrixToTensor(AlgorithmicTangent, D);

        // Perform extra kinematical operations required by materials of this
        // class at large strains for computation of the spatial modulus 'a'
        // Eq. (14.95), Computational Plasticity, de Souza Neto
        Fourth_Order_Tensor L;
        ComputeDerivativeIsotropicTensorFunction(log2, dlog2, L,
            left_cauchy_green_tensor_trial,
            left_cauchy_green_tensor_trial_pri,
            left_cauchy_green_tensor_trial_eigprj);

        Fourth_Order_Tensor B;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(B);
        const Matrix eye = IdentityMatrix(3);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    for (int l = 0; l < 3; ++l)
                        B[i][j](k, l) = eye(i, k) * left_cauchy_green_tensor_trial(j, l)
                                      + eye(j, k) * left_cauchy_green_tensor_trial(i, l);

        Fourth_Order_Tensor DL;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(DL);
        SD_MathUtils<double>::ProductFourthOrderTensor(1.0/J, D, L, DL);

        Fourth_Order_Tensor A;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
        SD_MathUtils<double>::ProductFourthOrderTensor(1.0, DL, B, A);

        Vector& StressVector = rValues.GetStressVector();
        Matrix StressTensor(3, 3);
        SD_MathUtils<double>::StressVectorToTensor(StressVector, StressTensor);

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    for (int l = 0; l < 3; ++l)
                        A[i][j](k, l) -= StressTensor(i, l) * eye(j, k);

        SD_MathUtils<double>::TensorToMatrix(A, AlgorithmicTangent);
    }
}

//**********************************************************************
int MultiplicativeFiniteStrainBridgingConstitutiveLaw::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    return mpConstitutiveLaw->Check( props, geom, CurrentProcessInfo );
}

} // Namespace Kratos
