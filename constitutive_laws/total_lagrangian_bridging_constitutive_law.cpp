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
#include "constitutive_laws/total_lagrangian_bridging_constitutive_law.h"
#include "structural_application_variables.h"

// #define DEBUG_CONSTITUTIVE_LAW

namespace Kratos
{

//**********************************************************************
TotalLagrangianBridgingConstitutiveLaw::TotalLagrangianBridgingConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//**********************************************************************
TotalLagrangianBridgingConstitutiveLaw::TotalLagrangianBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : ConstitutiveLaw(), mpConstitutiveLaw(pConstitutiveLaw)
{
}

//**********************************************************************
TotalLagrangianBridgingConstitutiveLaw::~TotalLagrangianBridgingConstitutiveLaw()
{
}

//**********************************************************************
bool TotalLagrangianBridgingConstitutiveLaw::Has( const Variable<int>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool TotalLagrangianBridgingConstitutiveLaw::Has( const Variable<double>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool TotalLagrangianBridgingConstitutiveLaw::Has( const Variable<Vector>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool TotalLagrangianBridgingConstitutiveLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
int& TotalLagrangianBridgingConstitutiveLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
double& TotalLagrangianBridgingConstitutiveLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Vector& TotalLagrangianBridgingConstitutiveLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Matrix& TotalLagrangianBridgingConstitutiveLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::SetValue( const Variable<ConstitutiveLaw::Pointer>& rThisVariable, ConstitutiveLaw::Pointer rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == CONSTITUTIVE_LAW)
        mpConstitutiveLaw = rValue;
    else
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->InitializeMaterial(props, geom, ShapeFunctionsValues);

    mLastF.resize(3, 3, false);
    mCurrentF.resize(3, 3, false);
    noalias(mLastF) = IdentityMatrix(3);
    noalias(mCurrentF) = mLastF;
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->ResetMaterial(props, geom, ShapeFunctionsValues);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);

    noalias(mLastF) = mCurrentF;
}

//**********************************************************************
void TotalLagrangianBridgingConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    const GeometryType& rGeometry = rValues.GetElementGeometry();
    const Properties& rProperties = rValues.GetMaterialProperties();
    const ProcessInfo& rProcessInfo = rValues.GetProcessInfo();
    const unsigned int dim = rGeometry.WorkingSpaceDimension();
    const unsigned int strain_size = rValues.GetStrainVector().size();

    #ifdef DEBUG_CONSTITUTIVE_LAW
    int ElemId, GaussId;
    mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
    mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, GaussId);
    #endif

    /* compute the current logarithmic strain */
    const Matrix& F = rValues.GetDeformationGradientF();
    if (F.size1() == 3)
    {
        noalias(mCurrentF) = F;
    }
    else if (F.size1() == 2)
    {
        mCurrentF.clear();
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                mCurrentF(i, j) = F(i, j);
        mCurrentF(2, 2) = 1.0;
    }
    else
        KRATOS_ERROR << "Invalid F size";

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        KRATOS_WATCH(mCurrentF)
    }
    #endif

    Matrix right_cauchy_green_tensor(3, 3); // C
    noalias(right_cauchy_green_tensor) = prod(trans(mCurrentF), mCurrentF);

    std::vector<double> right_cauchy_green_tensor_pri(3);
    std::vector<Matrix> right_cauchy_green_tensor_eigprj(3);
    EigenUtility::calculate_principle_stresses(
            right_cauchy_green_tensor(0, 0),
            right_cauchy_green_tensor(1, 1),
            right_cauchy_green_tensor(2, 2),
            right_cauchy_green_tensor(0, 1),
            right_cauchy_green_tensor(1, 2),
            right_cauchy_green_tensor(0, 2),
            right_cauchy_green_tensor_pri[0],
            right_cauchy_green_tensor_pri[1],
            right_cauchy_green_tensor_pri[2],
            right_cauchy_green_tensor_eigprj[0],
            right_cauchy_green_tensor_eigprj[1],
            right_cauchy_green_tensor_eigprj[2]);

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        std::cout << "right_cauchy_green_tensor_pri:";
        for (unsigned int i = 0; i < 3; ++i)
            std::cout << " " << right_cauchy_green_tensor_pri[i];
        std::cout << std::endl;
    }
    #endif

    Matrix current_strain_tensor(3, 3);
    SD_MathUtils<double>::ComputeIsotropicTensorFunction(EigenUtility::logd2, current_strain_tensor,
        right_cauchy_green_tensor_pri, right_cauchy_green_tensor_eigprj);

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        KRATOS_WATCH(current_strain_tensor)
    }
    #endif

    /* compute the right stretch tensor (U) and its inverse. Note that C = U^2 */
    Matrix right_stretch_tensor(3, 3);
    SD_MathUtils<double>::ComputeIsotropicTensorFunction(EigenUtility::sqrt, right_stretch_tensor,
        right_cauchy_green_tensor_pri, right_cauchy_green_tensor_eigprj);

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        KRATOS_WATCH(right_stretch_tensor)
    }
    #endif

    Matrix right_stretch_tensor_inversed(3, 3);
    double detU;
    MathUtils<double>::InvertMatrix(right_stretch_tensor, right_stretch_tensor_inversed, detU);

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        KRATOS_WATCH(right_stretch_tensor_inversed)
        KRATOS_WATCH(detU)
    }
    #endif

    /* compute the rotation tensor (R) */
    Matrix rotation_tensor(3, 3);
    noalias(rotation_tensor) = prod(mCurrentF, right_stretch_tensor_inversed);

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        KRATOS_WATCH(rotation_tensor)
    }
    #endif

    /* compute the corotational logarithmic strain */
    Matrix current_rotated_strain_tensor(3, 3);
    noalias(current_rotated_strain_tensor) = prod(trans(rotation_tensor),
            Matrix(prod(current_strain_tensor, rotation_tensor)));

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == 1)
    {
        KRATOS_WATCH(current_rotated_strain_tensor)
    }
    #endif

    /* call the material integration subroutine */
    // here we do not use the strain/stress vector and tangent matrix directly from the element
    // to prevent memory leaking. We create new ConstitutiveLaw's Parameters to compute the small
    // strain constitutive law.
    Vector NewStrainVector(strain_size), NewStressVector(strain_size);
    Matrix NewTangent(strain_size, strain_size);

    SD_MathUtils<double>::StrainTensorToVector(current_rotated_strain_tensor, NewStrainVector);

    ConstitutiveLaw::Parameters const_params;
    const_params.SetStrainVector(NewStrainVector);
    const_params.SetStressVector(NewStressVector);
    const_params.SetConstitutiveMatrix(NewTangent);
    const_params.SetProcessInfo(rProcessInfo);
    const_params.SetMaterialProperties(rProperties);
    const_params.SetElementGeometry(rGeometry);
    ConstitutiveLaw::StressMeasure stress_measure = ConstitutiveLaw::StressMeasure_Cauchy;

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(const_params);

    /* obtain the stress */
    Matrix current_rotated_stress_tensor(3, 3);
    SD_MathUtils<double>::StressVectorToTensor(NewStressVector, current_rotated_stress_tensor);

    /* transform the stress */
    Matrix current_pk2_stress_tensor(3, 3);
    double J = detU;
    noalias(current_pk2_stress_tensor) = J * prod(right_stretch_tensor_inversed,
            Matrix(prod(current_rotated_stress_tensor, right_stretch_tensor_inversed)));

    Vector& StressVector = rValues.GetStressVector();
    SD_MathUtils<double>::StressTensorToVector(current_pk2_stress_tensor, StressVector);

    /* obtain and transform the tangent */
    if (rValues.IsSetConstitutiveMatrix())
    {
        Fourth_Order_Tensor D;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(D);

        // obtain the tangent from the small strain constitutive law. It must be from
        // a fourth order tensor to not missing the out-of-plane component in plane strain
        // analysis
        if (strain_size != 6)
        {
            if (mpConstitutiveLaw->Has(THREED_ALGORITHMIC_TANGENT))
            {
                Matrix Dmat(6, 6);
                mpConstitutiveLaw->GetValue(THREED_ALGORITHMIC_TANGENT, Dmat);
                SD_MathUtils<double>::MatrixToTensor(Dmat, D);
            }
            else
            {
                // TODO sometimes the error line below will cause non-convergence in 2D analysis.
                // The reason is currently unknown. If that occurs repeatedly, this indicates some
                // problem with memory leaking.
                KRATOS_ERROR << "THREED_ALGORITHMIC_TANGENT is required from constitutive_law for non-3D problem";
            }
        }
        else
        {
            SD_MathUtils<double>::MatrixToTensor(NewTangent, D);
        }

        Fourth_Order_Tensor DSDE;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(DSDE);

        const Matrix& Ui = right_stretch_tensor_inversed;
        for (unsigned int i = 0; i < 3; ++i)
        {
            for (unsigned int j = 0; j < 3; ++j)
            {
                for (unsigned int k = 0; k < 3; ++k)
                {
                    for (unsigned int l = 0; l < 3; ++l)
                    {
                        for (unsigned int m = 0; m < 3; ++m)
                        {
                            for (unsigned int n = 0; n < 3; ++n)
                            {
                                for (unsigned int p = 0; p < 3; ++p)
                                {
                                    for (unsigned int q = 0; q < 3; ++q)
                                    {
                                        DSDE[i][j](k, l) += J * D[m][n](p, q)
                                            * Ui(i, m) * Ui(j, n) * Ui(k, p) * Ui(l, q);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        SD_MathUtils<double>::TensorToMatrix(DSDE, AlgorithmicTangent);

        #ifdef DEBUG_CONSTITUTIVE_LAW
        if (ElemId == 1)
        {
            KRATOS_WATCH(AlgorithmicTangent)
        }
        #endif
    }
}

//**********************************************************************
int TotalLagrangianBridgingConstitutiveLaw::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    return mpConstitutiveLaw->Check( props, geom, CurrentProcessInfo );
}

} // Namespace Kratos

#undef DEBUG_CONSTITUTIVE_LAW
