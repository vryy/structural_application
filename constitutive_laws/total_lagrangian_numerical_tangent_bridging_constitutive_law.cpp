/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 12 Sep 2024 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "constitutive_laws/total_lagrangian_numerical_tangent_bridging_constitutive_law.h"
#include "structural_application_variables.h"

// #define DEBUG_CONSTITUTIVE_LAW
#define DEBUG_ELEMENT_ID 3161
#define DEBUG_POINT_ID 0

namespace Kratos
{

//**********************************************************************
TotalLagrangianNumericalTangentBridgingConstitutiveLaw::TotalLagrangianNumericalTangentBridgingConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//**********************************************************************
TotalLagrangianNumericalTangentBridgingConstitutiveLaw::TotalLagrangianNumericalTangentBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : ConstitutiveLaw(), mpConstitutiveLaw(pConstitutiveLaw)
{
}

//**********************************************************************
TotalLagrangianNumericalTangentBridgingConstitutiveLaw::~TotalLagrangianNumericalTangentBridgingConstitutiveLaw()
{
}

//**********************************************************************
bool TotalLagrangianNumericalTangentBridgingConstitutiveLaw::Has( const Variable<int>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool TotalLagrangianNumericalTangentBridgingConstitutiveLaw::Has( const Variable<double>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool TotalLagrangianNumericalTangentBridgingConstitutiveLaw::Has( const Variable<Vector>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool TotalLagrangianNumericalTangentBridgingConstitutiveLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
int& TotalLagrangianNumericalTangentBridgingConstitutiveLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
double& TotalLagrangianNumericalTangentBridgingConstitutiveLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Vector& TotalLagrangianNumericalTangentBridgingConstitutiveLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Matrix& TotalLagrangianNumericalTangentBridgingConstitutiveLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::SetValue( const Variable<ConstitutiveLaw::Pointer>& rThisVariable, ConstitutiveLaw::Pointer rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == CONSTITUTIVE_LAW)
        mpConstitutiveLaw = rValue;
    else
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->InitializeMaterial(props, geom, ShapeFunctionsValues);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->ResetMaterial(props, geom, ShapeFunctionsValues);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void TotalLagrangianNumericalTangentBridgingConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rValues)
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

    /* compute the strain tensor */
    const Matrix& F = rValues.GetDeformationGradientF();

    // compute the Green-Lagrange strain tensor
    const Matrix eye = IdentityMatrix(3);
    const Matrix E = 0.5 * (prod(trans(F), F) - eye);
    Vector NewStrainVector(strain_size), NewStressVector(strain_size);
    SD_MathUtils<double>::StrainTensorToVector(E, NewStrainVector);

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if(ElemId == DEBUG_ELEMENT_ID && GaussId == DEBUG_POINT_ID)
    {
        KRATOS_WATCH(F)
        KRATOS_WATCH(NewStrainVector)
    }
    #endif

    // compute the stress
    ConstitutiveLaw::Parameters const_params;
    const_params.SetDeformationGradientF(F);
    const_params.SetStrainVector(NewStrainVector);
    const_params.SetStressVector(NewStressVector);
    const_params.SetProcessInfo(rProcessInfo);
    const_params.SetMaterialProperties(rProperties);
    const_params.SetElementGeometry(rGeometry);

    mpConstitutiveLaw->CalculateMaterialResponsePK2(const_params);
    Vector& StressVector = rValues.GetStressVector();
    noalias(StressVector) = NewStressVector;

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if(ElemId == DEBUG_ELEMENT_ID && GaussId == DEBUG_POINT_ID)
    {
        KRATOS_WATCH(NewStressVector)
    }
    #endif

    /* compute the numerical tangent */
    Matrix perturbed_F = F, inversed_F = F;
    double detF;
    MathUtils<double>::InvertMatrix(F, inversed_F, detF);
    Vector PerturbedStrainVector(strain_size), PerturbedStressVector(strain_size);
    const double epsilon = rProperties.Has(PERTURBATION_FACTOR) ? rProperties[PERTURBATION_FACTOR] : 1e-8;

    Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();

    // obtain the mapping for perturbed deformation gradient
    const auto indices = [=]() -> std::vector<std::pair<int, int> > {
        if (strain_size == 3) // plane strain
        {
            return {{0, 0}, {1, 1}, {0, 1}};
        }
        else if ((strain_size == 4)) // axisymmetric
        {
            return {{0, 0}, {1, 1}, {2, 2}, {0, 1}};
        }
        else if ((strain_size == 6)) // 3D
        {
            return {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}};
        }
        else
            KRATOS_ERROR << "Invalid strain size " << strain_size;
    }
    (); // execute the lambda immediately

    // compute the numerical tangent
    for (int i = 0; i < strain_size; ++i)
    {
        const int c = indices[i].first;
        const int d = indices[i].second;

        noalias(perturbed_F) = F;
        for (int j = 0; j < F.size1(); ++j)
        {
            perturbed_F(j, c) += 0.5 * epsilon * inversed_F(d, j);  // here we do the transpose
            perturbed_F(j, d) += 0.5 * epsilon * inversed_F(c, j);  // here we do the transpose
        }
        const Matrix PerturbedStrainTensor = 0.5 * (prod(trans(perturbed_F), perturbed_F) - eye);
        SD_MathUtils<double>::StrainTensorToVector(PerturbedStrainTensor, PerturbedStrainVector);

        // compute the perturbed stress
        const_params.SetDeformationGradientF(perturbed_F);
        const_params.SetStrainVector(PerturbedStrainVector);
        const_params.SetStressVector(PerturbedStressVector);

        mpConstitutiveLaw->CalculateMaterialResponsePK2(const_params);

        for (int j = 0; j < strain_size; ++j)
        {
            AlgorithmicTangent(i, j) = (PerturbedStressVector(j) - NewStressVector(j)) / epsilon;
        }
    }

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if(ElemId == DEBUG_ELEMENT_ID && GaussId == DEBUG_POINT_ID)
    {
        KRATOS_WATCH(AlgorithmicTangent)
    }
    #endif

    // compute the stress again to restore the internal state of the constitutive law
    const_params.SetStrainVector(NewStrainVector);
    const_params.SetStressVector(NewStressVector);
    const_params.SetProcessInfo(rProcessInfo);
    const_params.SetMaterialProperties(rProperties);
    const_params.SetElementGeometry(rGeometry);

    mpConstitutiveLaw->CalculateMaterialResponsePK2(const_params);
}

//**********************************************************************
int TotalLagrangianNumericalTangentBridgingConstitutiveLaw::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    return mpConstitutiveLaw->Check( props, geom, CurrentProcessInfo );
}

} // Namespace Kratos

#undef DEBUG_CONSTITUTIVE_LAW
#undef DEBUG_ELEMENT_ID
#undef DEBUG_POINT_ID
