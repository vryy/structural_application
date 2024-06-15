/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 5 Oct 2022 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes

// External includes

// Project includes
#include "constitutive_laws/finite_strain_bridging_constitutive_law.h"
#include "structural_application_variables.h"

namespace Kratos
{

//**********************************************************************
FiniteStrainBridgingConstitutiveLaw::FiniteStrainBridgingConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//**********************************************************************
FiniteStrainBridgingConstitutiveLaw::FiniteStrainBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : ConstitutiveLaw(), mpConstitutiveLaw(pConstitutiveLaw)
{
}

//**********************************************************************
FiniteStrainBridgingConstitutiveLaw::~FiniteStrainBridgingConstitutiveLaw()
{
}

//**********************************************************************
bool FiniteStrainBridgingConstitutiveLaw::Has( const Variable<int>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool FiniteStrainBridgingConstitutiveLaw::Has( const Variable<double>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool FiniteStrainBridgingConstitutiveLaw::Has( const Variable<Vector>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
bool FiniteStrainBridgingConstitutiveLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return mpConstitutiveLaw->Has(rThisVariable);
}

//**********************************************************************
int& FiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
double& FiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Vector& FiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Matrix& FiniteStrainBridgingConstitutiveLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    if (rThisVariable == THREED_ALGORITHMIC_TANGENT)
    {
        if (rValue.size1() != 9 || rValue.size2() != 9)
            rValue.resize(9, 9, false);
        this->ComputeTangent( rValue );
        return rValue;
    }
    else if (rThisVariable == CAUCHY_STRESS_TENSOR)
    {
        if (rValue.size1() != 3 || rValue.size2() != 3)
            rValue.resize(3, 3, false);
        noalias(rValue) = m_stress_n1;
        return rValue;
    }

    return mpConstitutiveLaw->GetValue(rThisVariable, rValue);
}

//**********************************************************************
Matrix& FiniteStrainBridgingConstitutiveLaw::CalculateValue( Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    if (rThisVariable == THREED_ALGORITHMIC_TANGENT)
    {
        if (rValue.size1() != 9 || rValue.size2() != 9)
            rValue.resize(9, 9, false);
        this->ComputeTangent( rParameterValues, rValue );
        return rValue;
    }

    return mpConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, rValue);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                            const array_1d<double, 3>& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::SetValue( const Variable<ConstitutiveLaw::Pointer>& rThisVariable, ConstitutiveLaw::Pointer rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == CONSTITUTIVE_LAW)
        mpConstitutiveLaw = rValue;
    else
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::InitializeMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
{
    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetStrainSize(dim);

    mpConstitutiveLaw->InitializeMaterial(props, geom, ShapeFunctionsValues);

    m_F_n.resize(3, 3, false);
    m_F_n1.resize(3, 3, false);
    noalias(m_F_n) = IdentityMatrix(3);
    noalias(m_F_n1) = m_F_n;
    m_J_n1 = 1.0;

    m_stress_n1.resize(3, 3, false);
    noalias(m_stress_n1) = ZeroMatrix(3, 3);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::ResetMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    mpConstitutiveLaw->ResetMaterial(props, geom, ShapeFunctionsValues);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->InitializeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeNonLinearIteration(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->FinalizeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);

    noalias(m_F_n) = m_F_n1;
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::UpdateDeformationGradient(Matrix& F_new, double& J_new, const Parameters& rValues) const
{
    const Matrix& F = rValues.GetDeformationGradientF();
    if (F.size1() == 3)
    {
        noalias(F_new) = F;
    }
    else if (F.size1() == 2)
    {
        F_new.clear();
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                F_new(i, j) = F(i, j);
        F_new(2, 2) = 1.0;
    }
    J_new = rValues.GetDeterminantF();
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::ComputeTangent(Matrix& AlgorithmicTangent) const
{
    Fourth_Order_Tensor A;
    this->ComputeTangent(A);
    SD_MathUtils<double>::TensorToUnsymmetricMatrix(A, AlgorithmicTangent);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::ComputeTangent(const Parameters& rValues, Matrix& AlgorithmicTangent) const
{
    Fourth_Order_Tensor A;
    this->ComputeTangent(rValues, A);
    SD_MathUtils<double>::TensorToUnsymmetricMatrix(A, AlgorithmicTangent);
}

//**********************************************************************
void FiniteStrainBridgingConstitutiveLaw::ComputeStressDerivatives(Fourth_Order_Tensor& D,
        const Parameters& rValues, const Matrix& F, double epsilon) const
{
    // compute reference stress
    Matrix stress_tensor(3, 3);
    this->StressIntegration(rValues, F, stress_tensor);

    Matrix newF(3, 3), new_stress_tensor(3, 3), aux(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newF) = F;
            newF(k, l) += epsilon;

            // integrate the new stress
            this->StressIntegration(rValues, newF, new_stress_tensor);

            noalias(aux) = (new_stress_tensor - stress_tensor) / epsilon;

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    D[i][j][k][l] = aux(i, j);
                }
            }
        }
    }
}

//**********************************************************************
int FiniteStrainBridgingConstitutiveLaw::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo ) const
{
    return mpConstitutiveLaw->Check( props, geom, CurrentProcessInfo );
}

} // Namespace Kratos
