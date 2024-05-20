/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 31 Oct 2022 $
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
#include "constitutive_laws/hyperelastic_finite_strain_bridging_constitutive_law.h"
#include "structural_application_variables.h"

// #define DEBUG_CONSTITUTIVE_LAW
#define DEBUG_ELEMENT_ID 1

#ifdef DEBUG_CONSTITUTIVE_LAW
#include <iomanip>
#endif

namespace Kratos
{

//**********************************************************************
void HyperelasticFiniteStrainBridgingConstitutiveLaw::StressIntegration(const Parameters& rValues,
        const Matrix& F, Matrix& stress_tensor) const
{
    const auto& geom = rValues.GetElementGeometry();

    const auto strain_measure = mpConstitutiveLaw->GetStrainMeasure();
    if (strain_measure != ConstitutiveLaw::StrainMeasure_GreenLagrange)
        KRATOS_ERROR << "The strain measure is not StrainMeasure_GreenLagrange";

    #ifdef DEBUG_CONSTITUTIVE_LAW
    int ElemId, GaussId;
    mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
    mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, GaussId);
    // std::cout << std::setprecision(16);
    // ElemId = 1;
    // GaussId = 0;
    #endif

    #ifdef DEBUG_CONSTITUTIVE_LAW
    // if (ElemId == DEBUG_ELEMENT_ID && GaussId == 0)
    if (ElemId == DEBUG_ELEMENT_ID)
    {
        KRATOS_WATCH(F)
    }
    #endif

    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetStrainSize(dim);

    // compute the Green-Lagrange strain
    Matrix C;
    Vector StrainVector(strain_size);
    if (dim == 2)
    {
        Matrix F2d(2, 2);
        noalias(F2d) = subrange(F, 0, 2, 0, 2);
        C = prod( trans( F2d ), F2d );
    }
    else if (dim == 3)
    {
        C = prod( trans( F ), F );
    }
    this->ComputeStrain(StrainVector, C);

    // integrate the constitutive law, obtaining stress
    Vector StressVector(strain_size);
    Matrix Dmat(strain_size, strain_size);
    //
    ConstitutiveLaw::Parameters const_params;
    const_params.SetStrainVector(StrainVector);
    const_params.SetDeformationGradientF(F);
    const_params.SetStressVector(StressVector);
    const_params.SetConstitutiveMatrix(Dmat);
    const_params.SetProcessInfo(rValues.GetProcessInfo());
    const_params.SetMaterialProperties(rValues.GetMaterialProperties());
    const_params.SetElementGeometry(rValues.GetElementGeometry());

    // integrate the PK2 stress
    mpConstitutiveLaw->CalculateMaterialResponsePK2(const_params);

    // obtain PK2 stress
    Matrix pk2_stress(3, 3);
    mpConstitutiveLaw->GetValue(PK2_STRESS_TENSOR, pk2_stress);

    // transform to Cauchy stress
    this->ComputeStress(stress_tensor, pk2_stress);
}

//**********************************************************************
void HyperelasticFiniteStrainBridgingConstitutiveLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    // integrate the Cauchy stress
    const Matrix& F = rValues.GetDeformationGradientF();
    this->UpdateDeformationGradient(m_F_n1, m_J_n1, rValues);
    this->StressIntegration(rValues, m_F_n1, m_stress_n1);

    // export the stress
    if (rValues.IsSetStressVector())
    {
        Vector& CauchyStressVector = rValues.GetStressVector();
        SD_MathUtils<double>::StressTensorToVector(m_stress_n1, CauchyStressVector);
    }

    // export the tangent
    if (rValues.IsSetConstitutiveMatrix())
    {
        Matrix& AlgorithmicTangent = rValues.GetConstitutiveMatrix();
        BaseType::ComputeTangent(AlgorithmicTangent);
    }
}

//**********************************************************************
void HyperelasticFiniteStrainBridgingConstitutiveLaw::ComputeTangent(Fourth_Order_Tensor& A) const
{
    if (!mpConstitutiveLaw->Has(THREED_ALGORITHMIC_TANGENT))
        KRATOS_ERROR << "Constitutive law is not able to return THREED_ALGORITHMIC_TANGENT";

    // obtain the tangent from the small strain constitutive law. It must be from
    // a fourth order tensor to not missing the out-of-plane component in plane strain
    // analysis
    Matrix Dmat(6, 6);
    mpConstitutiveLaw->GetValue(THREED_ALGORITHMIC_TANGENT, Dmat);

    Fourth_Order_Tensor D;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(D);
    SD_MathUtils<double>::MatrixToTensor(Dmat, D);

    // compute the tensor a
    const Matrix& F = m_F_n1;
    double J = m_J_n1;
    // KRATOS_WATCH(J)

    const Matrix& sigma = m_stress_n1;

    double aux;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    const Matrix eye = IdentityMatrix(3);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                for (int l = 0; l < 3; ++l)
                {
                    aux = eye(i, k) * sigma(l, j);

                    for (int m = 0; m < 3; ++m)
                        for (int u = 0; u < 3; ++u)
                            for (int p = 0; p < 3; ++p)
                                for (int q = 0; q < 3; ++q)
                                    aux += 0.5/J * F(i, m) * F(j, u) * D[m][u][p][q]
                                            * (F(l, p)*F(k, q) + F(l, q)*F(k, p));

                    A[i][j][k][l] = aux;
                }
            }
        }
    }
}

//**********************************************************************
void HyperelasticFiniteStrainBridgingConstitutiveLaw::ComputeStrain( Vector& StrainVector, const Matrix& C ) const
{
    KRATOS_TRY

    const unsigned int dimension = C.size1();

    if ( dimension == 2 )
    {
        if ( StrainVector.size() != 3 ) StrainVector.resize( 3, false );

        StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        StrainVector[2] = C( 0, 1 );
    }

    if ( dimension == 3 )
    {
        if ( StrainVector.size() != 6 ) StrainVector.resize( 6, false );

        StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        StrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

        StrainVector[3] = C( 0, 1 ); // xy

        StrainVector[4] = C( 1, 2 ); // yz

        StrainVector[5] = C( 0, 2 ); // xz
    }

    KRATOS_CATCH( "" )
}

//**********************************************************************
void HyperelasticFiniteStrainBridgingConstitutiveLaw::ComputeStress(Matrix& stress_tensor, const Matrix& PK2_stress) const
{
    noalias(stress_tensor) = 1.0/m_J_n1 * prod(m_F_n1, Matrix(prod(PK2_stress, trans(m_F_n1))));
}

//**********************************************************************
// void HyperelasticFiniteStrainBridgingConstitutiveLaw::ComputeStressDerivatives(Fourth_Order_Tensor& D,
//         const Parameters& rValues, const Matrix& F, double epsilon) const
// {
//     // integrate the stress
//     Matrix stress_tensor(3, 3);
//     this->StressIntegration(rValues, F, stress_tensor);

//     Matrix newF(3, 3), new_stress_tensor(3, 3), aux(3, 3);
//     for (unsigned int k = 0; k < 3; ++k)
//     {
//         for (unsigned int l = 0; l < 3; ++l)
//         {
//             noalias(newF) = F;
//             newF(k, l) += epsilon;

//             // integrate the new stress
//             this->StressIntegration(rValues, newF, new_stress_tensor);

//             noalias(aux) = (new_stress_tensor - stress_tensor) / epsilon;

//             for (unsigned int i = 0; i < 3; ++i)
//             {
//                 for (unsigned int j = 0; j < 3; ++j)
//                 {
//                     D[i][j][k][l] = aux(i, j);
//                 }
//             }
//         }
//     }
// }

} // Namespace Kratos

#undef DEBUG_CONSTITUTIVE_LAW
