/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 30 Nov 2023 $
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
#include "constitutive_laws/hypoelastic_finite_strain_bridging_constitutive_law.h"
#include "structural_application_variables.h"

// #define DEBUG_CONSTITUTIVE_LAW
#define DEBUG_ELEMENT_ID 1
#define DEBUG_POINT_ID 0

#ifdef DEBUG_CONSTITUTIVE_LAW
#include <iomanip>
#endif

namespace Kratos
{

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::SetValue( const Variable<array_1d<double, 3 > >& rThisVariable, const array_1d<double, 3 > & rValue,
        const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == INTEGRATION_POINT_LOCAL)
    {
        noalias(mIntPoint) = rValue;
    }
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::InitializeMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
    BaseType::InitializeMaterial(props, geom, ShapeFunctionsValues);

    m_stress_n.resize(3, 3, false);
    noalias(m_stress_n) = ZeroMatrix(3, 3);
    m_J_n = BaseType::m_J_n1;

    mLastDisp.resize( geom.size(), 3, false );
    for ( unsigned int node = 0; node < geom.size(); ++node )
        noalias( row( mLastDisp, node ) ) = geom[node].GetSolutionStepValue( DISPLACEMENT );
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    BaseType::FinalizeSolutionStep(props, geom, ShapeFunctionsValues, CurrentProcessInfo);
    noalias(m_stress_n) = BaseType::m_stress_n1;
    m_J_n = BaseType::m_J_n1;

    for ( unsigned int node = 0; node < geom.size(); ++node )
        noalias( row( mLastDisp, node ) ) = geom[node].GetSolutionStepValue( DISPLACEMENT );
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeDDu(Matrix& DDu, const GeometryType& rGeometry,
    const array_1d<double, 3>& rPoint ) const
{
    const unsigned int dim = rGeometry.WorkingSpaceDimension();

    // shape function values and local gradients
    Matrix DN_De;
    rGeometry.ShapeFunctionsLocalGradients( DN_De, rPoint );

    Vector N;
    rGeometry.ShapeFunctionsValues( N, rPoint );

    // current and last displacement
    Matrix CurrentDisp( rGeometry.size(), 3 );
    // Matrix LastDisp( rGeometry.size(), 3 );
    for ( unsigned int node = 0; node < rGeometry.size(); ++node )
    {
        // noalias( row( LastDisp, node ) ) = rGeometry[node].GetSolutionStepValue( DISPLACEMENT, 1 );
        noalias( row( CurrentDisp, node ) ) = rGeometry[node].GetSolutionStepValue( DISPLACEMENT );
    }

    // Jacobian at mid-point
    Matrix DeltaDisp( rGeometry.size(), 3 );
    for ( unsigned int node = 0; node < rGeometry.size(); ++node )
    {
        noalias( row( DeltaDisp, node ) ) = rGeometry[node].Coordinates()
                - rGeometry[node].GetInitialPosition() - row( CurrentDisp, node );
    }

    // deformation gradient at point n+1
    Matrix J;
    rGeometry.Jacobian( J, rPoint, DeltaDisp );

    Matrix InvJ(dim, dim);
    double DetJ;
    MathUtils<double>::InvertMatrix( J, InvJ, DetJ );

    Matrix DN_Dx = prod( DN_De, InvJ );

    Matrix Gx;
    CalculateG( Gx, N, DN_Dx, rGeometry );

    CalculateDu( dim, DDu, Gx, CurrentDisp - mLastDisp );
}

template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeDDuMidpoint(Matrix& DDu, const GeometryType& rGeometry,
    const array_1d<double, 3>& rPoint ) const
{
    const unsigned int dim = rGeometry.WorkingSpaceDimension();

    // shape function values and local gradients
    Matrix DN_De;
    rGeometry.ShapeFunctionsLocalGradients( DN_De, rPoint );

    Vector N;
    rGeometry.ShapeFunctionsValues( N, rPoint );

    // current and last displacement
    Matrix CurrentDisp( rGeometry.size(), 3 );
    // Matrix LastDisp( rGeometry.size(), 3 );
    for ( unsigned int node = 0; node < rGeometry.size(); ++node )
    {
        // noalias( row( LastDisp, node ) ) = rGeometry[node].GetSolutionStepValue( DISPLACEMENT, 1 );
        noalias( row( CurrentDisp, node ) ) = rGeometry[node].GetSolutionStepValue( DISPLACEMENT );
    }

    // Jacobian at mid-point
    Matrix DeltaDisp( rGeometry.size(), 3 );
    for ( unsigned int node = 0; node < rGeometry.size(); ++node )
    {
        noalias( row( DeltaDisp, node ) ) = rGeometry[node].Coordinates() - rGeometry[node].GetInitialPosition()
                - 0.5 * (row( mLastDisp, node ) + row( CurrentDisp, node ));
    }

    // gradient operator at mid-point
    Matrix Jhalf;
    rGeometry.Jacobian( Jhalf, rPoint, DeltaDisp );

    Matrix InvJhalf(dim, dim);
    double DetJhalf;
    MathUtils<double>::InvertMatrix( Jhalf, InvJhalf, DetJhalf );

    Matrix DN_Dx_half = prod( DN_De, InvJhalf );

    Matrix Gx_half;
    CalculateG( Gx_half, N, DN_Dx_half, rGeometry );

    // deformation gradient at mid-point
    CalculateDu( dim, DDu, Gx_half, CurrentDisp - mLastDisp );
}

template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeDDuMidpoint(Matrix& DDu_half, const Matrix& DDu,
    const GeometryType& rGeometry, const array_1d<double, 3>& rPoint ) const
{
    Fourth_Order_Tensor M;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(M);
    this->ComputeM(M, rGeometry, rPoint);

    DDu_half.clear();
    SD_MathUtils<double>::ContractFourthOrderTensor(1.0, M, DDu, DDu_half);
}

template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeM(Fourth_Order_Tensor& M, const GeometryType& rGeometry,
    const array_1d<double, 3>& rPoint) const
{
    const unsigned int dim = rGeometry.WorkingSpaceDimension();

    // shape function values and local gradients
    Matrix DN_De;
    rGeometry.ShapeFunctionsLocalGradients( DN_De, rPoint );

    Vector N;
    rGeometry.ShapeFunctionsValues( N, rPoint );

    // current and last displacement
    Matrix CurrentPos( rGeometry.size(), 3 ), LastPos( rGeometry.size(), 3 );
    for ( unsigned int node = 0; node < rGeometry.size(); ++node )
    {
        noalias( row( LastPos, node ) ) = rGeometry[node].GetInitialPosition() + rGeometry[node].GetSolutionStepValue( DISPLACEMENT, 1 );
        noalias( row( CurrentPos, node ) ) = rGeometry[node].GetInitialPosition() + rGeometry[node].GetSolutionStepValue( DISPLACEMENT );
    }

    // Jacobian at mid-point
    Matrix DeltaDisp( rGeometry.size(), 3 );
    for ( unsigned int node = 0; node < rGeometry.size(); ++node )
    {
        noalias( row( DeltaDisp, node ) ) = rGeometry[node].Coordinates()
                    - 0.5 * (row( LastPos, node ) + row( CurrentPos, node ));
    }

    // gradient operator at mid-point
    Matrix Jhalf;
    rGeometry.Jacobian( Jhalf, rPoint, DeltaDisp );

    Matrix InvJhalf(dim, dim);
    double DetJhalf;
    MathUtils<double>::InvertMatrix( Jhalf, InvJhalf, DetJhalf );

    Matrix DN_Dx_half = prod( DN_De, InvJhalf );

    Matrix Gx_half;
    CalculateG( Gx_half, N, DN_Dx_half, rGeometry );

    // gradient d_xn_half and d_x(n+1)_half
    Matrix Dx_0_half(3, 3), Dx_1_half(3, 3);
    CalculateDu( dim, Dx_1_half, Gx_half, CurrentPos );
    CalculateDu( dim, Dx_0_half, Gx_half, LastPos );

    // tensor operator M
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    M[i][j][k][l] += Dx_0_half(i, k) * Dx_1_half(l, j);
}

template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::CalculateDu( const unsigned int dim,
        Matrix& DDu, const Matrix& G_Operator, const Matrix& CurrentDisp ) const
{
    SD_MathUtils<double>::CalculateF<true>( dim, DDu, G_Operator, CurrentDisp );
}

template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::CalculateG( Matrix& G_Operator,
        const Vector& N, const Matrix& DN_DX, const GeometryType& rGeometry ) const
{
    const unsigned int dim = DN_DX.size2();
    const unsigned int number_of_nodes = N.size();

    if (G_Operator.size1() != dim*dim || G_Operator.size2() != dim*number_of_nodes)
        G_Operator.resize( dim*dim, dim*number_of_nodes, false );

    SD_MathUtils<double>::CalculateG( dim, G_Operator, DN_DX );
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeStrain(Vector& StrainVector, const Matrix& DDu_half, const double beta) const
{
    // compute rotation matrix
    Matrix DDu_half_skew = 0.5 * (DDu_half - trans(DDu_half));
    Matrix DDu_half_sym = 0.5 * (DDu_half + trans(DDu_half));

    const Matrix eye = IdentityMatrix(3);

    Matrix Aux = eye - beta * DDu_half_skew;
    Matrix Auxi(3, 3);
    double detAux;
    MathUtils<double>::InvertMatrix( Aux, Auxi, detAux );

    Matrix qDelta = prod( Auxi, eye + (1 - beta) * DDu_half_skew);

    // compute new strain
    Matrix elastic_strain_tensor(3, 3), elastic_strain_tensor_trial(3, 3);

    mpConstitutiveLaw->GetValue(ELASTIC_STRAIN_TENSOR, elastic_strain_tensor); // TODO change to ELASTIC_STRAIN_TENSOR_OLD (need to also change the material output)

    noalias(elastic_strain_tensor_trial) = prod(qDelta, Matrix(prod(elastic_strain_tensor, trans(qDelta))));

    if constexpr (THWSchemeType == 1)
    {
        noalias(elastic_strain_tensor_trial) += DDu_half_sym;
    }
    else if constexpr (THWSchemeType == 2)
    {
        Matrix qdelta = 0.5 * prod( qDelta, eye + trans(qDelta) );
        noalias(elastic_strain_tensor_trial) += prod(qdelta, Matrix(prod(DDu_half_sym, trans(qdelta))));
    }

    // create the strain vector as input to the small strain constitutive law
    Vector IncrementalStrainVector(StrainVector.size());
    SD_MathUtils<double>::StrainTensorToVector(elastic_strain_tensor_trial - elastic_strain_tensor, IncrementalStrainVector);
    mpConstitutiveLaw->GetValue(STRAIN, StrainVector);
    noalias(StrainVector) += IncrementalStrainVector;
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::StressIntegration(const Parameters& rValues,
        const Matrix& DDu, Matrix& stress_tensor) const
{
    const auto strain_measure = mpConstitutiveLaw->GetStrainMeasure();
    if (strain_measure != ConstitutiveLaw::StrainMeasure_Infinitesimal)
        KRATOS_ERROR << "The strain measure os the constitutive law " << mpConstitutiveLaw->Info()
                     << " is not StrainMeasure_Infinitesimal";

    const GeometryType& geom = rValues.GetElementGeometry();
    const unsigned int dim = geom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetStrainSize(dim);

    // compute the equivalent (small) strain vector
    Vector StrainVector(strain_size);
    this->ComputeStrain(StrainVector, DDu, 0.5); // TODO parameterize beta (currently 0.5)
    #ifdef DEBUG_CONSTITUTIVE_LAW
    int ElemId, PointId;
    mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
    mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, PointId);
    // if (ElemId == DEBUG_ELEMENT_ID && PointId == DEBUG_POINT_ID)
    // {
    //     KRATOS_WATCH(DDu)
    //     KRATOS_WATCH(StrainVector)
    // }
    #endif

    // integrate the constitutive law, obtaining stress
    Vector StressVector(strain_size);
    // Matrix Dmat(strain_size, strain_size);
    //
    ConstitutiveLaw::Parameters const_params;
    const_params.SetStrainVector(StrainVector);
    // const_params.SetDeformationGradientF(DDu);
    const_params.SetStressVector(StressVector);
    // const_params.SetConstitutiveMatrix(Dmat);
    const_params.SetProcessInfo(rValues.GetProcessInfo());
    const_params.SetMaterialProperties(rValues.GetMaterialProperties());
    const_params.SetElementGeometry(geom);

    // for small strain we always integrate the Kirchhoff/Cauchy stress using Cauchy stress integration routine
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(const_params);

    // obtain the stress (1:Cauchy stress; 2: Kirchhoff stress)
    if (!mpConstitutiveLaw->Has(CAUCHY_STRESS_TENSOR))
        KRATOS_ERROR << "Constitutive law is not able to return CAUCHY_STRESS_TENSOR";
    mpConstitutiveLaw->GetValue(CAUCHY_STRESS_TENSOR, stress_tensor);
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    // integrate the Cauchy/Kirchhoff stress
    const GeometryType& geom = rValues.GetElementGeometry();

    const Matrix& F = rValues.GetDeformationGradientF();
    this->UpdateDeformationGradient(m_F_n1, m_J_n1, rValues);

    Matrix DDu_half(3, 3); // TODO moving this calculation to element
    this->ComputeDDuMidpoint( DDu_half, geom, mIntPoint );
    this->StressIntegration(rValues, DDu_half, m_stress_n1);

    if constexpr (TStressType == 2)
    {
        // transform to Cauchy stress
        m_stress_n1 /= m_J_n1;
    }

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
        BaseType::ComputeTangent(rValues, AlgorithmicTangent);

        #ifdef DEBUG_CONSTITUTIVE_LAW
        int ElemId, PointId;
        mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
        mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, PointId);
        if (ElemId == DEBUG_ELEMENT_ID && PointId == DEBUG_POINT_ID)
        {
            KRATOS_WATCH(AlgorithmicTangent)
            KRATOS_WATCH("------------------------")
        }
        #endif
    }
}

template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeTangent(const Parameters& rValues, Fourth_Order_Tensor& AA) const
{
    const GeometryType& rGeometry = rValues.GetElementGeometry();

    const Matrix eye = IdentityMatrix(3);

    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(AA);

    // compute rotation matrix
    Matrix DDu_half(3, 3);
    this->ComputeDDuMidpoint( DDu_half, rGeometry, mIntPoint );

    Matrix DDu_half_skew = 0.5 * (DDu_half - trans(DDu_half));
    Matrix DDu_half_sym = 0.5 * (DDu_half + trans(DDu_half));

    const double beta = 0.5; // TODO parameterize beta
    Matrix Aux = eye - beta * DDu_half_skew;
    Matrix Auxi(3, 3);
    double detAux;
    MathUtils<double>::InvertMatrix( Aux, Auxi, detAux );

    Matrix qDelta = prod( Auxi, eye + (1 - beta) * DDu_half_skew);

    // obtain stress from previous step
    Matrix stress_n(3, 3);
    if constexpr (TStressType == 1)       // Cauchy stress
        noalias(stress_n) = m_stress_n;
    else if constexpr (TStressType == 2)  // Kirchhoff stress
        noalias(stress_n) = m_J_n*m_stress_n;

    // compute tensor A, B
    Fourth_Order_Tensor A, B;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(A);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(B);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                for (int l = 0; l < 3; ++l)
                {
                    for (int m = 0; m < 3; ++m)
                        A[i][j][k][l] += eye(i, k)*stress_n(l, m)*qDelta(j, m) + qDelta(i, m)*stress_n(m, l)*eye(j, k);

                    B[i][j][k][l] = 0.5*Auxi(i, k)*(eye(l, j) + qDelta(l, j));
                }
            }
        }
    }

    // compute Mskew, Msym
    Fourth_Order_Tensor M, Mskew, Msym;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(M);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Mskew);
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Msym);

    this->ComputeM(M, rGeometry, mIntPoint);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                for (int l = 0; l < 3; ++l)
                {
                    Msym[i][j][k][l] = 0.5 * (M[i][j][k][l] + M[j][i][k][l]);
                    Mskew[i][j][k][l] = 0.5 * (M[i][j][k][l] - M[j][i][k][l]);
                }
            }
        }
    }

    #ifdef DEBUG_CONSTITUTIVE_LAW
    Matrix DDu(3, 3);
    this->ComputeDDu(DDu, rGeometry, mIntPoint);

    int ElemId, PointId;
    mpConstitutiveLaw->GetValue(PARENT_ELEMENT_ID, ElemId);
    mpConstitutiveLaw->GetValue(INTEGRATION_POINT_INDEX, PointId);
    if (ElemId == DEBUG_ELEMENT_ID && PointId == DEBUG_POINT_ID)
    {
        Fourth_Order_Tensor NumericalM;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(NumericalM);

        this->ComputeNumericalM(NumericalM, DDu, rGeometry, mIntPoint, 1e-7);

        KRATOS_WATCH(M)
        KRATOS_WATCH(NumericalM)

        SD_MathUtils<double>::AddFourthOrderTensor(-1.0, M, NumericalM);
        const double diffM = SD_MathUtils<double>::normTensor(NumericalM);
        KRATOS_WATCH(diffM)

        ///////////

        Fourth_Order_Tensor NumericalMskew;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(NumericalMskew);

        this->ComputeNumericalMskew(NumericalMskew, DDu, rGeometry, mIntPoint, 1e-7);

        KRATOS_WATCH(Mskew)
        KRATOS_WATCH(NumericalMskew)

        SD_MathUtils<double>::AddFourthOrderTensor(-1.0, Mskew, NumericalMskew);
        const double diffMskew = SD_MathUtils<double>::normTensor(NumericalMskew);
        KRATOS_WATCH(diffMskew)

        ///////////

        Fourth_Order_Tensor DADL, Numerical_DADL;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(DADL);
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Numerical_DADL);

        this->ComputeNumericalDADL(Numerical_DADL, DDu, rGeometry, mIntPoint, 1e-8);

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    for (int l = 0; l < 3; ++l)
                        DADL[i][j][k][l] = -0.5*Mskew[i][j][k][l];

        KRATOS_WATCH(DADL)
        KRATOS_WATCH(Numerical_DADL)
        SD_MathUtils<double>::AddFourthOrderTensor(-1.0, DADL, Numerical_DADL);
        const double diff_DADL = SD_MathUtils<double>::normTensor(Numerical_DADL);
        KRATOS_WATCH(diff_DADL)

        ///////////

        Fourth_Order_Tensor DinvADL, Numerical_DinvADL;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(DinvADL);
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Numerical_DinvADL);

        this->ComputeNumericalDinvADL(Numerical_DinvADL, DDu, rGeometry, mIntPoint, 1e-7);

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    for (int l = 0; l < 3; ++l)
                        for (int m = 0; m < 3; ++m)
                            for (int n = 0; n < 3; ++n)
                                DinvADL[i][j][k][l] += 0.5*Auxi(i, m)*Mskew[m][n][k][l]*(Auxi(n, j));

        KRATOS_WATCH(DinvADL)
        KRATOS_WATCH(Numerical_DinvADL)
        SD_MathUtils<double>::AddFourthOrderTensor(-1.0, DinvADL, Numerical_DinvADL);
        const double diff_DinvADL = SD_MathUtils<double>::normTensor(Numerical_DinvADL);
        KRATOS_WATCH(diff_DinvADL)

        ///////////

        Fourth_Order_Tensor DqDeltaDL, Numerical_DqDeltaDL;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(DqDeltaDL);
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Numerical_DqDeltaDL);

        this->ComputeNumericalDqDeltaDL(Numerical_DqDeltaDL, DDu, rGeometry, mIntPoint, 1e-9);
        SD_MathUtils<double>::ProductFourthOrderTensor(1.0, B, Mskew, DqDeltaDL);

        KRATOS_WATCH(DqDeltaDL)
        KRATOS_WATCH(Numerical_DqDeltaDL)

        SD_MathUtils<double>::AddFourthOrderTensor(-1.0, DqDeltaDL, Numerical_DqDeltaDL);
        const double diff_DqDeltaDL = SD_MathUtils<double>::normTensor(Numerical_DqDeltaDL);
        KRATOS_WATCH(diff_DqDeltaDL)

        ///////////

        Fourth_Order_Tensor DqDeltaTnqDeltaDL, Numerical_DqDeltaTnqDeltaDL;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(DqDeltaTnqDeltaDL);
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Numerical_DqDeltaTnqDeltaDL);

        this->ComputeNumericalDqDeltaTnqDeltaDL(Numerical_DqDeltaTnqDeltaDL, DDu, stress_n, rGeometry, mIntPoint, 1e-9);
        SD_MathUtils<double>::ProductFourthOrderTensor(1.0, A, B, DqDeltaTnqDeltaDL);
        SD_MathUtils<double>::ProductFourthOrderTensor(1.0, DqDeltaTnqDeltaDL, Mskew); // A:B:Mskew

        KRATOS_WATCH(stress_n)
        KRATOS_WATCH(DqDeltaTnqDeltaDL)
        KRATOS_WATCH(Numerical_DqDeltaTnqDeltaDL)

        SD_MathUtils<double>::AddFourthOrderTensor(-1.0, DqDeltaTnqDeltaDL, Numerical_DqDeltaTnqDeltaDL);
        const double diff_DqDeltaTnqDeltaDL = SD_MathUtils<double>::normTensor(Numerical_DqDeltaTnqDeltaDL);
        KRATOS_WATCH(diff_DqDeltaTnqDeltaDL)
    }
    #endif

    // obtain the tangent from the small strain constitutive law. It must be from
    // a fourth order tensor to not missing the out-of-plane component in plane strain
    // analysis
     if (!mpConstitutiveLaw->Has(THREED_ALGORITHMIC_TANGENT))
        KRATOS_ERROR << "Constitutive law is not able to return THREED_ALGORITHMIC_TANGENT";
    Matrix Dmat(6, 6);
    mpConstitutiveLaw->GetValue(THREED_ALGORITHMIC_TANGENT, Dmat);

    Fourth_Order_Tensor Dm;
    SD_MathUtils<double>::CalculateFourthOrderZeroTensor(Dm);
    SD_MathUtils<double>::MatrixToTensor(Dmat, Dm);

    // compute tensor AA
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0, A, B, AA);
    SD_MathUtils<double>::ProductFourthOrderTensor(1.0, AA, Mskew); // AA is now A:B:Mskew

    if constexpr (THWSchemeType == 1)       // original Hughes-Winget scheme
    {
        SD_MathUtils<double>::ProductFourthOrderTensor(1.0, Dm, Msym, AA); // AA is now A:B:Mskew + Ce:Msym
    }
    else if constexpr (THWSchemeType == 2)  // modified Hughes-Winget scheme
    {
        Fourth_Order_Tensor C, D, E;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(C);
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(D);
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(E);

        Matrix qdelta = 0.5 * prod( qDelta, eye + trans(qDelta) );

        // compute tensor C, D, E
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    for (int l = 0; l < 3; ++l)
                    {
                        C[i][j][k][l] = 0.5*(eye(i, k)*eye(j, l) + eye(i, k)*qDelta(j, l) + qDelta(i, l)*eye(j, k));

                        for (int m = 0; m < 3; ++m)
                            D[i][j][k][l] += eye(i, k)*DDu_half_sym(l, m)*qdelta(j, m) + qdelta(i, m)*DDu_half_sym(m, l)*eye(j, k);

                        E[i][j][k][l] = qdelta(i, k)*qdelta(j, l);
                    }
                }
            }
        }

        SD_MathUtils<double>::ProductFourthOrderTensor(1.0, D, C);
        SD_MathUtils<double>::ProductFourthOrderTensor(1.0, D, B);
        SD_MathUtils<double>::ProductFourthOrderTensor(1.0, D, Mskew);
        SD_MathUtils<double>::ProductFourthOrderTensor(1.0, E, Msym, D); // D is now (D:C:B:Mskew + E:Msym)

        SD_MathUtils<double>::ProductFourthOrderTensor(1.0, Dm, D, AA); // AA is now A:B:Mskew + Ce:(D:C:B:Mskew + E:Msym)
    }

    ///////////////////////////////////////

    #ifdef DEBUG_CONSTITUTIVE_LAW
    if (ElemId == DEBUG_ELEMENT_ID && PointId == DEBUG_POINT_ID)
    {
        // debugging
        Fourth_Order_Tensor num_AA;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(num_AA);
        Matrix DDu(3, 3);
        this->ComputeDDu(DDu, rGeometry, mIntPoint);
        const double epsilon = 1e-8;
        this->ComputeNumericalTangent(num_AA, rValues, DDu, epsilon);
        KRATOS_WATCH(AA)
        KRATOS_WATCH(num_AA)
        SD_MathUtils<double>::AddFourthOrderTensor(-1.0, AA, num_AA);
        const double diffAA = SD_MathUtils<double>::normTensor(num_AA);
        KRATOS_WATCH(diffAA)
        KRATOS_WATCH("------------------------------------")
        // end debugging
    }
    #endif

    // transform the tangent to compatible element tangent

    if constexpr (TStressType == 1)
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    for (int l = 0; l < 3; ++l)
                        AA[i][j][k][l] += (m_stress_n1(i, j) * eye(k, l) - m_stress_n1(i, l) * eye(j, k));
    }
    else if constexpr (TStressType == 2)
    {
        // const double J = MathUtils<double>::Det(m_F_n1);
        const double J = m_J_n1; // TODO check this if it is OK to use with Fbar
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    for (int l = 0; l < 3; ++l)
                        AA[i][j][k][l] = AA[i][j][k][l] / J - m_stress_n1(i, l) * eye(j, k);
    }
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeNumericalTangent(Fourth_Order_Tensor& D,
        const Parameters& rValues, const Matrix& DDu, double epsilon) const
{
    const GeometryType& rGeometry = rValues.GetElementGeometry();

    // compute reference stress
    Matrix stress_tensor(3, 3), DDu_half(3, 3);
    this->ComputeDDuMidpoint(DDu_half, DDu, rGeometry, mIntPoint);
    this->StressIntegration(rValues, DDu_half, stress_tensor);

    Matrix newDDu(3, 3), new_stress_tensor(3, 3), aux(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newDDu) = DDu;
            newDDu(k, l) += epsilon;

            // integrate the new stress
            // for some reason the out-of-plane components of the numerical tangent is always double the correct values
            // it is probably because the stress tensor is symmetric. Do not know how to make it correctly by 4/2024 -> TODO
            this->ComputeDDuMidpoint(DDu_half, newDDu, rGeometry, mIntPoint);
            this->StressIntegration(rValues, DDu_half, new_stress_tensor);

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

    // compute stress again to restore the correct constitutive law state
    this->ComputeDDuMidpoint(DDu_half, DDu, rGeometry, mIntPoint);
    this->StressIntegration(rValues, DDu_half, stress_tensor);
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeNumericalM( Fourth_Order_Tensor& M, const Matrix& DDu,
        const GeometryType& rGeometry, const array_1d<double, 3>& rPoint, const double epsilon ) const
{
    // compute reference state
    Matrix DDu_half(3, 3);
    this->ComputeDDuMidpoint(DDu_half, DDu, rGeometry, rPoint);

    Matrix newDDu(3, 3), newDDu_half(3, 3), aux(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newDDu) = DDu;
            newDDu(k, l) += epsilon;

            // integrate the new state
            this->ComputeDDuMidpoint(newDDu_half, newDDu, rGeometry, rPoint);

            noalias(aux) = (newDDu_half - DDu_half) / epsilon;

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    M[i][j][k][l] = aux(i, j);
                }
            }
        }
    }
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeNumericalMskew( Fourth_Order_Tensor& M, const Matrix& DDu,
        const GeometryType& rGeometry, const array_1d<double, 3>& rPoint, const double epsilon ) const
{
    // compute reference state
    Matrix DDu_half(3, 3);
    this->ComputeDDuMidpoint(DDu_half, DDu, rGeometry, rPoint);
    Matrix DDu_half_skew = 0.5 * (DDu_half - trans(DDu_half));

    Matrix newDDu(3, 3), newDDu_half_skew(3, 3), aux(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newDDu) = DDu;
            newDDu(k, l) += epsilon;

            // integrate the new state
            this->ComputeDDuMidpoint(DDu_half, newDDu, rGeometry, rPoint);

            noalias(newDDu_half_skew) = 0.5 * (DDu_half - trans(DDu_half));

            noalias(aux) = (newDDu_half_skew - DDu_half_skew) / epsilon;

            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j)
                {
                    M[i][j][k][l] = aux(i, j);
                }
            }
        }
    }
}

//**********************************************************************
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeNumericalDqDeltaDL( Fourth_Order_Tensor& D, const Matrix& DDu, const GeometryType& rGeometry,
            const array_1d<double, 3>& rPoint, const double epsilon ) const
{
    // compute reference state
    Matrix DDu_half(3, 3);
    this->ComputeDDuMidpoint( DDu_half, DDu, rGeometry, rPoint );

    Matrix DDu_half_skew = 0.5 * (DDu_half - trans(DDu_half));

    const double beta = 0.5; // TODO parameterize beta
    const Matrix eye = IdentityMatrix(3);
    Matrix Aux = eye - beta * DDu_half_skew;
    Matrix Auxi(3, 3);
    double detAux;
    MathUtils<double>::InvertMatrix( Aux, Auxi, detAux );

    const Matrix qDelta = prod( Auxi, eye + (1 - beta) * DDu_half_skew);

    Matrix newDDu(3, 3), new_qDelta(3, 3), aux(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newDDu) = DDu;
            newDDu(k, l) += epsilon;

            // integrate the new state
            this->ComputeDDuMidpoint( DDu_half, newDDu, rGeometry, rPoint );

            noalias(DDu_half_skew) = 0.5 * (DDu_half - trans(DDu_half));

            noalias(Aux) = eye - beta * DDu_half_skew;
            MathUtils<double>::InvertMatrix( Aux, Auxi, detAux );

            noalias(new_qDelta) = prod( Auxi, eye + (1 - beta) * DDu_half_skew);

            noalias(aux) = (new_qDelta - qDelta) / epsilon;

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
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeNumericalDqDeltaTnqDeltaDL( Fourth_Order_Tensor& D, const Matrix& DDu, const Matrix& stress_n, const GeometryType& rGeometry,
            const array_1d<double, 3>& rPoint, const double epsilon ) const
{
    // compute reference state
    Matrix DDu_half(3, 3);
    this->ComputeDDuMidpoint( DDu_half, DDu, rGeometry, rPoint );

    Matrix DDu_half_skew = 0.5 * (DDu_half - trans(DDu_half));

    const double beta = 0.5; // TODO parameterize beta
    const Matrix eye = IdentityMatrix(3);
    Matrix Aux = eye - beta * DDu_half_skew;
    Matrix Auxi(3, 3);
    double detAux;
    MathUtils<double>::InvertMatrix( Aux, Auxi, detAux );

    Matrix qDelta = prod( Auxi, eye + (1 - beta) * DDu_half_skew);
    Matrix stress_n_rotated(3, 3);
    noalias(stress_n_rotated) = prod(qDelta, Matrix(prod(stress_n, trans(qDelta))));

    Matrix new_stress_n_rotated(3, 3), newDDu(3, 3), aux(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newDDu) = DDu;
            newDDu(k, l) += epsilon;

            // integrate the new state
            this->ComputeDDuMidpoint( DDu_half, newDDu, rGeometry, rPoint );

            noalias(DDu_half_skew) = 0.5 * (DDu_half - trans(DDu_half));

            noalias(Aux) = eye - beta * DDu_half_skew;
            MathUtils<double>::InvertMatrix( Aux, Auxi, detAux );

            noalias(qDelta) = prod( Auxi, eye + (1 - beta) * DDu_half_skew);

            noalias(new_stress_n_rotated) = prod(qDelta, Matrix(prod(stress_n, trans(qDelta))));

            noalias(aux) = (new_stress_n_rotated - stress_n_rotated) / epsilon;

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
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeNumericalDADL( Fourth_Order_Tensor& D, const Matrix& DDu, const GeometryType& rGeometry,
        const array_1d<double, 3>& rPoint, const double epsilon ) const
{
    // compute reference state
    Matrix DDu_half(3, 3);
    this->ComputeDDuMidpoint( DDu_half, DDu, rGeometry, rPoint );

    Matrix DDu_half_skew = 0.5 * (DDu_half - trans(DDu_half));

    const double beta = 0.5; // TODO parameterize beta
    const Matrix eye = IdentityMatrix(3);
    Matrix Aux = eye - beta * DDu_half_skew;

    Matrix newDDu(3, 3), new_Aux(3, 3), aux(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newDDu) = DDu;
            newDDu(k, l) += epsilon;

            // integrate the new state
            this->ComputeDDuMidpoint( DDu_half, newDDu, rGeometry, rPoint );

            noalias(DDu_half_skew) = 0.5 * (DDu_half - trans(DDu_half));

            noalias(new_Aux) = eye - beta * DDu_half_skew;

            noalias(aux) = (new_Aux - Aux) / epsilon;

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
template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>::ComputeNumericalDinvADL( Fourth_Order_Tensor& D,
        const Matrix& DDu, const GeometryType& rGeometry,
        const array_1d<double, 3>& rPoint, const double epsilon ) const
{
    // compute reference state
    Matrix DDu_half(3, 3);
    this->ComputeDDuMidpoint( DDu_half, DDu, rGeometry, rPoint );

    Matrix DDu_half_skew = 0.5 * (DDu_half - trans(DDu_half));

    const double beta = 0.5; // TODO parameterize beta
    const Matrix eye = IdentityMatrix(3);
    Matrix Aux = eye - beta * DDu_half_skew;
    Matrix Auxi(3, 3);
    double detAux;
    MathUtils<double>::InvertMatrix( Aux, Auxi, detAux );

    Matrix newDDu(3, 3), new_Auxi(3, 3), aux(3, 3);
    for (unsigned int k = 0; k < 3; ++k)
    {
        for (unsigned int l = 0; l < 3; ++l)
        {
            noalias(newDDu) = DDu;
            newDDu(k, l) += epsilon;

            // compute the new state
            this->ComputeDDuMidpoint( DDu_half, newDDu, rGeometry, rPoint );

            noalias(DDu_half_skew) = 0.5 * (DDu_half - trans(DDu_half));

            noalias(Aux) = eye - beta * DDu_half_skew;
            MathUtils<double>::InvertMatrix( Aux, new_Auxi, detAux );

            noalias(aux) = (new_Auxi - Auxi) / epsilon;

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
//**********************************************************************

template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw<THWSchemeType, TStressType>::CalculateDu( const unsigned int dim,
        Matrix& DDu, const Matrix& G_Operator, const Matrix& CurrentDisp ) const
{
    SD_MathUtils<double>::CalculateFaxi<true>( DDu, G_Operator, CurrentDisp );
}

template<int THWSchemeType, int TStressType>
void HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw<THWSchemeType, TStressType>::CalculateG( Matrix& G_Operator,
        const Vector& N, const Matrix& DN_DX, const GeometryType& rGeometry ) const
{
    const unsigned int number_of_nodes = N.size();

    if (G_Operator.size1() != 5 || G_Operator.size2() != 2*number_of_nodes)
        G_Operator.resize( 5, 2*number_of_nodes, false );

    SD_MathUtils<double>::CalculateGaxi( G_Operator, rGeometry, N, DN_DX );
}

//**********************************************************************

template class HypoelasticFiniteStrainBridgingConstitutiveLaw<1, 1>;
template class HypoelasticFiniteStrainBridgingConstitutiveLaw<1, 2>;
template class HypoelasticFiniteStrainBridgingConstitutiveLaw<2, 1>;
template class HypoelasticFiniteStrainBridgingConstitutiveLaw<2, 2>;

template class HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw<1, 1>;
template class HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw<1, 2>;
template class HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw<2, 1>;
template class HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw<2, 2>;

} // Namespace Kratos

#undef DEBUG_CONSTITUTIVE_LAW
#undef DEBUG_ELEMENT_ID
#undef DEBUG_POINT_ID
