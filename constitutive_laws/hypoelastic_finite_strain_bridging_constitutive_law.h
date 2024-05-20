/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 30 Nov 2023 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

#if !defined(KRATOS_HYPOELASTIC_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED )
#define  KRATOS_HYPOELASTIC_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED

// System includes

// External includes

// Project includes
#include "constitutive_laws/finite_strain_bridging_constitutive_law.h"


namespace Kratos
{

/**
 * A wrapper constitutive law employing (modified) Hughes-Winget mid-point integration scheme for finite strain using hypoelastic approach.
 * TStressType is the stress measure used in the integration algorithm. Cauchy and Kirchhoff stress are supported.
 * Reference:
 *  +   tr20_finite_strain.pdf
 * Remarks:
 * +    the underlying constitutive law must be able to return ELASTIC_STRAIN_TENSOR (ee_n), STRAIN (e_n) and CAUCHY_STRESS_TENSOR
 * +    TStressType: 1: Cauchy stress, 2: Kirchhoff stress
 * +    THWSchemeType: 1: Original Hughes-Winget scheme, 2: Modified Hughes-Winget scheme
 */
template<int THWSchemeType = 2, int TStressType = 2>
class HypoelasticFiniteStrainBridgingConstitutiveLaw : public FiniteStrainBridgingConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef FiniteStrainBridgingConstitutiveLaw BaseType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::Fourth_Order_Tensor Fourth_Order_Tensor;

    /**
     * Counted pointer of HypoelasticFiniteStrainBridgingConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(HypoelasticFiniteStrainBridgingConstitutiveLaw);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    HypoelasticFiniteStrainBridgingConstitutiveLaw() {}

    /**
     * Constructor with nested constitutive law.
     */
    HypoelasticFiniteStrainBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : BaseType(pConstitutiveLaw)
    {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new HypoelasticFiniteStrainBridgingConstitutiveLaw(mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~HypoelasticFiniteStrainBridgingConstitutiveLaw()
    {}

    /**
     * Operators
     */

    /**
     * Operations
     */

    void SetValue( const Variable<array_1d<double, 3 > >& rThisVariable, const array_1d<double, 3 > & rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;

    void InitializeMaterial( const Properties& props,
                             const GeometryType& geom,
                             const Vector& ShapeFunctionsValues ) override;

    void FinalizeSolutionStep( const Properties& props,
                               const GeometryType& geom, //this is just to give the array of nodes
                               const Vector& ShapeFunctionsValues,
                               const ProcessInfo& CurrentProcessInfo ) override;

    /**
     * Computes the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters& rValues) override;

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const override
    {
        if (THWSchemeType == 1 && TStressType == 1)
            return "HypoelasticFiniteStrainBridgingConstitutiveLaw<OriginalHW, Cauchy>";
        else if (THWSchemeType == 1 && TStressType == 2)
            return "HypoelasticFiniteStrainBridgingConstitutiveLaw<OriginalHW, Kirchhoff>";
        else if (THWSchemeType == 2 && TStressType == 1)
            return "HypoelasticFiniteStrainBridgingConstitutiveLaw<ModifiedHW, Cauchy>";
        else if (THWSchemeType == 2 && TStressType == 2)
            return "HypoelasticFiniteStrainBridgingConstitutiveLaw<ModifiedHW, Kirchhoff>";
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

protected:
    /**
     * Member Variables
     */

    array_1d<double, 3> mIntPoint;

    /// Compute the term \nabla_{n+1/2} \Delta u
    void ComputeDDuMidpoint( Matrix& DDu_half, const GeometryType& rGeometry,
            const array_1d<double, 3>& rPoint ) const;

    /// Compute strain
    void ComputeStrain(Vector& StrainVector, const Matrix& DDu_half, const double beta = 0.5) const;

    /// Integrate the stress; depending on the TStressType, the output is Cauchy stress (TStressType=1) or
    /// Kirchhoff stress (TStressType=2)
    void StressIntegration(const Parameters& rValues, const Matrix& F, Matrix& stress_tensor) const override;

    /// Compute the consistent tangent (tensor)
    void ComputeTangent(const Parameters& rValues, Fourth_Order_Tensor& A) const override;

    /// Compute the deformation gradient (but minus I)
    virtual void CalculateDu( const unsigned int dim, Matrix& DDu, const Matrix& G_Operator, const Matrix& CurrentDisp ) const;

    /// Compute G operator
    virtual void CalculateG( Matrix& G_Operator, const Vector& N, const Matrix& DN_DX, const GeometryType& rGeometry ) const;

private:

    Matrix m_stress_n; // Cauchy stress
    double m_J_n;
    Matrix mLastDisp;

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    ///@}

    /**
     * Static Member Variables
     */

    /**
     * Member Variables
     */

    /// Compute the numerical tangent (tensor)
    void ComputeNumericalTangent(Fourth_Order_Tensor& A,
        const Parameters& rValues, const Matrix& DDu, double epsilon) const;

    /*****************************************************************
     *                      AUXILIARY FUNCTIONS                      *
     *****************************************************************/

    /// Compute the term \nabla_{n+1} \Delta u
    void ComputeDDu( Matrix& DDu, const GeometryType& rGeometry,
            const array_1d<double, 3>& rPoint ) const;

    /// Compute the term \nabla_{n+1/2} \Delta u from \nabla_{n+1} \Delta u
    void ComputeDDuMidpoint( Matrix& DDu_half, const Matrix& DDu,
            const GeometryType& rGeometry, const array_1d<double, 3>& rPoint ) const;

    /// Compute tensor M
    void ComputeM( Fourth_Order_Tensor& M, const GeometryType& rGeometry, const array_1d<double, 3>& rPoint ) const;
    void ComputeNumericalM( Fourth_Order_Tensor& M, const Matrix& DDu, const GeometryType& rGeometry,
            const array_1d<double, 3>& rPoint, const double epsilon ) const;
    void ComputeNumericalMskew( Fourth_Order_Tensor& M, const Matrix& DDu, const GeometryType& rGeometry,
            const array_1d<double, 3>& rPoint, const double epsilon ) const;

    void ComputeNumericalDADL( Fourth_Order_Tensor& D, const Matrix& DDu, const GeometryType& rGeometry,
            const array_1d<double, 3>& rPoint, const double epsilon ) const;

    void ComputeNumericalDinvADL( Fourth_Order_Tensor& D, const Matrix& DDu, const GeometryType& rGeometry,
            const array_1d<double, 3>& rPoint, const double epsilon ) const;

    void ComputeNumericalDqDeltaDL( Fourth_Order_Tensor& D, const Matrix& DDu, const GeometryType& rGeometry,
            const array_1d<double, 3>& rPoint, const double epsilon ) const;

    void ComputeNumericalDqDeltaTnqDeltaDL( Fourth_Order_Tensor& D, const Matrix& DDu, const Matrix& stress_n, const GeometryType& rGeometry,
            const array_1d<double, 3>& rPoint, const double epsilon ) const;

}; // Class HypoelasticFiniteStrainBridgingConstitutiveLaw

/**
 * Variant of HypoelasticFiniteStrainBridgingConstitutiveLaw for axisymmetric problem
 */
template<int THWSchemeType, int TStressType>
class HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw : public HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>
{
public:
    /**
     * Type Definitions
     */
    typedef HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType> BaseType;
    typedef typename BaseType::GeometryType GeometryType;

    /**
     * Counted pointer of HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw);

    /**
     * Default constructor.
     */
    HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw() : BaseType()
    {}

    /**
     * Constructor with nested constitutive law.
     */
    HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : BaseType(pConstitutiveLaw)
    {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw(BaseType::mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

protected:

    unsigned int GetStrainSize(unsigned int dim) const override
    {
        return 4;
    }

    /// Compute the deformation gradient (but minus I)
    void CalculateDu( const unsigned int dim, Matrix& DDu, const Matrix& G_Operator, const Matrix& CurrentDisp ) const override;

    /// Compute G operator
    void CalculateG( Matrix& G_Operator, const Vector& N, const Matrix& DN_DX, const GeometryType& rGeometry ) const override;

}; // HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw

} // namespace Kratos.

#endif // KRATOS_HYPOELASTIC_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED  defined
