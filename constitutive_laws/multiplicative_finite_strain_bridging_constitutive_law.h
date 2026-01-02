/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 5 Oct 2022 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

#if !defined(KRATOS_MULTIPLICATIVE_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED )
#define  KRATOS_MULTIPLICATIVE_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED

// System includes

// External includes

// Project includes
#include "constitutive_laws/finite_strain_bridging_constitutive_law.h"


namespace Kratos
{

/**
 * A wrapper constitutive law that enables the use of infinitesimal strain constitutive law in the finite strain setting.
 * This is using the formulation for Isotropic elastic/elasto-plastic materials with logarithmic finite strain extension
 * TStressType is the stress measure used in the integration algorithm. Right now, Cauchy and Kirchhoff stress are supported.
 * If the underlying constitutive law uses GL strain and PK2 stress, this constitutive can convert the tangent that is compatible
 * with the FiniteStrain element
 * Reference:
 *  +   De Souza Neto, Computational Plasticity, Box 14.3
 *  +   tr20_finite_strain.pdf
 */
template<int TStressType = 2> // 1: Cauchy stress, 2: Kirchhoff stress
class KRATOS_API(STRUCTURAL_APPLICATION) MultiplicativeFiniteStrainBridgingConstitutiveLaw : public FiniteStrainBridgingConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef FiniteStrainBridgingConstitutiveLaw BaseType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::Fourth_Order_Tensor Fourth_Order_Tensor;

    /**
     * Counted pointer of MultiplicativeFiniteStrainBridgingConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(MultiplicativeFiniteStrainBridgingConstitutiveLaw);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    MultiplicativeFiniteStrainBridgingConstitutiveLaw() : BaseType() {}

    /**
     * Constructor with nested constitutive law.
     */
    MultiplicativeFiniteStrainBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : BaseType(pConstitutiveLaw)
    {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new MultiplicativeFiniteStrainBridgingConstitutiveLaw(mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~MultiplicativeFiniteStrainBridgingConstitutiveLaw()
    {}

    /**
     * Operators
     */

    /**
     * Operations
     */

    void InitializeMaterial( const Properties& props,
                             const GeometryType& geom,
                             const Vector& ShapeFunctionsValues ) override;

    void InitializeSolutionStep( const Properties& props,
                                 const GeometryType& geom, //this is just to give the array of nodes
                                 const Vector& ShapeFunctionsValues,
                                 const ProcessInfo& CurrentProcessInfo ) override;

    /**
     * Computes the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters& rValues) override;

    /**
     * Computes the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    // void CalculateMaterialResponseKirchhoff (Parameters& rValues) override;

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const override;

protected:
    /**
     * Member Variables
     */

    Matrix m_Be_trial;

    /// Compute an equivalent strain for integrating with the small strain constitutive law
    void ComputeStrain(Vector& StrainVector, Matrix& Be_trial, const Matrix& F) const;

    /// Compute the consistent tangent (tensor A)
    void ComputeTangent(Fourth_Order_Tensor& A) const override;

    /// Compute necessary terms to derive the consistent tangent for infinitesimal strain
    void ComputeTangentTerms(Fourth_Order_Tensor& D, Fourth_Order_Tensor& L, Fourth_Order_Tensor& B) const;

private:

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

    /**
     * Un accessible methods
     */

    /// Integrate the stress; depending on the TStressType, the output is Cauchy stress (TStressType=1) or
    /// Kirchhoff stress (TStressType=2)
    void StressIntegration(const Parameters& rValues, const Matrix& F, Matrix& stress_tensor, Matrix& Be_trial) const;

    /// Compute the numerical derivatives of Be^trial with respect to F
    void ComputeBetrialDerivatives(Fourth_Order_Tensor& B, const Matrix& F, double epsilon) const;

    /// Compute the numerical derivatives of (1: Cauchy; 2: Kirchhoff) stress with respect to F
    /// Remark: the current implementation does not work with Fbar
    void ComputeStressDerivatives(Fourth_Order_Tensor& D, const Parameters& rValues, const Matrix& F, double epsilon) const;
}; // Class MultiplicativeFiniteStrainBridgingConstitutiveLaw

/**
 * Variant of MultiplicativeFiniteStrainBridgingConstitutiveLaw for axisymmetric problem
 */
template<int TStressType>
class MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw : public MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>
{
public:
    /**
     * Type Definitions
     */
    typedef MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType> BaseType;

    /**
     * Counted pointer of MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw);

    /**
     * Default constructor.
     */
    MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw() : BaseType()
    {}

    /**
     * Constructor with nested constitutive law.
     */
    MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : BaseType(pConstitutiveLaw)
    {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw(BaseType::mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

protected:

    unsigned int GetStrainSize(unsigned int dim) const override
    {
        return 4;
    }

}; // MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw

} // namespace Kratos.

#endif // KRATOS_MULTIPLICATIVE_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED  defined
