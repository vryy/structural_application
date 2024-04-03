/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 31 Oct 2023 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

#if !defined(KRATOS_HYPERELASTIC_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED )
#define  KRATOS_HYPERELASTIC_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "constitutive_laws/multiplicative_finite_strain_bridging_constitutive_law.h"


namespace Kratos
{

/**
 * This is the wrapper that allows the hyperelastic constitutive law, i.e. using GL strain and PK2 stress
 * to work with the FiniteStrain element. Basically it transforms the hyperelastic tangent to the compatible
 * form of FiniteStrain.
 * Reference:
 *  +   tr20_finite_strain.pdf
 */
class HyperelasticFiniteStrainBridgingConstitutiveLaw : public FiniteStrainBridgingConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef FiniteStrainBridgingConstitutiveLaw BaseType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::Fourth_Order_Tensor Fourth_Order_Tensor;

    /**
     * Counted pointer of HyperelasticFiniteStrainBridgingConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(HyperelasticFiniteStrainBridgingConstitutiveLaw);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    HyperelasticFiniteStrainBridgingConstitutiveLaw() {}

    /**
     * Constructor with nested constitutive law.
     */
    HyperelasticFiniteStrainBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : BaseType(pConstitutiveLaw)
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new HyperelasticFiniteStrainBridgingConstitutiveLaw(mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~HyperelasticFiniteStrainBridgingConstitutiveLaw()
    {}

    /**
     * Operators
     */

    /**
     * Operations
     */

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
        return "HyperelasticFiniteStrainBridgingConstitutiveLaw";
    }

protected:
    /**
     * Member Variables
     */

    /// Integrate the Cauchy stress
    void StressIntegration(const Parameters& rValues, const Matrix& F, Matrix& stress_tensor) const override;

    /// Compute the consistent tangent (A tensor)
    void ComputeTangent(Fourth_Order_Tensor& A) const override;

private:

    ///@}
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

    /**
     * Static Member Variables
     */

    /**
     * Member Variables
     */

    /**
     * Un accessible methods
     */

    /// Compute an equivalent strain for integrating with the small strain constitutive law
    void ComputeStrain(Vector& StrainVector, const Matrix& F) const;

    /// Compute Cauchy stress from PK2 stress
    void ComputeStress(Matrix& stress_tensor, const Matrix& PK2_stress) const;
}; // Class HyperelasticFiniteStrainBridgingConstitutiveLaw

/**
 * Variant of HyperelasticFiniteStrainBridgingConstitutiveLaw for axisymmetric problem
 */
class HyperelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw : public HyperelasticFiniteStrainBridgingConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef HyperelasticFiniteStrainBridgingConstitutiveLaw BaseType;

    /**
     * Counted pointer of HyperelasticFiniteStrainBridgingConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(HyperelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw);

    /**
     * Default constructor.
     */
    HyperelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw() : BaseType()
    {}

    /**
     * Constructor with nested constitutive law.
     */
    HyperelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : BaseType(pConstitutiveLaw)
    {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new HyperelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw(BaseType::mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

protected:

    unsigned int GetStrainSize(unsigned int dim) const override
    {
        return 4;
    }

}; // HyperelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw

} // namespace Kratos.

#endif // KRATOS_HYPERELASTIC_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED  defined
