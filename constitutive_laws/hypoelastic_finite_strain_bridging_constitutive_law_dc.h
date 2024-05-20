/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 22 Apr 2024 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

#if !defined(KRATOS_HYPOELASTIC_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_DC_INCLUDED )
#define  KRATOS_HYPOELASTIC_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_DC_INCLUDED

// System includes

// External includes

// Project includes
#include "constitutive_laws/hypoelastic_finite_strain_bridging_constitutive_law.h"


namespace Kratos
{

/**
 * A wrapper constitutive law that enables the use of infinitesimal strain constitutive law in the finite strain setting.
 * This is using the formulation for Isotropic elastic/elasto-plastic materials with logarithmic finite strain extension
 * Reference: De Souza Neto, Computational Plasticity, Box 14.3
 * In different to load control, the material update is performed during FinalizeNonLinearIteration.
 */
template<int THWSchemeType, int TStressType>
class HypoelasticFiniteStrainBridgingConstitutiveLawDC : public HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>
{
public:
    /**
     * Type Definitions
     */
    typedef HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType> BaseType;
    typedef typename BaseType::BaseType SuperType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::Fourth_Order_Tensor Fourth_Order_Tensor;
    typedef typename BaseType::Parameters Parameters;

    /**
     * Counted pointer of HypoelasticFiniteStrainBridgingConstitutiveLawDC
     */
    KRATOS_CLASS_POINTER_DEFINITION(HypoelasticFiniteStrainBridgingConstitutiveLawDC);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    HypoelasticFiniteStrainBridgingConstitutiveLawDC() {}

    /**
     * Constructor with nested constitutive law.
     */
    HypoelasticFiniteStrainBridgingConstitutiveLawDC(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : BaseType(pConstitutiveLaw)
    {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new HypoelasticFiniteStrainBridgingConstitutiveLawDC(BaseType::mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~HypoelasticFiniteStrainBridgingConstitutiveLawDC()
    {}

    /**
     * Operators
     */

    /**
     * Operations
     */

    bool Has( const Variable<Matrix>& rThisVariable ) override;

    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;

    void FinalizeNonLinearIteration( const Properties& props,
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
        return "HypoelasticFiniteStrainBridgingConstitutiveLawDC";
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;

private:

    ///@}

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
}; // Class HypoelasticFiniteStrainBridgingConstitutiveLawDC

/**
 * Variant of HypoelasticFiniteStrainBridgingConstitutiveLaw for axisymmetric problem
 */
template<int THWSchemeType, int TStressType>
class HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC : public HypoelasticFiniteStrainBridgingConstitutiveLawDC<THWSchemeType, TStressType>
{
public:
    /**
     * Type Definitions
     */
    typedef HypoelasticFiniteStrainBridgingConstitutiveLawDC<THWSchemeType, TStressType> BaseType;
    typedef typename BaseType::GeometryType GeometryType;

    /**
     * Counted pointer of HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC
     */
    KRATOS_CLASS_POINTER_DEFINITION(HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC);

    /**
     * Default constructor.
     */
    HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC() : BaseType()
    {}

    /**
     * Constructor with nested constitutive law.
     */
    HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : BaseType(pConstitutiveLaw)
    {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC(BaseType::mpConstitutiveLaw->Clone()) );
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

}; // HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC

} // namespace Kratos.

#endif // KRATOS_HYPOELASTIC_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_DC_INCLUDED  defined
