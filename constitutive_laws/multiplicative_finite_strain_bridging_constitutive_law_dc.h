/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 20 Sep 2023 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

#if !defined(KRATOS_MULTIPLICATIVE_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_DC_INCLUDED )
#define  KRATOS_MULTIPLICATIVE_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_DC_INCLUDED

// System includes

// External includes

// Project includes
#include "constitutive_laws/multiplicative_finite_strain_bridging_constitutive_law.h"


namespace Kratos
{

/**
 * A wrapper constitutive law that enables the use of infinitesimal strain constitutive law in the finite strain setting.
 * This is using the formulation for Isotropic elastic/elasto-plastic materials with logarithmic finite strain extension
 * Reference: De Souza Neto, Computational Plasticity, Box 14.3
 * In different to load control, the material update is performed during FinalizeNonLinearIteration.
 */
template<int TStressType>
class MultiplicativeFiniteStrainBridgingConstitutiveLawDC : public MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>
{
public:
    /**
     * Type Definitions
     */
    typedef MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType> BaseType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef SD_MathUtils<double>::Fourth_Order_Tensor Fourth_Order_Tensor;
    typedef typename BaseType::Parameters Parameters;

    /**
     * Counted pointer of MultiplicativeFiniteStrainBridgingConstitutiveLawDC
     */
    KRATOS_CLASS_POINTER_DEFINITION(MultiplicativeFiniteStrainBridgingConstitutiveLawDC);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    MultiplicativeFiniteStrainBridgingConstitutiveLawDC();

    /**
     * Constructor with nested constitutive law.
     */
    MultiplicativeFiniteStrainBridgingConstitutiveLawDC(ConstitutiveLaw::Pointer pConstitutiveLaw);

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new MultiplicativeFiniteStrainBridgingConstitutiveLawDC(BaseType::mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~MultiplicativeFiniteStrainBridgingConstitutiveLawDC();

    /**
     * Operators
     */

    /**
     * Operations
     */

    bool Has( const Variable<Matrix>& rThisVariable ) override;

    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue ) override;

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
    std::string Info() const override
    {
        return "MultiplicativeFiniteStrainBridgingConstitutiveLawDC";
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

protected:
    /**
     * Member Variables
     */

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
}; // Class MultiplicativeFiniteStrainBridgingConstitutiveLawDC

/**
 * Variant of MultiplicativeFiniteStrainBridgingConstitutiveLaw for axisymmetric problem
 */
template<int TStressType>
class MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLawDC : public MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>
{
public:
    /**
     * Type Definitions
     */
    typedef MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType> BaseType;

    /**
     * Counted pointer of MultiplicativeFiniteStrainBridgingConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLawDC);

    /**
     * Default constructor.
     */
    MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLawDC() : BaseType()
    {}

    /**
     * Constructor with nested constitutive law.
     */
    MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLawDC(ConstitutiveLaw::Pointer pConstitutiveLaw)
    : BaseType(pConstitutiveLaw)
    {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLawDC(BaseType::mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

protected:

    unsigned int GetStrainSize(unsigned int dim) const override
    {
        return 4;
    }

}; // MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLawDC

} // namespace Kratos.

#endif // KRATOS_MULTIPLICATIVE_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_DC_INCLUDED  defined
