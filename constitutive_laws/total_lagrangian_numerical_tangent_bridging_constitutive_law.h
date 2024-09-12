/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 12 Sep 2024 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

#if !defined(KRATOS_TOTAL_LAGRANGIAN_NUMERICAL_TANGENT_BRIDGING_CONSTITUTIVE_LAW_INCLUDED )
#define  KRATOS_TOTAL_LAGRANGIAN_NUMERICAL_TANGENT_BRIDGING_CONSTITUTIVE_LAW_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

/**
 * A wrapper constitutive law that approximates the tangent. This works with finite strain
 * constitutive law. It does not need to provide the tangent. This will be computed approximately
 * using a finite difference scheme.
 * Remarks:
 *  +   only work with TotalLagrangian element.
 *  +   only work with Constitutive law based on PK2.
 * +    it can be chained up with other bridging law to provide better tangent, e.g. TotalLagrangianNumericalTangentBridgingConstitutiveLaw(TotalLagrangianBridgingConstitutiveLaw(VonMises3dImplicit()))
 * Reference:
 *  + Miehe, Numerical computation of algorithmic (consistent) tangent moduli in large-strain computational inelasticity, CMAME 1995.
 */
class TotalLagrangianNumericalTangentBridgingConstitutiveLaw : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;
    typedef BaseType::SizeType SizeType;

    /**
     * Counted pointer of TotalLagrangianNumericalTangentBridgingConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(TotalLagrangianNumericalTangentBridgingConstitutiveLaw);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    TotalLagrangianNumericalTangentBridgingConstitutiveLaw();

    /**
     * Constructor with nested constitutive law.
     */
    TotalLagrangianNumericalTangentBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw);

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new TotalLagrangianNumericalTangentBridgingConstitutiveLaw(mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~TotalLagrangianNumericalTangentBridgingConstitutiveLaw();

    /**
     * Operators
     */

    /**
     * Operations
     */

    /**
     * returns the size of the strain vector of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    SizeType GetStrainSize() const final
    {
        return mpConstitutiveLaw->GetStrainSize();
    }

    bool IsIncremental() final
    {
        return mpConstitutiveLaw->IsIncremental();
    }

    ConstitutiveLaw::StrainMeasure GetStrainMeasure() final
    {
        return ConstitutiveLaw::StrainMeasure_GreenLagrange;
    }

    ConstitutiveLaw::StressMeasure GetStressMeasure() final
    {
        return ConstitutiveLaw::StressMeasure_PK2;
    }

    void GetLawFeatures(Features& rFeatures) final
    {
        rFeatures.SetStrainMeasure(this->GetStrainMeasure());
    }

    bool Has( const Variable<int>& rThisVariable ) final;
    bool Has( const Variable<double>& rThisVariable ) final;
    bool Has( const Variable<Vector>& rThisVariable ) final;
    bool Has( const Variable<Matrix>& rThisVariable );

    int& GetValue( const Variable<int>& rThisVariable, int& rValue ) final;
    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) final;
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue ) final;
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue ) final;

    void SetValue( const Variable<int>& rThisVariable, const int& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final;
    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final;
    void SetValue( const Variable<array_1d<double, 3 > >& rThisVariable, const array_1d<double, 3 > & rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final;
    void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final;
    void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final;
    void SetValue( const Variable<ConstitutiveLaw::Pointer>& rThisVariable, ConstitutiveLaw::Pointer rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final;

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& props,
                             const GeometryType& geom,
                             const Vector& ShapeFunctionsValues ) final;

    /**
     * As this constitutive law describes only linear elastic material properties
     * this function is rather useless and in fact does nothing
     */
    void InitializeSolutionStep( const Properties& props,
                                 const GeometryType& geom, //this is just to give the array of nodes
                                 const Vector& ShapeFunctionsValues,
                                 const ProcessInfo& CurrentProcessInfo ) final;

    void InitializeNonLinearIteration( const Properties& props,
                                       const GeometryType& geom, //this is just to give the array of nodes
                                       const Vector& ShapeFunctionsValues,
                                       const ProcessInfo& CurrentProcessInfo ) final;

    void ResetMaterial( const Properties& props,
                        const GeometryType& geom,
                        const Vector& ShapeFunctionsValues ) final;

    void FinalizeNonLinearIteration( const Properties& props,
                                     const GeometryType& geom, //this is just to give the array of nodes
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo ) final;

    void FinalizeSolutionStep( const Properties& props,
                               const GeometryType& geom, //this is just to give the array of nodes
                               const Vector& ShapeFunctionsValues,
                               const ProcessInfo& CurrentProcessInfo ) final;

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    int Check( const Properties& props,
               const GeometryType& geom,
               const ProcessInfo& CurrentProcessInfo ) const final;

    /**
     * Computes the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(Parameters& rValues) final;

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "TotalLagrangianNumericalTangentBridgingConstitutiveLaw";
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const final
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

    ConstitutiveLaw::Pointer mpConstitutiveLaw;

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw );
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

    /**
     * Assignment operator.
     */
    //TotalLagrangianNumericalTangentBridgingConstitutiveLaw& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
    /**
     * Copy constructor.
     */
    //TotalLagrangianNumericalTangentBridgingConstitutiveLaw(const IsotropicPlaneStressWrinklingNew& rOther);
}; // Class TotalLagrangianNumericalTangentBridgingConstitutiveLaw

} // namespace Kratos.

#endif // KRATOS_TOTAL_LAGRANGIAN_NUMERICAL_TANGENT_BRIDGING_CONSTITUTIVE_LAW_INCLUDED  defined
