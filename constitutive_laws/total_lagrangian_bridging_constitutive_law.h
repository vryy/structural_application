/*
LICENSE: see material_point_application/LICENSE.txt
*/
/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 16 Sep 2023 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

#if !defined(KRATOS_TOTAL_LAGRANGIAN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED )
#define  KRATOS_TOTAL_LAGRANGIAN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

/**
 * A wrapper constitutive law that enables the use of infinitesimal strain
 * constitutive law with the TotalLagrangian element.
 * Basically the constitutive law perform the conversion of the Green-Lagrange
 * strain input to the logarithmic strain (Hencky). The Cauchy stress output
 * of the constitutive law is then converted to the 2nd PK stress. The tangent
 * is converted (approximately) accordingly.
 * Reference:
 *  + Calculix's umat_abaqusnl.f
 *  + ccx_2.21.pdf
 */
class TotalLagrangianBridgingConstitutiveLaw : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;
    typedef BaseType::SizeType SizeType;
    typedef SD_MathUtils<double>::Fourth_Order_Tensor Fourth_Order_Tensor;

    /**
     * Counted pointer of TotalLagrangianBridgingConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(TotalLagrangianBridgingConstitutiveLaw);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    TotalLagrangianBridgingConstitutiveLaw();

    /**
     * Constructor with nested constitutive law.
     */
    TotalLagrangianBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw);

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new TotalLagrangianBridgingConstitutiveLaw(mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~TotalLagrangianBridgingConstitutiveLaw();

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
        return "TotalLagrangianBridgingConstitutiveLaw";
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

    Matrix mLastF;
    Matrix mCurrentF;

    /**
     * Un accessible methods
     */

    /**
     * Assignment operator.
     */
    //TotalLagrangianBridgingConstitutiveLaw& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
    /**
     * Copy constructor.
     */
    //TotalLagrangianBridgingConstitutiveLaw(const IsotropicPlaneStressWrinklingNew& rOther);
}; // Class TotalLagrangianBridgingConstitutiveLaw

} // namespace Kratos.

#endif // KRATOS_TOTAL_LAGRANGIAN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED  defined
