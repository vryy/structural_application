/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 5 Oct 2022 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

#if !defined(KRATOS_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED )
#define  KRATOS_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "custom_utilities/sd_math_utils.h"


namespace Kratos
{

/**
 * A wrapper/bridging constitutive law to use with the FiniteStrain element
 */
class KRATOS_API(STRUCTURAL_APPLICATION) FiniteStrainBridgingConstitutiveLaw : virtual public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;
    typedef BaseType::SizeType SizeType;
    typedef SD_MathUtils<double>::Fourth_Order_Tensor Fourth_Order_Tensor;

    /**
     * Counted pointer of FiniteStrainBridgingConstitutiveLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(FiniteStrainBridgingConstitutiveLaw);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    FiniteStrainBridgingConstitutiveLaw();

    /**
     * Constructor with nested constitutive law.
     */
    FiniteStrainBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw);

    ConstitutiveLaw::Pointer Clone() const override
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /**
     * Destructor.
     */
    virtual ~FiniteStrainBridgingConstitutiveLaw();

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
    SizeType GetStrainSize() const override
    {
        return mpConstitutiveLaw->GetStrainSize();
    }

    bool IsIncremental() override
    {
        return mpConstitutiveLaw->IsIncremental();
    }

    ConstitutiveLaw::StrainMeasure GetStrainMeasure() override
    {
        return ConstitutiveLaw::StrainMeasure_Deformation_Gradient;
    }

    ConstitutiveLaw::StressMeasure GetStressMeasure() override
    {
        return ConstitutiveLaw::StressMeasure_Cauchy;
    }

    void GetLawFeatures(Features& rFeatures) override
    {
        /// Since the constitutive law is bridging one, it supports all types of strain measure required by the elemment.
        /// Nevertheless, it mainly support the StrainMeasure_Deformation_Gradient
        rFeatures.SetStrainMeasure(ConstitutiveLaw::StrainMeasure_Infinitesimal);
        rFeatures.SetStrainMeasure(ConstitutiveLaw::StrainMeasure_GreenLagrange);
        rFeatures.SetStrainMeasure(ConstitutiveLaw::StrainMeasure_Almansi);
        rFeatures.SetStrainMeasure(ConstitutiveLaw::StrainMeasure_Hencky_Material);
        rFeatures.SetStrainMeasure(ConstitutiveLaw::StrainMeasure_Hencky_Spatial);
        rFeatures.SetStrainMeasure(ConstitutiveLaw::StrainMeasure_Deformation_Gradient);
        rFeatures.SetStrainMeasure(ConstitutiveLaw::StrainMeasure_Right_CauchyGreen);
        rFeatures.SetStrainMeasure(ConstitutiveLaw::StrainMeasure_Left_CauchyGreen);
    }

    bool Has( const Variable<int>& rThisVariable ) override;
    bool Has( const Variable<double>& rThisVariable ) override;
    bool Has( const Variable<Vector>& rThisVariable ) override;
    bool Has( const Variable<Matrix>& rThisVariable ) override;

    int& GetValue( const Variable<int>& rThisVariable, int& rValue ) override;
    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue ) override;
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue ) override;

    Matrix& CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

    void SetValue( const Variable<int>& rThisVariable, const int& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<array_1d<double, 3 > >& rThisVariable, const array_1d<double, 3 > & rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<ConstitutiveLaw::Pointer>& rThisVariable, ConstitutiveLaw::Pointer rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;

    void InitializeMaterial( const Properties& props,
                             const GeometryType& geom,
                             const Vector& ShapeFunctionsValues ) override;

    void InitializeSolutionStep( const Properties& props,
                                 const GeometryType& geom, //this is just to give the array of nodes
                                 const Vector& ShapeFunctionsValues,
                                 const ProcessInfo& CurrentProcessInfo ) override;

    void InitializeNonLinearIteration( const Properties& props,
                                       const GeometryType& geom, //this is just to give the array of nodes
                                       const Vector& ShapeFunctionsValues,
                                       const ProcessInfo& CurrentProcessInfo ) override;

    void ResetMaterial( const Properties& props,
                        const GeometryType& geom,
                        const Vector& ShapeFunctionsValues ) override;

    void FinalizeNonLinearIteration( const Properties& props,
                                     const GeometryType& geom, //this is just to give the array of nodes
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo ) override;

    void FinalizeSolutionStep( const Properties& props,
                               const GeometryType& geom, //this is just to give the array of nodes
                               const Vector& ShapeFunctionsValues,
                               const ProcessInfo& CurrentProcessInfo ) override;

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
               const ProcessInfo& CurrentProcessInfo ) const override;

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const override
    {
        std::stringstream ss;
        ss << "FiniteStrainBridgingConstitutiveLaw("
           << mpConstitutiveLaw->Info() << ")";
        return ss.str();
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
    void PrintData(std::ostream& rOStream) const override
    {
        mpConstitutiveLaw->PrintData(rOStream);
    }

protected:
    /**
     * Member Variables
     */

    ConstitutiveLaw::Pointer mpConstitutiveLaw;

    Matrix m_F_n;
    Matrix m_F_n1;
    double m_J_n1;

    Matrix m_stress_n1; // Cauchy stress

    /// Update the internal deformation gradient
    void UpdateDeformationGradient(Matrix& F, double& J, const Parameters& rValues) const;

    /// Integrate the stress
    virtual void StressIntegration(const Parameters& rValues, const Matrix& F, Matrix& stress_tensor) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Compute the consistent tangent (tensor A)
    virtual void ComputeTangent(Fourth_Order_Tensor& A) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Compute the consistent tangent (tensor A)
    virtual void ComputeTangent(const Parameters& rValues, Fourth_Order_Tensor& A) const
    {
        this->ComputeTangent(A);
    }

    /// Compute the consistent tangent (matrix)
    virtual void ComputeTangent(Matrix& AlgorithmicTangent) const;

    /// Compute the consistent tangent (matrix)
    virtual void ComputeTangent(const Parameters& rValues, Matrix& AlgorithmicTangent) const;

    /// Get the strain size
    virtual unsigned int GetStrainSize(unsigned int dim) const
    {
        return dim*(dim+1) / 2;
    }

    /// Compute the numerical derivatives of stress with respect to F
    /// Remark: the current implementation does not work with Fbar
    void ComputeStressDerivatives(Fourth_Order_Tensor& D, const Parameters& rValues, const Matrix& F, double epsilon) const;

private:

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
}; // Class FiniteStrainBridgingConstitutiveLaw

} // namespace Kratos.

#endif // KRATOS_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED  defined
