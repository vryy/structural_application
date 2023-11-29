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
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

/**
 * A wrapper constitutive law that enables the use of infinitesimal strain constitutive law in the finite strain setting.
 * This is using the formulation for Isotropic elastic/elasto-plastic materials with logarithmic finite strain extension
 * TStressType is the stress measure used in the integration algorithm. Right now, Cauchy and Kirchhoff stress are supported.
 * Reference:
 *  +   De Souza Neto, Computational Plasticity, Box 14.3
 *  +   tr20_finite_strain.pdf
 */
template<int TStressType = 2> // 1: Cauchy stress, 2: Kirchhoff stress
class MultiplicativeFiniteStrainBridgingConstitutiveLaw : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;
    typedef BaseType::SizeType SizeType;
    typedef SD_MathUtils<double>::Fourth_Order_Tensor Fourth_Order_Tensor;

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
    MultiplicativeFiniteStrainBridgingConstitutiveLaw();

    /**
     * Constructor with nested constitutive law.
     */
    MultiplicativeFiniteStrainBridgingConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw);

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone( new MultiplicativeFiniteStrainBridgingConstitutiveLaw(mpConstitutiveLaw->Clone()) );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~MultiplicativeFiniteStrainBridgingConstitutiveLaw();

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

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& props,
                             const GeometryType& geom,
                             const Vector& ShapeFunctionsValues ) override;

    /**
     * As this constitutive law describes only linear elastic material properties
     * this function is rather useless and in fact does nothing
     */
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
        if (TStressType == 1)
            return "MultiplicativeFiniteStrainBridgingConstitutiveLaw<Cauchy>";
        else if (TStressType == 2)
            return "MultiplicativeFiniteStrainBridgingConstitutiveLaw<Kirchhoff>";
        else
            KRATOS_ERROR << "Unsupported stress type " << TStressType;
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

    ConstitutiveLaw::Pointer mpConstitutiveLaw;

    Matrix m_F_n;
    Matrix m_F_n1;
    double m_J_n1;

    Matrix m_left_elastic_cauchy_green_tensor_trial;
    Matrix m_stress_n1; // Cauchy stress

    double mPrestressFactor;
    Vector mPrestress;

    /// Compute an equivalent strain for integrating with the small strain constitutive law
    void ComputeLogarithmicStrain(Vector& StrainVector, Matrix& left_elastic_cauchy_green_tensor_trial, const Matrix& F) const;

    /// Integrate the stress; depending on the TStressType, the output is Cauchy stress (TStressType=1) or
    /// Kirchhoff stress (TStressType=2)
    void StressIntegration(const Parameters& rValues, const Matrix& F, Matrix& stress_tensor, Matrix& left_elastic_cauchy_green_tensor_trial) const;

    /// Compute the consistent tangent (tensor)
    virtual void ComputeTangent(Fourth_Order_Tensor& A) const;

    /// Compute the consistent tangent (matrix)
    virtual void ComputeTangent(Matrix& AlgorithmicTangent) const;

    /// Get the strain size
    virtual unsigned int GetStrainSize(unsigned int dim) const
    {
        return dim*(dim+1) / 2;
    }

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

    /// Update the internal deformation gradient
    void UpdateDeformationGradient(Matrix& F, double& J, const Parameters& rValues) const;

    /// Compute E from C
    void ComputeGreenLagrangeStrain(Vector& StrainVector, const Matrix& C) const;

    /// Compute Kirchhoff stress (TStressType=1) or Cauchy stress (TStressType=2) from PK2 stress
    void ComputeStress(Matrix& stress_tensor, const Matrix& PK2_stress) const;

    /// Compute tangent matrix A
    void ComputeInfinitesimalTangent(Fourth_Order_Tensor& A) const;

    /// Compute necessary terms to derive the consistent tangent for infinitesimal strain
    void ComputeInfinitesimalTangentTerms(Fourth_Order_Tensor& D, Fourth_Order_Tensor& L, Fourth_Order_Tensor& B) const;

    /// Compute the tangent according to couple PK2 stress and Green-Lagrange strain
    void ComputeGreenLagrangeTangent(Fourth_Order_Tensor& A) const;

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
     * Counted pointer of MultiplicativeFiniteStrainBridgingConstitutiveLaw
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
