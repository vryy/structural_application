/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 27 Nov 2023 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_STRUCTURAL_APP_CONSTITUTIVE_LAW_DC_WRAPPER_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APP_CONSTITUTIVE_LAW_DC_WRAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application_variables.h"


namespace Kratos
{

/**
 * Wrapper to use the constitutive law in Displacement (DC) Control analysis.
 * In different to load control, the material update is performed during FinalizeNonLinearIteration.
 * Hence, the tangent in the first iteration is from the last step. It improves
 * convergence particularly for DC analysis.
 * In order to successfully inherit from this wrapper, the underlying constitutive law
 * must provide StressIntegration and ComputeTangent with appropriate inputs.
 */
template<class TConstitutiveLawType>
class ConstitutiveLawDcWrapper : public TConstitutiveLawType
{
public:
    /**
     * Type Definitions
     */
    typedef TConstitutiveLawType BaseType;
    typedef typename BaseType::GeometryType GeometryType;

    KRATOS_CLASS_POINTER_DEFINITION(ConstitutiveLawDcWrapper);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    ConstitutiveLawDcWrapper() : BaseType() {}

    /**
     * Destructor.
     */
    virtual ~ConstitutiveLawDcWrapper() {}

    /**
     * Operators
     */

    /**
     * Operations
     */

    ConstitutiveLaw::Pointer Clone() const override
    {
        ConstitutiveLaw::Pointer p_clone(new ConstitutiveLawDcWrapper<TConstitutiveLawType>());
        return p_clone;
    }

    bool Has( const Variable<Vector>& rThisVariable ) override
    {
        if ( rThisVariable == CURRENT_STRAIN_VECTOR )
            return true;

        return BaseType::Has(rThisVariable);
    }

    void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override
    {
        if ( rThisVariable == CURRENT_STRAIN_VECTOR )
        {
            mCurrentStrain = rValue; // size could be changed to accommodate for axisymmetric analysis
            BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
            return;
        }

        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }

    /**
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param props the Properties instance of the current element
     * @param geom the geometry of the current element
     * @param ShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void ResetMaterial( const Properties& props,
                        const GeometryType& geom,
                        const Vector& ShapeFunctionsValues ) override
    {
        unsigned int dim = geom.WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;
        mCurrentStrain = ZeroVector(strain_size);
        BaseType::ResetMaterial(props, geom, ShapeFunctionsValues);
    }

    /**
     * Computes the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear strain measure is used)
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param CurrentProcessInfo current ProcessInfo instance
     * @param props the material's Properties object
     * @param geom the element's geometry
     * @param ShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     */
    /// DEPRECATED interface
    void CalculateMaterialResponse( const Vector& StrainVector,
                                    const Matrix& DeformationGradient,
                                    Vector& StressVector,
                                    Matrix& AlgorithmicTangent,
                                    const ProcessInfo& CurrentProcessInfo,
                                    const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues,
                                    bool CalculateStresses,
                                    int CalculateTangent,
                                    bool SaveInternalVariables) override
    {
        if (CurrentProcessInfo[SET_CALCULATE_REACTION])
        {
            if (CalculateStresses)
            {
                this->GetValue(STRESSES, StressVector);
                return;
            }
        }

        if (CalculateTangent)
            this->ComputeTangent(AlgorithmicTangent, CurrentProcessInfo, props);
        if (CalculateStresses)
            this->GetValue(STRESSES, StressVector);
    }

    /**
     * to be called at the end of each step iteration
     * (e.g. from Element::FinalizeNonLinearIteration)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeNonLinearIteration( const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo) override
    {
        double TOL;
        if (props.Has(LOCAL_ERROR_TOLERANCE))
            TOL = props[LOCAL_ERROR_TOLERANCE];
        else
            TOL = 1.0e-10;

        // Here we delay the stres update after the first iteration is finished. The tangent
        // from the previous step is used for the first iteration. This approach is used by HYPLAS.
        this->StressIntegration(mCurrentStrain, TOL, CurrentProcessInfo, props);
    }

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const override
    {
        return BaseType::Info() + "DC";
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
        BaseType::PrintData(rOStream);
    }

protected:

    /// Reset only the material state to initial state
    void ResetState() override
    {
        BaseType::ResetState();
        mCurrentStrain.clear();
    }

private:

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    void load(Serializer& rSerializer) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }

    /**
     * Member Variables
     */

    Vector mCurrentStrain;

    /**
     * Un accessible methods
     */
    /**
     * Assignment operator.
     */
    /**
     * Copy constructor.
     */
}; // Class ConstitutiveLawDcWrapper

}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_APP_CONSTITUTIVE_LAW_DC_WRAPPER_H_INCLUDED  defined
