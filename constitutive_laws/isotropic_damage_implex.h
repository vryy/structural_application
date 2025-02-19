/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */
/* *********************************************************
 *
 *   Last Modified by:    $Author: seyedali $
 *   Date:                $Date: 2014-05-15 00:00:00 $
 *   Revision:            $Revision: 0.54 $
 *
 * ***********************************************************/

#if !defined(KRATOS_ISOTROPIC_DAMAGE_IMPLEX_H_INCLUDED )
#define  KRATOS_ISOTROPIC_DAMAGE_IMPLEX_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

/**
 * Defines an IMPL-EX isotropic damage constitutive law in 3D space.
 * This material law is defined by the parameters E (Young's modulus),
 * NU (Poisson's ratio), Ft (Tensile Strength) and GF (Fracture Energy).
 * Remark: with reference to J. Oliver et al. (2008) and the
 * Java version of IMPL-EX algorithm by hbui
 */
class IsotropicDamageIMPLEX : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;
    /**
     * Counted pointer of Isotropic_Damage_Implex
     */
    typedef boost::shared_ptr<IsotropicDamageIMPLEX> Pointer;

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    IsotropicDamageIMPLEX();

    virtual boost::shared_ptr<ConstitutiveLaw> Clone() const
    {
        boost::shared_ptr<ConstitutiveLaw> p_clone( new IsotropicDamageIMPLEX() );
        return p_clone;
    }

    /**
     * Destructor.
     */
    virtual ~IsotropicDamageIMPLEX();

    /**
     * Operators
     */

    /**
     * Operations
     */

    ConstitutiveLaw::StrainMeasure GetStrainMeasure() override
    {
        return StrainMeasure_Infinitesimal;
    }

    ConstitutiveLaw::StressMeasure GetStressMeasure() override
    {
        return StressMeasure_Cauchy;
    }

    void GetLawFeatures(Features& rFeatures) final
    {
        rFeatures.SetStrainMeasure(this->GetStrainMeasure());
    }

    bool Has( const Variable<int>& rThisVariable );
    bool Has( const Variable<double>& rThisVariable ) override;
    bool Has( const Variable<Vector>& rThisVariable ) override;
    bool Has( const Variable<Matrix>& rThisVariable ) override;

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue ) override;
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue ) override;

    void SetValue( const Variable<bool>& rThisVariable, const bool& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<int>& rThisVariable, const int& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<array_1d<double, 3 > >& rThisVariable,
                   const array_1d<double, 3 > & rValue, const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
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
                                 const GeometryType& geom,
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
                               const GeometryType& geom,
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
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters& parameters) final;

    void CalculateMaterialResponse( const Vector& StrainVector,
                                    const Matrix& DeformationGradient,
                                    Vector& StressVector,
                                    Matrix& AlgorithmicTangent,
                                    const ProcessInfo& CurrentProcessInfo,
                                    const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues,
                                    bool CalculateStresses = true,
                                    int CalculateTangent = true,
                                    bool SaveInternalVariables = true
                                  ) override;

    /**
     * returns the size of the strain vector of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    SizeType GetStrainSize() const override
    {
        return 6;
    }

    /**
     * Calculates the elastic constitutive tensor
     * @param rResult elastic tangent operator
     */
    void CalculateElasticMatrix( Matrix& C, const double E, const double NU ) const;

    /**
     * Calculates the softening law based on an internal variable alpha
     * @param rResult softening function H
     */
    double SofteningLaw( const double alpha ) const;

    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    std::string Info() const override
    {
        return "IsotropicDamageIMPLEX";
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
    {}

protected:

    virtual double DamageFunction(const double kappa) const;

    virtual double DamageFunctionDerivative(const double kappa) const;

private:

    double mFt, mGf, mE, mNU, mE_0, mL, mE_f, mD;
    Vector mCurrentStrain;
    Vector mCurrentStress;
    double mAlpha, mAlpha_old, mAlpha_old_old, mdAlpha, mAlpha_alg;
    double mq, mq_old, mq_alg;
    double mDamage_alg;
    Matrix mC_alg;
    double mDeltaTime, mDeltaTime_old;

    double mInitialDamage;
    double mInitialEps;
    int mDamageFlag;

    int mElemId, mGaussId;

    /// Compute kappa, providing damage
    double ComputeKappa(const double d, const double TOL = 1e-10, const int max_iters = 30) const;

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
}; // Class IsotropicDamageIMPLEX

} // namespace Kratos.

#endif // KRATOS_ISOTROPIC_DAMAGE_IMPLEX_H_INCLUDED  defined
