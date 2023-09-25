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
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 21 Sep 2023 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined(KRATOS_STRUCTURAL_APPLICATION_ISOTROPIC_3D_DC_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APPLICATION_ISOTROPIC_3D_DC_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{
/**
 * Defines a linear elastic isotropic constitutive law in 3D space.
 * This material law is defined by the parameters E (Young's modulus)
 * and NU (Poisson ratio)
 * As there are no further parameters the functionality is limited
 * to linear elasticity.
 * For finite strain/total Lagrangian. This constitutive law behaves like
 * St-Venant Kirhhoff material.
 * This version of plane strain is defined specifically for displacement control
 * analysis
 */
#ifdef SD_APP_FORWARD_COMPATIBILITY
class KRATOS_API(STRUCTURAL_APPLICATION) Isotropic3DDC
#else
class Isotropic3DDC
#endif
: public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;
    /**
     * Counted pointer of Isotropic3DDC
     */
    KRATOS_CLASS_POINTER_DEFINITION(Isotropic3DDC);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    Isotropic3DDC();

    /**
     * Destructor.
     */
    virtual ~Isotropic3DDC();

    /**
     * Operators
     */

    /**
     * Operations
     */

    ConstitutiveLaw::Pointer Clone() const final
    {
         ConstitutiveLaw::Pointer p_clone(new Isotropic3DDC());
        return p_clone;
    }

    ConstitutiveLaw::StrainMeasure GetStrainMeasure() final
    {
        return StrainMeasure_Infinitesimal;
    }

    ConstitutiveLaw::StressMeasure GetStressMeasure() final
    {
        return StressMeasure_Cauchy;
    }

    void GetLawFeatures(Features& rFeatures) final
    {
        rFeatures.SetStrainMeasure(StrainMeasure_Infinitesimal);
    }

    bool Has( const Variable<int>& rThisVariable );
    bool Has( const Variable<double>& rThisVariable );
    bool Has( const Variable<Vector>& rThisVariable );
    bool Has( const Variable<Matrix>& rThisVariable );

    int& GetValue( const Variable<int>& rThisVariable, int& rValue );
    double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );

    void SetValue( const Variable<int>& rVariable,
                   const int& Value,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<double>& rVariable,
                   const double& Value,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Vector>& rThisVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Matrix>& rThisVariable,
                   const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& props,
                             const GeometryType& geom,
                             const Vector& ShapeFunctionsValues ) final;

    void ResetMaterial( const Properties& props,
                        const GeometryType& geom,
                        const Vector& ShapeFunctionsValues ) final;

    void InitializeNonLinearIteration( const Properties& rMaterialProperties,
                                       const GeometryType& rElementGeometry,
                                       const Vector& rShapeFunctionsValues,
                                       const ProcessInfo& rCurrentProcessInfo ) final;

    void FinalizeNonLinearIteration( const Properties& rMaterialProperties,
                                     const GeometryType& rElementGeometry,
                                     const Vector& rShapeFunctionsValues,
                                     const ProcessInfo& rCurrentProcessInfo ) final;

    void FinalizeSolutionStep ( const Properties& props,
                                const GeometryType& geom, //this is just to give the array of nodes
                                const Vector& ShapeFunctionsValues ,
                                const ProcessInfo& CurrentProcessInfo ) final;

    /**
     * Calculates the stresses for given strain state
     * @param StrainVector the current vector of strains
     * @param rResult the stress vector corresponding to the given strains
     */
    void CalculateStress(const Matrix& strain_tensor, Matrix& stress_tensor) const;

    /**
     * As this constitutive law describes only linear elastic material properties
     * this function is rather useless and in fact does nothing
     */
    /*            void InitializeSolutionStep( const Properties& props,
                        const GeometryType& geom, //this is just to give the array of nodes
                        const Vector& ShapeFunctionsValues ,
                        const ProcessInfo& CurrentProcessInfo);

                void FinalizeSolutionStep( const Properties& props,
                        const GeometryType& geom, //this is just to give the array of nodes
                        const Vector& ShapeFunctionsValues ,
                        const ProcessInfo& CurrentProcessInfo);
    */

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters& rValues) final;

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters& rValues) final;

    /// DEPRECATED interface
    void CalculateMaterialResponse( const Vector& StrainVector,
                                     const Matrix& DeformationGradient,
                                     Vector& StressVector,
                                     Matrix& AlgorithmicTangent,
                                     const ProcessInfo& CurrentProcessInfo,
                                     const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues,
                                     bool CalculateStresses = true,
                                     int CalculateTangent = 1,
                                     bool SaveInternalVariables = true );

    int Check(const Properties& props,
              const GeometryType& geom,
              const ProcessInfo& CurrentProcessInfo) const final;

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "Isotropic3DDC";
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
     * there are no protected class members
     */

private:
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
        rSerializer.save( "mE", mE );
        rSerializer.save( "mNU", mNU );
        rSerializer.save( "mDE", mDE );
        rSerializer.save( "m_stress_n1", m_stress_n1 );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
        rSerializer.load( "mE", mE );
        rSerializer.load( "mNU", mNU );
        rSerializer.load( "mDE", mDE );
        rSerializer.load( "m_stress_n1", m_stress_n1 );
    }

    /**
     * Static Member Variables
     */

    /**
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */
    void CalculateElasticMatrix(Matrix& C, const double& E, const double& NU) const;

    double mE, mNU, mDE;

    Matrix m_strain_n;
    Matrix m_strain_n1;
    Matrix m_stress_n;
    Matrix m_stress_n1;

    double mPrestressFactor;
    Vector mPreStress;


    /**
     * Un accessible methods
     */
    /**
     * Assignment operator.
     */
    //Isotropic3DDC& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
    /**
     * Copy constructor.
     */
    //Isotropic3DDC(const IsotropicPlaneStressWrinklingNew& rOther);
}; // Class Isotropic3DDC
}  // namespace Kratos.
#endif // KRATOS_STRUCTURAL_APPLICATION_ISOTROPIC_3D_DC_H_INCLUDED  defined
