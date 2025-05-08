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
*   Last Modified by:    $Author: kazem $
*   Date:                $Date: 2009-01-16 10:49:53 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_STRUCTURAL_APPLICATION_PLANE_STRAIN_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APPLICATION_PLANE_STRAIN_H_INCLUDED

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
 * St-Venant Kirhhoff material. (It's not true, it's different material)
 */
template<class TNodeType>
class KRATOS_API(STRUCTURAL_APPLICATION) PlaneStrainImpl : public ConstitutiveLawImpl<TNodeType>
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLawImpl<TNodeType> BaseType;

    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::DataType DataType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;

    typedef typename MatrixVectorTypeSelector<DataType>::ZeroVectorType ZeroVectorType;
    typedef typename MatrixVectorTypeSelector<DataType>::ZeroMatrixType ZeroMatrixType;

    /**
     * Counted pointer of PlaneStrainImpl
     */
    KRATOS_CLASS_POINTER_DEFINITION(PlaneStrainImpl);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    PlaneStrainImpl();

    /**
     * Destructor.
     */
    ~PlaneStrainImpl() override;

    /**
     * Operators
     */

    /**
     * Operations
     */

    typename BaseType::Pointer Clone() const override
    {
        typename BaseType::Pointer p_clone(new PlaneStrainImpl());
        return p_clone;
    }

    typename BaseType::StrainMeasure GetStrainMeasure() override
    {
        return BaseType::StrainMeasure_Infinitesimal;
    }

    typename BaseType::StressMeasure GetStressMeasure() override
    {
        return BaseType::StressMeasure_Cauchy;
    }

    void GetLawFeatures(typename BaseType::Features& rFeatures) final
    {
        rFeatures.SetStrainMeasure(this->GetStrainMeasure());
    }

    std::size_t GetStrainSize() const final;

    bool Has( const Variable<int>& rThisVariable ) override;
    bool Has( const Variable<DataType>& rThisVariable ) override;
    bool Has( const Variable<VectorType>& rThisVariable ) override;
    bool Has( const Variable<MatrixType>& rThisVariable ) override;

    int& GetValue( const Variable<int>& rThisVariable, int& rValue ) override;
    DataType& GetValue( const Variable<DataType>& rThisVariable, DataType& rValue ) override;
    VectorType& GetValue( const Variable<VectorType>& rThisVariable, VectorType& rValue ) override;
    MatrixType& GetValue( const Variable<MatrixType>& rThisVariable, MatrixType& rValue ) override;

    void SetValue( const Variable<int>& rVariable,
                   const int& Value,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<DataType>& rVariable,
                   const DataType& Value,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<VectorType>& rThisVariable,
                   const VectorType& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<MatrixType>& rThisVariable,
                   const MatrixType& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
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

    /**
     * Calculates the constitutive matrix for a given strain vector
     * @param StrainVector the current vector of strains the constitutive
     * matrix is to be generated for
     * @param rResult MatrixType the result will be stored in
     */
    void CalculateConstitutiveMatrix(MatrixType& rResult) const;

    /**
     * Calculates the stresses for given strain state
     * @param StrainVector the current vector of strains
     * @param rResult the stress vector corresponding to the given strains
     */
    void CalculateStress(const VectorType& StrainVector, VectorType& rResult);

    /**
     * Calculates the stresses for given strain state
     * @param StrainVector the current vector of strains
     * @param rResult the stress vector corresponding to the given strains
     */
    void CalculateStress(const DataType E, const DataType NU, const VectorType& StrainVector, VectorType& rResult) const;

    /**
     * Calculates the strain for given stress state
     * @param StrainVector the current vector of strains
     * @param rResult the stress vector corresponding to the given strains
     */
    void CalculateStrain(const DataType E, const DataType NU, const VectorType& StressVector, VectorType& rResult ) const;

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
     * Calculates the cauchy stresses. For a given deformation and stress state
     * the cauchy stress vector is calculated
     * @param Cauchy_StressVector the result vector of cauchy stresses (will be overwritten)
     * @param F the current deformation gradient
     * @param PK2_StressVector the current second Piola-Kirchhoff-Stress vector
     * @param GreenLagrangeStrainVector the current Green-Lagrangian strains
     */
    void CalculateCauchyStresses( VectorType& Cauchy_StressVector,
                                  const MatrixType& F,
                                  const VectorType& PK2_StressVector,
                                  const VectorType& GreenLagrangeStrainVector ) const;

    /**
     * converts a strain vector styled variable into its form, which the
     * deviatoric parts are no longer multiplied by 2
     */
    void Calculate(const Variable<MatrixType >& rVariable, MatrixType& rResult,
                   const ProcessInfo& rCurrentProcessInfo) const;

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy (typename BaseType::Parameters& rValues) final;

    /// DEPRECATED interface
    void CalculateMaterialResponse( const VectorType& StrainVector,
                                    const MatrixType& DeformationGradient,
                                    VectorType& StressVector,
                                    MatrixType& AlgorithmicTangent,
                                    const ProcessInfo& CurrentProcessInfo,
                                    const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues,
                                    bool CalculateStresses = true,
                                    int CalculateTangent = 1,
                                    bool SaveInternalVariables = true ) override;

    int Check(const Properties& props,
              const GeometryType& geom,
              const ProcessInfo& CurrentProcessInfo) const final;

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const override
    {
        return "PlaneStrain";
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
     * there are several protected class members
     */

    VectorType mCurrentStress;

private:

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save( "mE", mE );
        rSerializer.save( "mNU", mNU );
        rSerializer.save( "mDE", mDE );
        rSerializer.save( "mCurrentStress", mCurrentStress );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load( "mE", mE );
        rSerializer.load( "mNU", mNU );
        rSerializer.load( "mDE", mDE );
        rSerializer.load( "mCurrentStress", mCurrentStress );
    }

    ///@}

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
    void CalculateElasticMatrix(MatrixType& C, const DataType E, const DataType NU) const;

    DataType mE, mNU, mDE;
    DataType mPrestressFactor;
    VectorType mPreStress;

    /**
     * Un accessible methods
     */
    /**
     * Assignment operator.
     */
    //PlaneStrainImpl& operator=(const PlaneStrainImpl& rOther);
    /**
     * Copy constructor.
     */
    //PlaneStrainImpl(const PlaneStrainImpl& rOther);
}; // Class PlaneStrainImpl

typedef PlaneStrainImpl<RealNode> PlaneStrain;
typedef PlaneStrainImpl<ComplexNode> ComplexPlaneStrain;
typedef PlaneStrainImpl<GComplexNode> GComplexPlaneStrain;

}  // namespace Kratos.

#endif // KRATOS_ISOTROPIC_2D_H_INCLUDED  defined
