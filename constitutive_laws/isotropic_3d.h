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
 *   Last Modified by:    $Author: janosch $
 *   Date:                $Date: 2008-01-25 08:37:31 $
 *   Revision:            $Revision: 1.10 $
 *
 * ***********************************************************/

#if !defined(KRATOS_ISOTROPIC_3D_H_INCLUDED )
#define  KRATOS_ISOTROPIC_3D_H_INCLUDED

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
 */
template<class TNodeType>
class KRATOS_API(STRUCTURAL_APPLICATION) Isotropic3DImpl : public ConstitutiveLawImpl<TNodeType>
{
    public:
        /**
         * Type Definitions
         */
        typedef ConstitutiveLawImpl<TNodeType> BaseType;

        typedef typename BaseType::GeometryType GeometryType;
        typedef typename BaseType::DataType DataType;
        typedef typename BaseType::SizeType SizeType;
        typedef typename BaseType::VectorType VectorType;
        typedef typename BaseType::MatrixType MatrixType;

        typedef typename MatrixVectorTypeSelector<DataType>::ZeroVectorType ZeroVectorType;
        typedef typename MatrixVectorTypeSelector<DataType>::ZeroMatrixType ZeroMatrixType;

        /**
         * Counted pointer of Isotropic3D
         */
        KRATOS_CLASS_POINTER_DEFINITION(Isotropic3DImpl);

        /**
         * Life Cycle
         */
        /**
         * Default constructor.
         */
        Isotropic3DImpl();

        /**
         * Destructor.
         */
        ~Isotropic3DImpl() override;

        /**
         * Operators
         */

        /**
         * Operations
         */

        typename BaseType::Pointer Clone() const override
        {
            typename BaseType::Pointer p_clone( new Isotropic3DImpl() );
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

        bool Has( const Variable<int>& rThisVariable ) override;
        bool Has( const Variable<DataType>& rThisVariable ) override;
        bool Has( const Variable<VectorType>& rThisVariable ) override;
        bool Has( const Variable<MatrixType>& rThisVariable ) override;

        int& GetValue( const Variable<int>& rThisVariable, int& rValue ) override;
        DataType& GetValue( const Variable<DataType>& rThisVariable, DataType& rValue ) override;
        VectorType& GetValue( const Variable<VectorType>& rThisVariable, VectorType& rValue ) override;
        MatrixType& GetValue( const Variable<MatrixType>& rThisVariable, MatrixType& rValue ) override;

        void SetValue( const Variable<int>& rThisVariable, const int& rValue,
                       const ProcessInfo& rCurrentProcessInfo ) override;
        void SetValue( const Variable<DataType>& rThisVariable, const DataType& rValue,
                       const ProcessInfo& rCurrentProcessInfo ) override;
        void SetValue( const Variable<array_1d<DataType, 3 > >& rThisVariable,
                       const array_1d<DataType, 3 > & rValue, const ProcessInfo& rCurrentProcessInfo ) override;
        void SetValue( const Variable<VectorType>& rThisVariable, const VectorType& rValue,
                       const ProcessInfo& rCurrentProcessInfo ) override;
        void SetValue( const Variable<MatrixType>& rThisVariable, const MatrixType& rValue,
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
                                      const VectorType& GreenLagrangeStrainVector ) override;

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
         * Computes the material response in terms of Cauchy stresses and constitutive tensor
         * @see Parameters
         */
        void CalculateMaterialResponseCauchy (typename BaseType::Parameters& rValues) final;

        /// DEPRECATED function
        void CalculateMaterialResponse( const VectorType& StrainVector,
                                        const MatrixType& DeformationGradient,
                                        VectorType& StressVector,
                                        MatrixType& AlgorithmicTangent,
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
        SizeType GetStrainSize() const final
        {
            return 6;
        }

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
        void CalculateStrain(const DataType E, const DataType NU, const VectorType& StressVector, VectorType& rResult) const;

        /**
         * calculates the linear elastic constitutive matrix in terms of Young's modulus and
         * Poisson ratio
         * @param E the Young's modulus
         * @param NU the Poisson ratio
         * @return the linear elastic constitutive matrix
         */
        static void CalculateElasticMatrix( MatrixType& C, const DataType E, const DataType NU );

        /**
         * Input and output
         */

        /**
         * Turn back information as a string.
         */
        std::string Info() const override
        {
            return "Isotropic3DImpl";
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

        void save( Serializer& rSerializer ) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
            rSerializer.save( "Prestress", mPrestress );
            rSerializer.save( "PrestressFactor", mPrestressFactor );
            rSerializer.save( "CurrentStress", mCurrentStress );
            rSerializer.save( "mE", mE );
            rSerializer.save( "mNU", mNU );
            rSerializer.save( "mDE", mDE );
        }

        void load( Serializer& rSerializer ) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
            rSerializer.load( "Prestress", mPrestress );
            rSerializer.load( "PrestressFactor", mPrestressFactor );
            rSerializer.load( "CurrentStress", mCurrentStress );
            rSerializer.load( "mE", mE );
            rSerializer.load( "mNU", mNU );
            rSerializer.load( "mDE", mDE );
        }

        ///@}

        /**
         * Static Member Variables
         */

        /**
         * Calculates the stresses for given strain state
         * @param StrainVector the current vector of strains
         * @param rResult the stress vector corresponding to the given strains
         */
        void CalculateStress( const VectorType& StrainVector, const MatrixType& AlgorithmicTangent, VectorType& rResult );

        VectorType mPrestress;
        DataType mPrestressFactor;
        DataType mE, mNU, mDE;

        /**
         * Un accessible methods
         */
        /**
         * Assignment operator.
         */
        //Isotropic3DImpl& operator=(const Isotropic3DImpl& rOther);
        /**
         * Copy constructor.
         */
        //Isotropic3DImpl(const Isotropic3DImpl& rOther);
}; // Class Isotropic3DImpl

typedef Isotropic3DImpl<RealNode> Isotropic3D;
typedef Isotropic3DImpl<ComplexNode> ComplexIsotropic3D;
typedef Isotropic3DImpl<GComplexNode> GComplexIsotropic3D;

} // namespace Kratos.

#endif // KRATOS_ISOTROPIC_3D_H_INCLUDED  defined
