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

class Isotropic3D : public ConstitutiveLaw
{
    public:
        /**
         * Type Definitions
         */
        typedef ConstitutiveLaw BaseType;
        /**
         * Counted pointer of Isotropic3D
         */
        KRATOS_CLASS_POINTER_DEFINITION(Isotropic3D);

        /**
         * Life Cycle
         */
        /**
         * Default constructor.
         */
        Isotropic3D();

        /**
         * Destructor.
         */
        virtual ~Isotropic3D();

        /**
         * Operators
         */

        /**
         * Operations
         */

        ConstitutiveLaw::Pointer Clone() const override
        {
            ConstitutiveLaw::Pointer p_clone( new Isotropic3D() );
            return p_clone;
        }

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
        void SetValue( const Variable<array_1d<double, 3 > >& rThisVariable,
                       const array_1d<double, 3 > & rValue, const ProcessInfo& rCurrentProcessInfo ) override;
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
        void CalculateCauchyStresses( Vector& Cauchy_StressVector,
                                      const Matrix& F,
                                      const Vector& PK2_StressVector,
                                      const Vector& GreenLagrangeStrainVector ) override;

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
        void CalculateMaterialResponseCauchy (Parameters& rValues) final;

        /// DEPRECATED function
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
        SizeType GetStrainSize() const final
        {
            return 6;
        }

        /**
         * Calculates the stresses for given strain state
         * @param StrainVector the current vector of strains
         * @param rResult the stress vector corresponding to the given strains
         */
        void CalculateStress(const double& E, const double& NU, const Vector& StrainVector, Vector& rResult) const;

        /**
         * Calculates the strain for given stress state
         * @param StrainVector the current vector of strains
         * @param rResult the stress vector corresponding to the given strains
         */
        void CalculateStrain(const double& E, const double& NU, const Vector& StressVector, Vector& rResult) const;

        /**
         * calculates the linear elastic constitutive matrix in terms of Young's modulus and
         * Poisson ratio
         * @param E the Young's modulus
         * @param NU the Poisson ratio
         * @return the linear elastic constitutive matrix
         */
        static void CalculateElasticMatrix( Matrix& C, const double& E, const double& NU );

        /**
         * Input and output
         */

        /**
         * Turn back information as a string.
         */
        std::string Info() const override
        {
            return "Isotropic3D";
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

        Vector mCurrentStress;

    private:

        ///@name Serialization
        ///@{

        friend class Serializer;

        void save( Serializer& rSerializer ) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw );
            rSerializer.save( "Prestress", mPrestress );
            rSerializer.save( "PrestressFactor", mPrestressFactor );
            rSerializer.save( "CurrentStress", mCurrentStress );
            rSerializer.save( "mE", mE );
            rSerializer.save( "mNU", mNU );
            rSerializer.save( "mDE", mDE );
        }

        void load( Serializer& rSerializer ) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw );
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
        void CalculateStress( const Vector& StrainVector, Matrix& AlgorithmicTangent, Vector& rResult );

        Vector mPrestress;
        double mPrestressFactor;
        double mE, mNU, mDE;

        /**
         * Un accessible methods
         */
        /**
         * Assignment operator.
         */
        //Isotropic3D& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
        /**
         * Copy constructor.
         */
        //Isotropic3D(const IsotropicPlaneStressWrinklingNew& rOther);
}; // Class Isotropic3D

} // namespace Kratos.

#endif // KRATOS_ISOTROPIC_3D_H_INCLUDED  defined
