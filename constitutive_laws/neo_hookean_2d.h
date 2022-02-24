/*
LICENSE: see soil_mechanics_application/LICENSE.txt
*/
/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 4 Mar 2016 $
 *   Revision:            $Revision: 1.0 $
 *
 * ***********************************************************/

#if !defined(KRATOS_NEO_HOOKEAN_2D_H_INCLUDED )
#define  KRATOS_NEO_HOOKEAN_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

/**
 * Defines a neo-hookean constitutive law for plane strain.
 * This material law is defined by the parameters E (Young's modulus)
 * and NU (Poisson ratio)
 * * TODO: Reference??
 */
class NeoHookean2D : public ConstitutiveLaw
{
    public:
        /**
         * Type Definitions
         */
        typedef ConstitutiveLaw BaseType;
        /**
         * Counted pointer of NeoHookean2D
         */
        KRATOS_CLASS_POINTER_DEFINITION(NeoHookean2D);

        /**
         * Life Cycle
         */
        /**
         * Default constructor.
         */
        NeoHookean2D();

        /**
         * Destructor.
         */
        virtual ~NeoHookean2D();

        /**
         * Operators
         */

        /**
         * Operations
         */

        ConstitutiveLaw::Pointer Clone() const final
        {
            ConstitutiveLaw::Pointer p_clone( new NeoHookean2D() );
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

        bool Has( const Variable<int>& rThisVariable );
        bool Has( const Variable<double>& rThisVariable );
        bool Has( const Variable<Vector>& rThisVariable );
        bool Has( const Variable<Matrix>& rThisVariable );

        int& GetValue( const Variable<int>& rThisVariable, int& rValue );
        double& GetValue( const Variable<double>& rThisVariable, double& rValue );
        Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
        Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );

        void SetValue( const Variable<int>& rThisVariable, const int& rValue,
                       const ProcessInfo& rCurrentProcessInfo );
        void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                       const ProcessInfo& rCurrentProcessInfo );
        void SetValue( const Variable<array_1d<double, 3 > >& rThisVariable,
                       const array_1d<double, 3 > & rValue, const ProcessInfo& rCurrentProcessInfo );
        void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                       const ProcessInfo& rCurrentProcessInfo );
        void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                       const ProcessInfo& rCurrentProcessInfo );

        /**
         * Material parameters are inizialized
         */
        void InitializeMaterial( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues );

        /**
         * As this constitutive law describes only linear elastic material properties
         * this function is rather useless and in fact does nothing
         */
        void InitializeSolutionStep( const Properties& props,
                                     const GeometryType& geom, //this is just to give the array of nodes
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo );

        void InitializeNonLinearIteration( const Properties& props,
                                           const GeometryType& geom, //this is just to give the array of nodes
                                           const Vector& ShapeFunctionsValues,
                                           const ProcessInfo& CurrentProcessInfo );

        void ResetMaterial( const Properties& props,
                            const GeometryType& geom,
                            const Vector& ShapeFunctionsValues );

        void FinalizeNonLinearIteration( const Properties& props,
                                         const GeometryType& geom, //this is just to give the array of nodes
                                         const Vector& ShapeFunctionsValues,
                                         const ProcessInfo& CurrentProcessInfo );

        void FinalizeSolutionStep( const Properties& props,
                                   const GeometryType& geom, //this is just to give the array of nodes
                                   const Vector& ShapeFunctionsValues,
                                   const ProcessInfo& CurrentProcessInfo );

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
                                      const Vector& GreenLagrangeStrainVector );

        /**
         * This function is designed to be called once to perform all the checks needed
         * on the input provided. Checks can be "expensive" as the function is designed
         * to catch user's errors.
         * @param props
         * @param geom
         * @param CurrentProcessInfo
         * @return
         */
        virtual int Check( const Properties& props,
                           const GeometryType& geom,
                           const ProcessInfo& CurrentProcessInfo );

        /**
         * Computes the material response in terms of Cauchy stresses and constitutive tensor
         * @see Parameters
         */
        void CalculateMaterialResponseCauchy (Parameters& rValues) final;

        /// DEPRECATED FUNCTION
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
                                      );

        /**
         * returns the size of the strain vector of the current constitutive law
         * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
         */
        SizeType GetStrainSize() const final
        {
            return 3;
        }

        /**
         * converts a strain vector styled variable into its form, which the
         * deviatoric parts are no longer multiplied by 2
         */
        //             void Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, const ProcessInfo& rCurrentProcessInfo);

        /**
         * Input and output
         */

        /**
         * Turn back information as a string.
         */
        std::string Info() const final
        {
            return "NeoHookean2D";
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

        virtual void save( Serializer& rSerializer ) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw );
            rSerializer.save( "Prestress", mPrestress );
            rSerializer.save( "PrestressFactor", mPrestressFactor );
            rSerializer.save( "CurrentStress", mCurrentStress );
            rSerializer.save( "mE", mE );
            rSerializer.save( "mNU", mNU );
            rSerializer.save( "mDE", mDE );
        }

        virtual void load( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw );
            rSerializer.load( "Prestress", mPrestress );
            rSerializer.load( "PrestressFactor", mPrestressFactor );
            rSerializer.load( "CurrentStress", mCurrentStress );
            rSerializer.load( "mE", mE );
            rSerializer.load( "mNU", mNU );
            rSerializer.load( "mDE", mDE );
        }

        /**
         * Static Member Variables
         */

        /**
         * Calculates the stresses for given strain state
         * @param StrainVector the current vector of strains
         * @param rResult the stress vector corresponding to the given strains
         */
        void CalculateStress( Vector& StressVector, const Vector& StrainVector );

        /**
         * calculates the linear elastic constitutive matrix in terms of Young's modulus and
         * Poisson ratio
         * @param E the Young's modulus
         * @param NU the Poisson ratio
         * @return the linear elastic constitutive matrix
         */
        void CalculateTangentMatrix( Matrix& C, const Vector& StrainVector );

        Vector mPrestress;
        double mPrestressFactor;
        Vector mCurrentStress;
        double mE, mNU, mDE;

        /**
         * Un accessible methods
         */
        /**
         * Assignment operator.
         */
        //NeoHookean2D& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
        /**
         * Copy constructor.
         */
        //NeoHookean2D(const IsotropicPlaneStressWrinklingNew& rOther);
}; // Class NeoHookean2D

} // namespace Kratos.
#endif // KRATOS_NEO_HOOKEAN_3D_H_INCLUDED  defined
