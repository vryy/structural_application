/*
LICENSE: see material_point_application/LICENSE.txt
*/
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
 * Reference: De Souza Neto, Computational Plasticity, Box 14.3
 */
class MultiplicativeFiniteStrainBridgingConstitutiveLaw : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;
    typedef BaseType::SizeType SizeType;
    typedef double (*unitary_func_t)(double);
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
        return ConstitutiveLaw::StrainMeasure_Deformation_Gradient;
    }

    ConstitutiveLaw::StressMeasure GetStressMeasure() final
    {
        return ConstitutiveLaw::StressMeasure_Cauchy;
    }

    void GetLawFeatures(Features& rFeatures) final
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
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters& rValues) final;

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "MultiplicativeFiniteStrainBridgingConstitutiveLaw";
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

    static double exp2(double x)
    {
        return exp(2.0*x);
    }

    static double log2(double x)
    {
        return 0.5*log(x);
    }

    static double dlog2(double x)
    {
        return 0.5/x;
    }

    /*
    * Compute the isotropic function of the type
    *       Y(X) = sum{ y(x_i) E_i }
    * WHERE Y AND X ARE SYMMETRIC TENSORS, x_i AND E_i ARE, RESPECTIVELY
    * THE EIGENVALUES AND EIGENPROJECTIONS OF X, AND y(.) IS A SCALAR
    * FUNCTION.
    * X must be 3 x 3 matrix
    * Y will be resized accordingly
    * e must be sorted
    * Rererence: Section A.5.2, Computational Plasticity, de Souza Neto.
    * TODO to move to SD_MathUtils
    */
    template<typename TMatrixType>
    static void ComputeIsotropicTensorFunction(
        unitary_func_t func,
        TMatrixType& Y,
        const std::vector<double>& e,
        const std::vector<Matrix>& eigprj,
        const double& TOL = 1.0e-10
    )
    {
        if ((Y.size1() != 3) || (Y.size2() != 3))
            Y.resize(3, 3, false);

        Y.clear();
        if (fabs(e[0] - e[1]) > TOL && fabs(e[1] - e[2]) > TOL)
        {
            for (int d = 0; d < 3; ++d)
                noalias(Y) += func(e[d]) * eigprj[d];
        }
        else
        {
            const Matrix eye = IdentityMatrix(3);

            if (fabs(e[0] - e[1]) < TOL && fabs(e[1] - e[2]) < TOL)
            {
                noalias(Y) = func(e[0]) * eye;
            }
            else
            {
                if (fabs(e[0] - e[1]) < TOL)
                {
                    noalias(Y) = func(e[2]) * eigprj[2]
                               + func(e[0]) * (eye - eigprj[0]);
                }
                else // (fabs(e[1] - e[2]) < TOL)
                {
                    noalias(Y) = func(e[0]) * eigprj[0]
                               + func(e[2]) * (eye - eigprj[2]);
                }
            }
        }
    }

    /*
    * Compute the derivative of the isotropic function of the type
    *       Y(X) = sum{ y(x_i) E_i }
    * WHERE Y AND X ARE SYMMETRIC TENSORS, x_i AND E_i ARE, RESPECTIVELY
    * THE EIGENVALUES AND EIGENPROJECTIONS OF X, AND y(.) IS A SCALAR
    * FUNCTION.
    * X must be 3 x 3 matrix
    * e and eigprj are the principle values and eigenprojections of X
    * dYdX will be resized accordingly
    * Rererence: Section A.5.2, Computational Plasticity, de Souza Neto.
    * TODO to move to SD_MathUtils
    */
    template<typename TMatrixType>
    static void ComputeDerivativeIsotropicTensorFunction(
        unitary_func_t func,
        unitary_func_t dfunc,
        Fourth_Order_Tensor& dYdX,
        const TMatrixType& X,
        const std::vector<double>& e,
        const std::vector<TMatrixType>& eigprj,
        const double& TOL = 1.0e-10
    )
    {
        Fourth_Order_Tensor dX2dX;
        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dX2dX);
        TMatrixType I = IdentityMatrix(3);
        for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
                for(int k = 0; k < 3; ++k)
                    for(int l = 0; l < 3; ++l)
                        dX2dX[i][j](k, l) = 0.5 * (I(i, k) * X(l, j)
                                                 + I(i, l) * X(k, j)
                                                 + X(i, k) * I(j, l)
                                                 + X(i, l) * I(k, j));

        Fourth_Order_Tensor Is;
        SD_MathUtils<double>::CalculateFourthOrderSymmetricTensor(Is);

        SD_MathUtils<double>::CalculateFourthOrderZeroTensor(dYdX);

        if (fabs(e[0] - e[1]) > TOL && fabs(e[1] - e[2]) > TOL)
        {
            constexpr int a = 0, b = 1, c = 2;

            double aux1 = func(e[a]) / ((e[a] - e[b]) * (e[a] - e[c]));
            double aux2 = -aux1 * (e[b] + e[c]);
            double aux3 = -aux1 * (e[a] - e[b] + e[a] - e[c]);
            double aux4 = -aux1 * (e[b] - e[c]);

            SD_MathUtils<double>::AddFourthOrderTensor(aux1, dX2dX, dYdX);
            SD_MathUtils<double>::AddFourthOrderTensor(aux2, Is, dYdX);
            SD_MathUtils<double>::OuterProductFourthOrderTensor(aux3, eigprj[a], eigprj[a], dYdX);
            SD_MathUtils<double>::OuterProductFourthOrderTensor(aux4, eigprj[b], eigprj[b], dYdX);
            SD_MathUtils<double>::OuterProductFourthOrderTensor(-aux4, eigprj[c], eigprj[c], dYdX);
            SD_MathUtils<double>::OuterProductFourthOrderTensor(dfunc(e[a]), eigprj[a], eigprj[a], dYdX);
        }
        else
        {
            if (fabs(e[0] - e[1]) < TOL && fabs(e[1] - e[2]) < TOL)
            {
                SD_MathUtils<double>::AddFourthOrderTensor(dfunc(e[0]), Is, dYdX);
            }
            else
            {
                const Matrix eye = IdentityMatrix(3);
                int a, c;

                if (fabs(e[0] - e[1]) < TOL)
                {
                    a = 2;
                    c = 0;
                }
                else // (fabs(e[1] - e[2]) < TOL)
                {
                    a = 0;
                    c = 2;
                }

                const double s1 = (func(e[a]) - func(e[c])) / pow(e[a] - e[c], 2) - dfunc(e[c]) / (e[a] - e[c]);
                const double s2 = 2.0 * e[c] * (func(e[a]) - func(e[c])) / pow(e[a] - e[c], 2) - dfunc(e[c]) * (e[a] + e[c]) / (e[a] - e[c]);
                const double s3 = 2.0 * (func(e[a]) - func(e[c])) / pow(e[a] - e[c], 3) - (dfunc(e[a]) + dfunc(e[c])) / pow(e[a] - e[c], 2);
                const double s4 = e[c] * s3;
                const double s5 = s4;
                const double s6 = e[c]*e[c] * s3;

                SD_MathUtils<double>::AddFourthOrderTensor(s1, dX2dX, dYdX);
                SD_MathUtils<double>::AddFourthOrderTensor(-s2, Is, dYdX);
                SD_MathUtils<double>::OuterProductFourthOrderTensor(-s3, X, X, dYdX);
                SD_MathUtils<double>::OuterProductFourthOrderTensor(s4, X, eye, dYdX);
                SD_MathUtils<double>::OuterProductFourthOrderTensor(s5, eye, X, dYdX);
                SD_MathUtils<double>::OuterProductFourthOrderTensor(-s6, eye, eye, dYdX);
            }
        }
    }

    /**
     * Assignment operator.
     */
    //MultiplicativeFiniteStrainBridgingConstitutiveLaw& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
    /**
     * Copy constructor.
     */
    //MultiplicativeFiniteStrainBridgingConstitutiveLaw(const IsotropicPlaneStressWrinklingNew& rOther);
}; // Class MultiplicativeFiniteStrainBridgingConstitutiveLaw

} // namespace Kratos.

#endif // KRATOS_MULTIPLICATIVE_FINITE_STRAIN_BRIDGING_CONSTITUTIVE_LAW_INCLUDED  defined
