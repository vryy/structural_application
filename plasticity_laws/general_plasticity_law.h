/*
==============================================================================
see soil_mechanics_appplication/LICENSE.txt
==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Giang Bui-Hoang $
//   Date:                $Date: 16 Feb 2020 $
//   Revision:            $Revision: 1.0 $
//
//



#if !defined(KRATOS_SOIL_MECHANICS_GENERAL_PLASTICITY_LAW_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_GENERAL_PLASTICITY_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/sd_math_utils.h"

// #define DEBUG_GENERAL_PLASTICITY_LAW

namespace Kratos
{

/**
 * Asbtract class for the general plasticity law.
 * It is noted that q denotes the stress-like internal variables, not deviatoric pressure
 * REF: gen_plas.pdf
 */
class KRATOS_API(STRUCTURAL_APPLICATION) GeneralPlasticityLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(GeneralPlasticityLaw);
    typedef SD_MathUtils<double>::Third_Order_Tensor Third_Order_Tensor;
    typedef SD_MathUtils<double>::Fourth_Order_Tensor Fourth_Order_Tensor;

    /**
     * Constructor.
     */
    GeneralPlasticityLaw()
    {}

    /**
     * Destructor.
     */
    virtual ~GeneralPlasticityLaw()
    {}

    /**
     * Operations
     */

    virtual bool Has( const Variable<int>& rThisVariable ) const
    {
        return false;
    }

    virtual bool Has( const Variable<double>& rThisVariable ) const
    {
        return false;
    }

    virtual bool Has( const Variable<Vector>& rThisVariable ) const
    {
        return false;
    }

    virtual bool Has( const Variable<Matrix>& rThisVariable ) const
    {
        return false;
    }

    virtual double& GetValue( const Matrix& stress, const Vector& q, const Vector& alpha,
        const Variable<double>& rThisVariable, double& rValue ) const
    {
        return rValue;
    }

    virtual Vector& GetValue( const Matrix& stress, const Vector& q, const Vector& alpha,
        const Variable<Vector>& rThisVariable, Vector& rValue ) const
    {
        return rValue;
    }

    virtual Matrix& GetValue( const Matrix& stress, const Vector& q, const Vector& alpha,
        const Variable<Matrix>& rThisVariable, Matrix& rValue ) const
    {
        return rValue;
    }

    virtual void SetValue( Matrix& stress, Vector& q, Vector& alpha,
        const Variable<double>& rThisVariable, const double& rValue,
        const ProcessInfo& rCurrentProcessInfo )
    {}

    virtual void SetValue( Matrix& stress, Vector& q, Vector& alpha,
        const Variable<Vector>& rThisVariable, const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo )
    {}

    virtual void SetValue( Matrix& stress, Vector& q, Vector& alpha,
        const Variable<Matrix>& rThisVariable, const Matrix& rValue,
        const ProcessInfo& rCurrentProcessInfo )
    {}

    /**
     * Material parameters are initialized
     */
    virtual void InitializeMaterial( Matrix& stress, Vector& q, Vector& alpha, const Properties& props )
    {}

    /**
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param props the Properties instance of the current element
     * @param geom the geometry of the current element
     * @param ShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void ResetMaterial( Matrix& stress, Vector& q, Vector& alpha, const Properties& props )
    {}

    ///////////////////////////////////////////////////////////
    ////////////////// PLASTIC INTERFACE //////////////////////
    ///////////////////////////////////////////////////////////

    /**
     * returns the number of internal variables associated with the constitutive law
     */
    virtual int NumberOfInternalVariables() const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Yield surface
    virtual double F(const Matrix& stress, const Vector& q, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Derivatives of yield surface w.r.t stress
    virtual void dFdSigma(Matrix& n, const Matrix& stress, const Vector& q, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Hessian of yield surface w.r.t stress
    virtual void d2FdSigma2(Fourth_Order_Tensor& dn_dsigma, const Matrix& stress, const Vector& q, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Derivatives of yield surface w.r.t internal variables
    virtual void dFdQ(Vector& dfdq, const Matrix& stress, const Vector& q, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Plastic potential
    virtual double G(const Matrix& stress, const Vector& q, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Hessian of plastic potential w.r.t stress
    virtual void d2GdSigma2(Fourth_Order_Tensor& dn_dsigma, const Matrix& stress, const Vector& q, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Derivatives of plastic potential w.r.t stress
    virtual void dGdSigma(Matrix& m, const Matrix& stress, const Vector& q, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Derivatives of plastic potential w.r.t stress and thermodynamics forces
    virtual void d2GdSigmadQ(Third_Order_Tensor& dmdq, const Matrix& stress, const Vector& q, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Derivatives of plastic potential w.r.t internal variables
    virtual void dGdQ(Vector& dgdq, const Matrix& stress, const Vector& q, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Hardening matrix. It is the derivatives of q w.r.t alpha
    virtual void dQdAlpha(Matrix& Hmat, const Matrix& stress, const Vector& q, const Vector& alpha, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Derivatives of strain-like parameter w.r.t plastic multiplier
    virtual void dAlphadPhi(Vector& dalpha_dphi, const Vector& alpha, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// The hardening rule in stress space, i.e. derivatives of stress-like parameter w.r.t plastic multiplier
    /// In the technical note gen_plas.pdf, it is computed as (-d)
    /// This can also be computed as dQdAlpha * dAlphadPhi
    virtual void dQdPhi(Vector& dq_dphi, const Matrix& stress, const Vector& q, const Vector& alpha, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Give/Compute a reference yield value. This value is used to check the yield condition during plastic integration.
    virtual double ReferenceYieldValue(const Matrix& stress, const Vector& q, const Vector& alpha, const ProcessInfo& CurrentProcessInfo, const Properties& props) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    ///////////////////////////////////////////////////////////
    ///////////// PLASTIC CALCULATION SUBROUTINES /////////////
    ///////////////////////////////////////////////////////////

    /// Compute the first yield point
    double ComputeFirstYieldPoint(const Matrix& stress0, const Matrix& delta_stress,
        const Vector& q, const int ndiv, const double NTOL, const double FTOL,
        const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    /// Return the stress to the yield surface
    int DriftCorrection(Matrix& stress, Vector& q, Vector& alpha, const Fourth_Order_Tensor& Ce,
        const double FTOL, const int max_iters,
        const ProcessInfo& CurrentProcessInfo, const Properties& props,
        const int debug_level=0) const;

    ///////////////////////////////////////////////////////////

    virtual int Check( const Properties& props,
                       const ProcessInfo& CurrentProcessInfo) const
    {
        return 0;
    }

    /**
     * Turn back information as a string.
     */
    virtual std::string Name() const
    {
        return "GeneralPlasticityLaw";
    }

    /**
     * Print information about this object.
     */
    virtual void Print(std::ostream& rOStream) const
    {
        rOStream << Name();
    }

    ///////////////////////////////////////////////////////////
    /////////////////// DEBUG SUBROUTINES /////////////////////
    ///////////////////////////////////////////////////////////

    /// Compute dFdSigma numerically using finite difference scheme
    void Num_dFdSigma(Matrix& n, const Matrix& stress, const Vector& q, const double epsilon,
            const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    /// Compute d2FdSigma2 numerically using finite difference scheme
    void Num_d2FdSigma2(Fourth_Order_Tensor& dn_dsigma, const Matrix& stress, const Vector& q, const double epsilon,
            const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    /// Compute dGdSigma numerically using finite difference scheme
    void Num_dGdSigma(Matrix& m, const Matrix& stress, const Vector& q, const double epsilon,
            const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    /// Compute d2GdSigma2 numerically using finite difference scheme
    void Num_d2GdSigma2(Fourth_Order_Tensor& dm_dsigma, const Matrix& stress, const Vector& q, const double epsilon,
            const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    /// Compute dFdQ numerically using finite difference scheme
    void Num_dFdQ(Vector& dfdq, const Matrix& stress, const Vector& q, const double epsilon,
            const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

    /// Compute d2GdSigmadQ numerically using finite difference scheme
    void Num_d2GdSigmadQ(Third_Order_Tensor& dmdq, const Matrix& stress, const Vector& q, const double epsilon,
            const ProcessInfo& CurrentProcessInfo, const Properties& props) const;

private:

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }

    ///@}

}; /* Class GeneralPlasticityLaw */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_GENERAL_PLASTICITY_LAW_H_INCLUDED  defined */
