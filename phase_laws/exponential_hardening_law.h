//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Dec 2018 $
//   Revision:            $Revision: 1.0 $
//
//



#if !defined(KRATOS_STRUCTURAL_APP_EXPONENTIAL_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APP_EXPONENTIAL_HARDENING_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "phase_laws/hardening_law.h"


namespace Kratos
{

/**
 * Exponential hardening law of the form q = c_inf - (c_inf-c_0)*exp(-rho*alpha)
 */
class KRATOS_API(STRUCTURAL_APPLICATION) ExponentialHardeningLaw : public HardeningLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(ExponentialHardeningLaw);

    typedef HardeningLaw BaseType;

    /**
     * Constructor.
     */
    ExponentialHardeningLaw();
    ExponentialHardeningLaw(const double C0, const double Cinf, const double Rho);

    /**
     * Destructor.
     */
    virtual ~ExponentialHardeningLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      ExponentialHardeningLaw::Pointer p_clone(new ExponentialHardeningLaw());
     *      return p_clone;
     */
    BaseType::Pointer Clone() const final;

    /**
     * Operations
     */
    bool Has( const Variable<double>& rThisVariable ) final;

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) final;

    void SetValue( const Variable<double>& rThisVariable, const double rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final;

    /// Get the value of the hardening function w.r.t consistent parameter
    double GetValue(const double phi) const final;

    /// Get the derivative of the hardening function w.r.t consistent parameter
    double GetDerivative(const double phi) const final;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ExponentialHardeningLaw";
    }

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    double mCinf, mC0;
    double mRho;

    void save(Serializer& rSerializer) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
        rSerializer.save("mCinf", mCinf);
        rSerializer.save("mC0", mC0);
        rSerializer.save("mRho", mRho);
    }

    void load(Serializer& rSerializer) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
        rSerializer.load("mCinf", mCinf);
        rSerializer.load("mC0", mC0);
        rSerializer.load("mRho", mRho);
    }

    ///@}

}; /* Class ExponentialHardeningLaw */

//template class KRATOS_API(KRATOS_CORE) KratosComponents<ExponentialHardeningLaw>;

//void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, ExponentialHardeningLaw const& ThisComponent);

} /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APP_LINEAR_HARDENING_LAW_H_INCLUDED  defined */
