//
//   Project Name:        KratosStructuralApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2 Sep 2026 $
//
//

#if !defined(KRATOS_STRUCTURAL_APP_MIXED_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APP_MIXED_HARDENING_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "phase_laws/hardening_law.h"


namespace Kratos
{

/**
 * Mixed hardening law of the form sigmay = sigma0 + Q*(1-exp(-b*alpha)) + H*alpha
 */
class KRATOS_API(STRUCTURAL_APPLICATION) MixedHardeningLaw : public HardeningLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(MixedHardeningLaw);

    typedef HardeningLaw BaseType;

    /**
     * Constructor.
     */
    MixedHardeningLaw();
    MixedHardeningLaw(const double Oy, const double H, const double Q, const double b);

    /**
     * Destructor.
     */
    ~MixedHardeningLaw() override;

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      MixedHardeningLaw::Pointer p_clone(new MixedHardeningLaw());
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
        return "MixedHardeningLaw";
    }

private:

    ///@name Serialization
    ///@{

    friend class Serializer;

    double mOy;
    double mH;
    double mQ;
    double mb;

    void save(Serializer& rSerializer) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("mOy", mOy);
        rSerializer.save("mH", mH);
        rSerializer.save("mQ", mQ);
        rSerializer.save("mb", mb);
    }

    void load(Serializer& rSerializer) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("mOy", mOy);
        rSerializer.load("mH", mH);
        rSerializer.load("mQ", mQ);
        rSerializer.load("mb", mb);
    }

    ///@}

}; /* Class MixedHardeningLaw */

} /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APP_MIXED_HARDENING_LAW_H_INCLUDED  defined */
