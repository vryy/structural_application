//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 1 Sep 2023 $
//   Revision:            $Revision: 1.0 $
//
//



#if !defined(KRATOS_STRUCTURAL_APP_POWER_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APP_POWER_HARDENING_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "phase_laws/hardening_law.h"


namespace Kratos
{

/**
 * Power hardening law of the form q = K*(alpha+e0)^n
 */
class PowerHardeningLaw : public HardeningLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(PowerHardeningLaw);

    typedef HardeningLaw BaseType;

    /**
     * Constructor.
     */
    PowerHardeningLaw();
    PowerHardeningLaw(const double& K, const double& e0, const double& n);

    /**
     * Destructor.
     */
    virtual ~PowerHardeningLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      PowerHardeningLaw::Pointer p_clone(new PowerHardeningLaw());
     *      return p_clone;
     */
    BaseType::Pointer Clone() const final;

    /**
     * Operations
     */
    bool Has( const Variable<double>& rThisVariable ) final;

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) final;

    double GetK() const {return mK;}
    double GetE0() const {return me0;}
    double GetN() const {return mn;}
    void SetK(const double value) {mK = value;}
    void SetE0(const double value) {me0 = value;}
    void SetN(const double value) {mn = value;}

    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final;

    /// Get the value of the hardening function w.r.t consistent parameter
    double GetValue(const double& phi) const final;

    /// Get the derivative of the hardening function w.r.t consistent parameter
    double GetDerivative(const double& phi) const final;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PowerHardeningLaw";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "K: " << mK << ", e0: " << me0 << ", n: " << mn;
    }

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    double mK, me0, mn;

    void save(Serializer& rSerializer) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
        rSerializer.save("mK", mK);
        rSerializer.save("me0", me0);
        rSerializer.save("mn", mn);
    }

    void load(Serializer& rSerializer) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
        rSerializer.load("mK", mK);
        rSerializer.load("me0", me0);
        rSerializer.load("mn", mn);
    }

    ///@}

}; /* Class PowerHardeningLaw */

//template class KRATOS_API(KRATOS_CORE) KratosComponents<PowerHardeningLaw>;

//void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, PowerHardeningLaw const& ThisComponent);

} /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APP_POWER_HARDENING_LAW_H_INCLUDED  defined */

