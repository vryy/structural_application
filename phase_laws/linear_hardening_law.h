//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Dec 2018 $
//   Revision:            $Revision: 1.0 $
//
//



#if !defined(KRATOS_STRUCTURAL_APP_LINEAR_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APP_LINEAR_HARDENING_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "phase_laws/hardening_law.h"


namespace Kratos
{

/**
 * Linear hardening law of the form q = oy + H*alpha
 */
class LinearHardeningLaw : public HardeningLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(LinearHardeningLaw);

    typedef HardeningLaw BaseType;

    /**
     * Constructor.
     */
    LinearHardeningLaw();
    LinearHardeningLaw(const double Oy, const double H);

    /**
     * Destructor.
     */
    virtual ~LinearHardeningLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      LinearHardeningLaw::Pointer p_clone(new LinearHardeningLaw());
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
        return "LinearHardeningLaw";
    }

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    double mOy;
    double mH;

    void save(Serializer& rSerializer) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
        rSerializer.save("mOy", mOy);
        rSerializer.save("mH", mH);
    }

    void load(Serializer& rSerializer) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
        rSerializer.load("mOy", mOy);
        rSerializer.load("mH", mH);
    }

    ///@}

}; /* Class LinearHardeningLaw */

//template class KRATOS_API(KRATOS_CORE) KratosComponents<LinearHardeningLaw>;

//void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, LinearHardeningLaw const& ThisComponent);

} /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APP_LINEAR_HARDENING_LAW_H_INCLUDED  defined */

