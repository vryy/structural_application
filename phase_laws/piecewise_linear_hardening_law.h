//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 9 Oct 2022 $
//   Revision:            $Revision: 1.0 $
//
//



#if !defined(KRATOS_STRUCTURAL_APP_PIECEWISE_LINEAR_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APP_PIECEWISE_LINEAR_HARDENING_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "phase_laws/hardening_law.h"


namespace Kratos
{

/**
 * Nonlinear hardening law characterized by piecewise curve
 */
class PiecewiseLinearHardeningLaw : public HardeningLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(PiecewiseLinearHardeningLaw);

    typedef HardeningLaw BaseType;

    /**
     * Constructor.
     */
    PiecewiseLinearHardeningLaw();

    /**
     * Destructor.
     */
    virtual ~PiecewiseLinearHardeningLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      PiecewiseLinearHardeningLaw::Pointer p_clone(new PiecewiseLinearHardeningLaw());
     *      return p_clone;
     */
    BaseType::Pointer Clone() const final;

    /**
     * Operations
     */

    void AddPoint(const double& x, const double& y)
    {
        mPx.push_back(x);
        mPy.push_back(y);
    }

    bool Has( const Variable<double>& rThisVariable ) final;

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) final;

    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final;

    /// Get the value of the hardening function w.r.t consistent parameter
    double GetValue(const double& phi) const final;

    /// Get the derivative of the hardening function w.r.t consistent parameter
    double GetDerivative(const double& phi) const final;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PiecewiseLinearHardeningLaw";
    }

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    std::vector<double> mPx;
    std::vector<double> mPy;

    void save(Serializer& rSerializer) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
        rSerializer.save("X", mPx);
        rSerializer.save("Y", mPy);
    }

    void load(Serializer& rSerializer) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
        rSerializer.load("X", mPx);
        rSerializer.load("Y", mPy);
    }

    ///@}

}; /* Class PiecewiseLinearHardeningLaw */

//template class KRATOS_API(KRATOS_CORE) KratosComponents<PiecewiseLinearHardeningLaw>;

//void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, PiecewiseLinearHardeningLaw const& ThisComponent);

} /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_APP_PIECEWISE_LINEAR_HARDENING_LAW_H_INCLUDED  defined */
