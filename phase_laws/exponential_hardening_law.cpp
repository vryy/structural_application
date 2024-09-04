//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Dec 2018 $
//   Revision:            $Revision: 1.0 $
//
//



#include "phase_laws/exponential_hardening_law.h"
#include "structural_application_variables.h"


namespace Kratos
{

ExponentialHardeningLaw::ExponentialHardeningLaw() : HardeningLaw(), mCinf(0.0), mC0(0.0), mRho(0.0)
{}

ExponentialHardeningLaw::ExponentialHardeningLaw(const double C0, const double Cinf, const double Rho)
: HardeningLaw(), mCinf(Cinf), mC0(C0), mRho(Rho)
{}

ExponentialHardeningLaw::~ExponentialHardeningLaw()
{}

HardeningLaw::Pointer ExponentialHardeningLaw::Clone() const
{
    return HardeningLaw::Pointer(new ExponentialHardeningLaw(mCinf, mC0, mRho));
}

bool ExponentialHardeningLaw::Has( const Variable<double>& rThisVariable )
{
}

double& ExponentialHardeningLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
}

void ExponentialHardeningLaw::SetValue( const Variable<double>& rThisVariable, const double rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

double ExponentialHardeningLaw::GetValue(const double phi) const
{
    return mCinf - (mCinf-mC0)*std::exp(-mRho*phi);
}

double ExponentialHardeningLaw::GetDerivative(const double phi) const
{
    return (mCinf-mC0)*mRho*std::exp(-mRho*phi);
}

} /* namespace Kratos.*/

