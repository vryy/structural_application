//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 1 Sep 2023 $
//   Revision:            $Revision: 1.0 $
//
//



#include "constitutive_laws/power_hardening_law.h"
#include "structural_application_variables.h"


namespace Kratos
{

PowerHardeningLaw::PowerHardeningLaw() : HardeningLaw(), mK(0.0), me0(1.0e-10), mn(0.0)
{}

PowerHardeningLaw::PowerHardeningLaw(const double& K, const double& e0, const double& n)
: HardeningLaw(), mK(K), me0(e0), mn(n)
{}

PowerHardeningLaw::~PowerHardeningLaw()
{}

HardeningLaw::Pointer PowerHardeningLaw::Clone() const
{
    return HardeningLaw::Pointer(new PowerHardeningLaw(mK, me0, mn));
}

bool PowerHardeningLaw::Has( const Variable<double>& rThisVariable )
{
}

double& PowerHardeningLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
}

void PowerHardeningLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

double PowerHardeningLaw::GetValue(const double& phi) const
{
    return mK*pow(phi + me0, mn);
}

double PowerHardeningLaw::GetDerivative(const double& phi) const
{
    return mK*mn*pow(phi + me0, mn-1);
}

} /* namespace Kratos.*/
