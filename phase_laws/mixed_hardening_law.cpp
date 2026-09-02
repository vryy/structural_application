//
//   Project Name:        KratosStructuralApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2 Sep 2026 $
//
//

#include "phase_laws/mixed_hardening_law.h"
#include "structural_application_variables.h"

namespace Kratos
{

MixedHardeningLaw::MixedHardeningLaw() : HardeningLaw(), mOy(0.0), mH(0.0), mQ(0.0), mb(0.0)
{}

MixedHardeningLaw::MixedHardeningLaw(const double Oy, const double H, const double Q, const double b)
: HardeningLaw(), mOy(Oy), mH(H), mQ(Q), mb(b)
{}

MixedHardeningLaw::~MixedHardeningLaw()
{}

HardeningLaw::Pointer MixedHardeningLaw::Clone() const
{
    return HardeningLaw::Pointer(new MixedHardeningLaw(mOy, mH, mQ, mb));
}

bool MixedHardeningLaw::Has( const Variable<double>& rThisVariable )
{
    if( rThisVariable == TENSILE_STRENGTH )
        return true;
    if( rThisVariable == COHESION )
        return true;
    if( rThisVariable == ISOTROPIC_HARDENING_MODULUS )
        return true;
    return false;
}

double& MixedHardeningLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if( rThisVariable == TENSILE_STRENGTH )
        rValue = mOy;
    if( rThisVariable == COHESION )
        rValue = mOy;
    if( rThisVariable == ISOTROPIC_HARDENING_MODULUS )
        rValue = mH;
    return rValue;
}

void MixedHardeningLaw::SetValue( const Variable<double>& rThisVariable, const double rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == TENSILE_STRENGTH )
        mOy = rValue;
    if( rThisVariable == COHESION )
        mOy = rValue;
    if( rThisVariable == ISOTROPIC_HARDENING_MODULUS )
        mH = rValue;
}

double MixedHardeningLaw::GetValue(const double phi) const
{
    return mOy + mH*phi + mQ*(1.0 - std::exp(-mb*phi));
}

double MixedHardeningLaw::GetDerivative(const double phi) const
{
    return mH + mQ*mb*std::exp(-mb*phi);
}

} /* namespace Kratos.*/
