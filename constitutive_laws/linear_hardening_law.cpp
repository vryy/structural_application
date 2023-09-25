//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Dec 2018 $
//   Revision:            $Revision: 1.0 $
//
//



#include "constitutive_laws/linear_hardening_law.h"
#include "structural_application_variables.h" // from structural_application


namespace Kratos
{

LinearHardeningLaw::LinearHardeningLaw() : HardeningLaw(), mOy(0.0), mH(0.0)
{}

LinearHardeningLaw::LinearHardeningLaw(const double& Oy, const double& H) : HardeningLaw(), mOy(Oy), mH(H)
{}

LinearHardeningLaw::~LinearHardeningLaw()
{}

HardeningLaw::Pointer LinearHardeningLaw::Clone() const
{
    return HardeningLaw::Pointer(new LinearHardeningLaw(mOy, mH));
}

bool LinearHardeningLaw::Has( const Variable<double>& rThisVariable )
{
    if( rThisVariable == TENSILE_STRENGTH )
        return true;
    if( rThisVariable == COHESION )
        return true;
    if( rThisVariable == ISOTROPIC_HARDENING_MODULUS )
        return true;
    return false;
}

double& LinearHardeningLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if( rThisVariable == TENSILE_STRENGTH )
        rValue = mOy;
    if( rThisVariable == COHESION )
        rValue = mOy;
    if( rThisVariable == ISOTROPIC_HARDENING_MODULUS )
        rValue = mH;
    return rValue;
}

void LinearHardeningLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
    if( rThisVariable == TENSILE_STRENGTH )
        mOy = rValue;
    if( rThisVariable == COHESION )
        mOy = rValue;
    if( rThisVariable == ISOTROPIC_HARDENING_MODULUS )
        mH = rValue;
}

double LinearHardeningLaw::GetValue(const double& phi) const
{
    return mOy + mH*phi;
}

double LinearHardeningLaw::GetDerivative(const double& phi) const
{
    return mH;
}

} /* namespace Kratos.*/

