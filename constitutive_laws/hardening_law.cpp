//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Dec 2018 $
//   Revision:            $Revision: 1.0 $
//
//



#include "constitutive_laws/hardening_law.h"


namespace Kratos
{

HardeningLaw::HardeningLaw()
{}

HardeningLaw::~HardeningLaw()
{}

HardeningLaw::Pointer HardeningLaw::Clone() const
{
    return HardeningLaw::Pointer(new HardeningLaw());
}

bool HardeningLaw::Has( const Variable<int>& rThisVariable )
{
}

bool HardeningLaw::Has( const Variable<double>& rThisVariable )
{
}

bool HardeningLaw::Has( const Variable<Vector>& rThisVariable )
{
}

bool HardeningLaw::Has( const Variable<Matrix>& rThisVariable )
{
}

int& HardeningLaw::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
}

double& HardeningLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
}

Vector& HardeningLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
}

Matrix& HardeningLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
}

void HardeningLaw::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

void HardeningLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

void HardeningLaw::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                           const array_1d<double, 3>& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

void HardeningLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

void HardeningLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                           const ProcessInfo& rCurrentProcessInfo )
{
}

double HardeningLaw::GetValue(const double& phi) const
{
    return 0.0;
}

double HardeningLaw::GetDerivative(const double& phi) const
{
    return 0.0;
}

} /* namespace Kratos.*/

