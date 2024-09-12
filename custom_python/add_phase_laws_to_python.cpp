/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hurga $
//   Date:                $Date: 2009-02-02 14:03:23 $
//   Revision:            $Revision: 1.5 $
//
//


// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "custom_python/add_phase_laws_to_python.h"
#include "phase_laws/hardening_law.h"
#include "phase_laws/linear_hardening_law.h"
#include "phase_laws/exponential_hardening_law.h"
#include "phase_laws/piecewise_linear_hardening_law.h"
#include "phase_laws/power_hardening_law.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void HardeningLaw_Assign(HardeningLaw& rDummy,
    Variable<HardeningLaw::Pointer>& rThisVariable, HardeningLaw::Pointer pLaw, Properties::Pointer pProperties)
{
    rDummy.Assign(rThisVariable, pLaw, pProperties);
}

void AddPhaseLawsToPython()
{
    class_<Variable<HardeningLaw::Pointer>, bases<VariableData>, boost::noncopyable >( "HardeningLawVariable", no_init );

    double(HardeningLaw::*pointer_to_GetValue)(const double) const = &HardeningLaw::GetValue;
    double(HardeningLaw::*pointer_to_GetDerivative)(const double) const = &HardeningLaw::GetDerivative;

    class_< HardeningLaw, bases< Flags >, boost::noncopyable >
    ( "HardeningLaw", init<>() )
    .def("GetValue", pointer_to_GetValue)
    .def("GetDerivative", pointer_to_GetDerivative)
    .def("Assign", &HardeningLaw_Assign)
    ;

    class_< LinearHardeningLaw, bases< HardeningLaw >, boost::noncopyable >
    ( "LinearHardeningLaw", init<>() )
    .def(init<const double, const double>())
    ;

    class_< ExponentialHardeningLaw, bases< HardeningLaw >, boost::noncopyable >
    ( "ExponentialHardeningLaw", init<>() )
    .def(init<const double, const double, const double>())
    ;

    class_< PiecewiseLinearHardeningLaw, bases< HardeningLaw >, boost::noncopyable >
    ( "PiecewiseLinearHardeningLaw", init<>() )
    .def("AddPoint", &PiecewiseLinearHardeningLaw::AddPoint)
    ;

    class_< PowerHardeningLaw, bases< HardeningLaw >, boost::noncopyable >
    ( "PowerHardeningLaw", init<>() )
    .def(init<const double, const double, const double>())
    .add_property("K", &PowerHardeningLaw::GetK, &PowerHardeningLaw::SetK)
    .add_property("e0", &PowerHardeningLaw::GetE0, &PowerHardeningLaw::SetE0)
    .add_property("n", &PowerHardeningLaw::GetN, &PowerHardeningLaw::SetN)
    ;
}

}  // namespace Python.

}  // namespace Kratos.
