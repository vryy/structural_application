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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"
#include "custom_python/add_constitutive_laws_to_python.h"
#include "constitutive_laws/dummy_constitutive_law.h"
#include "constitutive_laws/tutorial_damage_model.h"
#include "constitutive_laws/isotropic_3d.h"
#include "constitutive_laws/isotropic_3d_dc.h"
#include "constitutive_laws/neo_hookean_2d.h"
#include "constitutive_laws/neo_hookean_3d.h"
#include "constitutive_laws/plane_strain.h"
#include "constitutive_laws/plane_stress.h"
#include "constitutive_laws/cam_clay_3d.h"
#include "constitutive_laws/isotropic_damage_implex.h"
#include "constitutive_laws/st_venant_kirchhoff.h"
#include "constitutive_laws/multiplicative_finite_strain_bridging_constitutive_law.h"
#include "constitutive_laws/multiplicative_finite_strain_bridging_constitutive_law_dc.h"
#include "constitutive_laws/total_lagrangian_bridging_constitutive_law.h"
#include "constitutive_laws/hyperelastic_finite_strain_bridging_constitutive_law.h"
#include "constitutive_laws/hypoelastic_finite_strain_bridging_constitutive_law.h"
#include "constitutive_laws/hypoelastic_finite_strain_bridging_constitutive_law_dc.h"
#include "constitutive_laws/values_container_constitutive_law.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

typedef ConstitutiveLaw ConstitutiveLawBaseType;
typedef Properties::Pointer PropertiesPointer;

typedef std::vector<ConstitutiveLaw::Pointer> MaterialsContainer;
typedef ConstitutiveLaw::Pointer  ConstitutiveLawPointer;

void Push_Back_Constitutive_Laws( MaterialsContainer& ThisMaterialsContainer,
                                  ConstitutiveLawPointer ThisConstitutiveLaw )
{
    ThisMaterialsContainer.push_back( ThisConstitutiveLaw );
}

template<int TStressType>
void AddMultiplicativeFiniteStrainBridgingConstitutiveLaw(const std::string& postfix)
{
    typedef FiniteStrainBridgingConstitutiveLaw BaseType;

    std::string name;
    name = "MultiplicativeFiniteStrainBridgingConstitutiveLaw_" + postfix;
    class_< MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>, bases< BaseType >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    name = "MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw_" + postfix;
    class_< MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw<TStressType>, bases< MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType> >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    name = "MultiplicativeFiniteStrainBridgingConstitutiveLawDC_" + postfix;
    class_< MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>, bases< BaseType >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    name = "MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLawDC_" + postfix;
    class_< MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLawDC<TStressType>, bases< MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType> >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;
}

template<int THWSchemeType, int TStressType>
void AddHypoelasticFiniteStrainBridgingConstitutiveLaw(const std::string& postfix)
{
    typedef FiniteStrainBridgingConstitutiveLaw BaseType;

    std::string name;
    name = "HypoelasticFiniteStrainBridgingConstitutiveLaw_" + postfix;
    class_< HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType>, bases< BaseType >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    name = "HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw_" + postfix;
    class_< HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLaw<THWSchemeType, TStressType>, bases< HypoelasticFiniteStrainBridgingConstitutiveLaw<THWSchemeType, TStressType> >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    name = "HypoelasticFiniteStrainBridgingConstitutiveLawDC_" + postfix;
    class_< HypoelasticFiniteStrainBridgingConstitutiveLawDC<THWSchemeType, TStressType>, bases< BaseType >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    name = "HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC_" + postfix;
    class_< HypoelasticFiniteStrainAxisymmetricBridgingConstitutiveLawDC<THWSchemeType, TStressType>, bases< HypoelasticFiniteStrainBridgingConstitutiveLawDC<THWSchemeType, TStressType> >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;
}

void  AddConstitutiveLawsToPython()
{
    class_< MaterialsContainer >( "MaterialsContainer", init<>() )
    .def( "PushBack", Push_Back_Constitutive_Laws )
    ;

    class_< DummyConstitutiveLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "DummyConstitutiveLaw",
      init<>() )
    ;

    class_< TutorialDamageModel, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "TutorialDamageModel",
      init<>() )
    ;

    void(PlaneStrain::*PlaneStrain_CalculateStress)(const double&, const double&, const Vector&, Vector&) const = &PlaneStrain::CalculateStress;
    class_< PlaneStrain, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "PlaneStrain",
      init<>() )
    .def("CalculateStress", PlaneStrain_CalculateStress)
    ;

    class_< PlaneStress, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "PlaneStress",
      init<>() )
    ;

    class_< NeoHookean2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "NeoHookean2D", init<>() );

    class_< NeoHookean3D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "NeoHookean3D", init<>() );

    void(Isotropic3D::*Isotropic3D_CalculateStress)(const double&, const double&, const Vector&, Vector&) const = &Isotropic3D::CalculateStress;
    class_< Isotropic3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
    ( "Isotropic3D",
      init<>() )
    .def("CalculateStress", Isotropic3D_CalculateStress)
    ;

    class_< Isotropic3DDC, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
    ( "Isotropic3DDC",
      init<>() )
    ;

    class_< StVenantKirchhoff<2>, bases< PlaneStrain >,  boost::noncopyable >
    ( "StVenantKirchhoff_PlaneStrain",
      init<>() )
    ;

    class_< StVenantKirchhoff<3>, bases< Isotropic3D >,  boost::noncopyable >
    ( "StVenantKirchhoff_3D",
      init<>() )
    ;

    class_< IsotropicDamageIMPLEX, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "IsotropicDamageIMPLEX",
      init<>() )
    .def( init<>() )
    ;

    class_< CamClay3D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "CamClay3D",
      init<>() )
    ;

    class_< ValuesContainerConstitutiveLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "ValuesContainerConstitutiveLaw", init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    class_< FiniteStrainBridgingConstitutiveLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "FiniteStrainBridgingConstitutiveLaw", init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    AddMultiplicativeFiniteStrainBridgingConstitutiveLaw<1>("Cauchy");
    AddMultiplicativeFiniteStrainBridgingConstitutiveLaw<2>("Kirchhoff");
    AddHypoelasticFiniteStrainBridgingConstitutiveLaw<1, 1>("Cauchy_HW");
    AddHypoelasticFiniteStrainBridgingConstitutiveLaw<1, 2>("Kirchhoff_HW");
    AddHypoelasticFiniteStrainBridgingConstitutiveLaw<2, 1>("Cauchy");
    AddHypoelasticFiniteStrainBridgingConstitutiveLaw<2, 2>("Kirchhoff");

    class_< TotalLagrangianBridgingConstitutiveLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "TotalLagrangianBridgingConstitutiveLaw", init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    class_< HyperelasticFiniteStrainBridgingConstitutiveLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperelasticFiniteStrainBridgingConstitutiveLaw", init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;
}

}  // namespace Python.

}  // namespace Kratos.
