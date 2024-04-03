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


#if !defined(KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED


// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "custom_python/add_constitutive_laws_to_python.h"
#include "constitutive_laws/dummy_constitutive_law.h"
#include "constitutive_laws/tutorial_damage_model.h"
#include "constitutive_laws/isotropic_2d.h"
#include "constitutive_laws/isotropic_3d.h"
#include "constitutive_laws/isotropic_3d_dc.h"
#include "constitutive_laws/neo_hookean_2d.h"
#include "constitutive_laws/neo_hookean_3d.h"
#include "constitutive_laws/hyperelastic_3d.h"
#include "constitutive_laws/hyperelastic_2d.h"
// #include "constitutive_laws/viscoelastic_2d.h" // new VISCOELASTICITY
// #include "constitutive_laws/viscofibers_2d.h" // new VISCOELASTIC Fibers
// #include "constitutive_laws/viscofibers_hypermatrix_2d.h" // new VISCOELASTIC Fibers and Hyperelastic Matrix
#include "constitutive_laws/von_mises_3d.h"
#include "constitutive_laws/hypoelastic_2d.h"
#include "constitutive_laws/plane_strain.h"
#include "constitutive_laws/plane_stress.h"
#include "constitutive_laws/fluid_2d.h"
// #include "constitutive_laws/external_isotropic_3d.h"
#include "constitutive_laws/drucker_prager.h"
#include "constitutive_laws/cam_clay_3d.h"
//#include "constitutive_laws/isotropic_elastic_large_strain.h"
// #include "constitutive_laws/hooks_law.h"
#include "constitutive_laws/isotropic_planestress_wrinkling.h"
#include "constitutive_laws/isotropic_damage_2d.h"
#include "constitutive_laws/isotropic_rankine_damage_2d.h"
#include "constitutive_laws/isotropic_rankine_damage_3d.h"
#include "constitutive_laws/isotropic_damage_3d.h"
#include "constitutive_laws/isotropic_damage_implex.h"
#include "constitutive_laws/plasticity_2d.h"
#include "constitutive_laws/plane_stress_J2.h"
#include "constitutive_laws/brittle_material_2d.h"
#include "constitutive_laws/orthotropic_3d.h"
#include "constitutive_laws/st_venant_kirchhoff.h"
#include "constitutive_laws/multiplicative_finite_strain_bridging_constitutive_law.h"
#include "constitutive_laws/multiplicative_finite_strain_bridging_constitutive_law_dc.h"
#include "constitutive_laws/total_lagrangian_bridging_constitutive_law.h"
#include "constitutive_laws/hyperelastic_finite_strain_bridging_constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "python/add_mesh_to_python.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"
#include "fluency_criteria/fluency_criteria.h"
#include "constitutive_laws/von_mises_3d.h"
#include "constitutive_laws/mohr_coulomb_plane_strain.h"
#include "constitutive_laws/values_container_constitutive_law.h"
#include "constitutive_laws/hardening_law.h"
#include "constitutive_laws/linear_hardening_law.h"
#include "constitutive_laws/exponential_hardening_law.h"
#include "constitutive_laws/piecewise_linear_hardening_law.h"
#include "constitutive_laws/power_hardening_law.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

typedef ConstitutiveLaw ConstitutiveLawBaseType;
typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
typedef FluencyCriteria::Pointer FluencyCriteriaPointer;
typedef SofteningHardeningCriteria::Pointer SofteningHardeningCriteriaPointer;
typedef Properties::Pointer PropertiesPointer;

typedef std::vector<ConstitutiveLaw::Pointer> MaterialsContainer;
typedef ConstitutiveLaw::Pointer  ConstitutiveLawPointer;

void Push_Back_Constitutive_Laws( MaterialsContainer& ThisMaterialsContainer,
                                  ConstitutiveLawPointer ThisConstitutiveLaw )
{
    ThisMaterialsContainer.push_back( ThisConstitutiveLaw );
}

void HardeningLaw_Assign(HardeningLaw& rDummy,
    Variable<HardeningLaw::Pointer>& rThisVariable, HardeningLaw::Pointer pLaw, Properties::Pointer pProperties)
{
    rDummy.Assign(rThisVariable, pLaw, pProperties);
}

template<int TStressType>
void AddMultiplicativeFiniteStrainBridgingConstitutiveLaw(const std::string& postfix)
{
    std::string name;
    name = "MultiplicativeFiniteStrainBridgingConstitutiveLaw_" + postfix;
    class_< MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType>, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    name = "MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw_" + postfix;
    class_< MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLaw<TStressType>, bases< MultiplicativeFiniteStrainBridgingConstitutiveLaw<TStressType> >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    name = "MultiplicativeFiniteStrainBridgingConstitutiveLawDC_" + postfix;
    class_< MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType>, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( name.c_str(), init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    name = "MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLawDC_" + postfix;
    class_< MultiplicativeFiniteStrainAxisymmetricBridgingConstitutiveLawDC<TStressType>, bases< MultiplicativeFiniteStrainBridgingConstitutiveLawDC<TStressType> >, boost::noncopyable >
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

    class_< Isotropic2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "Isotropic2D",
      init<>() )
    // .def("Clone",              &Isotropic2D::Clone)
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

    class_< MohrCoulombPlaneStrain, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "MohrCoulombPlaneStrain",
      init<>() )
    ;

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

    class_< DruckerPrager, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "DruckerPrager",
      init<>() )
    ;

    class_< Orthotropic3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
    ( "Orthotropic3D",
      init<>() )
    ;

    class_< Isotropic_Damage_2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "IsotropicDamage2D",
      init<>() )
    .def( init<FluencyCriteriaPointer, SofteningHardeningCriteriaPointer, PropertiesPointer>() )
    ;

    class_< Isotropic_Damage_3D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "IsotropicDamage3D",
      init<>() )
    .def( init<FluencyCriteriaPointer, SofteningHardeningCriteriaPointer, PropertiesPointer>() )
    ;

    class_< IsotropicDamageIMPLEX, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "IsotropicDamageIMPLEX",
      init<>() )
    .def( init<>() )
    ;


    class_<Plasticity2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "Plasticity2D",
      init<>() )
    .def( init<FluencyCriteriaPointer, PropertiesPointer>() )
    ;

    class_<PlaneStressJ2, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "PlaneStressJ2",
      init<>() )
    ;

//             class_<Plasticity3D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
//             ("Plasticity3D",
//             init<>() )
//                 .def(init<FluencyCriteriaPointer,SofteningHardeningCriteriaPointer, PropertiesPointer>())
//                  ;





    class_<BrittleMaterial2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "BrittleMaterial2D",
      init<>() )
    .def( init<FluencyCriteriaPointer, PropertiesPointer>() )
    ;


    class_<IsotropicRankineDamage2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "IsotropicRankineDamage2D",
      init<>() )
    .def( init<>() )
    ;

    class_<IsotropicRankineDamage3D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "IsotropicRankineDamage3D",
      init<>() )
    .def( init<>() )
    ;

    class_< VonMises3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
    ( "VonMises3D",
      init<>() )
    ;

    class_< Hypoelastic2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "Hypoelastic2D",
      init<>() )
    ;

    class_< Fluid2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "Fluid2D",
      init<>() )
    ;

    // class_< ExternalIsotropic3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
    // ( "ExternalIsotropic3D",
    //   init<>() )
    // ;

    // class_< HooksLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    // ( "HooksLaw",
    //   init<>() )
    // ;

    class_< IsotropicPlaneStressWrinkling, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "IsotropicPlaneStressWrinkling",
      init<>() )
    ;

    class_< Hyperelastic3D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "Hyperelastic3D",
      init<>() )
    ;


    class_< Hyperelastic2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "Hyperelastic2D",
      init<>() )
    ;

    class_< CamClay3D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "CamClay3D",
      init<>() )
    ;

    class_< ValuesContainerConstitutiveLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "ValuesContainerConstitutiveLaw", init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    AddMultiplicativeFiniteStrainBridgingConstitutiveLaw<1>("Cauchy");
    AddMultiplicativeFiniteStrainBridgingConstitutiveLaw<2>("Kirchhoff");

    class_< TotalLagrangianBridgingConstitutiveLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "TotalLagrangianBridgingConstitutiveLaw", init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    class_< HyperelasticFiniteStrainBridgingConstitutiveLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperelasticFiniteStrainBridgingConstitutiveLaw", init<>() )
    .def(init<ConstitutiveLawBaseType::Pointer>())
    ;

    /*
           class_< Viscofibers2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                   ("Viscofibers2D",
                    init<>() )
                   ;


           class_< Viscoelastic2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                   ("Viscoelastic2D",
                    init<>() )
                   ;


    class_< Viscofibers_Hypermatrix2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                   ("Viscofibers_Hypermatrix2D",
                    init<>() )
                   ;
     */
//    class_<Plane_Stress_Damage_Orthotropic_2D  , bases< ConstitutiveLawBaseType >, boost::noncopyable >
//    ("PlaneStressDamageOrthotropic2D",
//    init<>() )
//    //.def(init<FluencyCriteriaType const&>())
//                         .def(init<FluencyCriteriaPointer>())
//    ;
    /*
       class_<ComposeMaterial , bases< ConstitutiveLawBaseType >, boost::noncopyable >
       ("ComposeMaterial",
       init<>() )
                            .def(init<MaterialsContainer>())
       ;*/

    class_<Variable<HardeningLaw::Pointer>, bases<VariableData>, boost::noncopyable >( "HardeningLawVariable", no_init );

    double(HardeningLaw::*pointer_to_GetValue)(const double&) const = &HardeningLaw::GetValue;

    class_< HardeningLaw, bases< Flags >, boost::noncopyable >
    ( "HardeningLaw", init<>() )
    .def("GetValue", pointer_to_GetValue)
    .def("GetDerivative", &HardeningLaw::GetDerivative)
    .def("Assign", &HardeningLaw_Assign)
    ;

    class_< LinearHardeningLaw, bases< HardeningLaw >, boost::noncopyable >
    ( "LinearHardeningLaw", init<>() )
    .def(init<const double&, const double&>())
    ;

    class_< ExponentialHardeningLaw, bases< HardeningLaw >, boost::noncopyable >
    ( "ExponentialHardeningLaw", init<>() )
    .def(init<const double&, const double&, const double&>())
    ;

    class_< PiecewiseLinearHardeningLaw, bases< HardeningLaw >, boost::noncopyable >
    ( "PiecewiseLinearHardeningLaw", init<>() )
    .def("AddPoint", &PiecewiseLinearHardeningLaw::AddPoint)
    ;

    class_< PowerHardeningLaw, bases< HardeningLaw >, boost::noncopyable >
    ( "PowerHardeningLaw", init<>() )
    .def(init<const double&, const double&, const double&>())
    ;
}
}  // namespace Python.
}  // namespace Kratos.
#endif // KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED defined
