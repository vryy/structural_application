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


// Project includes
#include "custom_python3/add_constitutive_laws_to_python.h"
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "constitutive_laws/dummy_constitutive_law.h"
// #include "constitutive_laws/tutorial_damage_model.h"
// #include "constitutive_laws/isotropic_2d.h"
#include "constitutive_laws/isotropic_3d.h"
// #include "constitutive_laws/neo_hookean_2d.h"
// #include "constitutive_laws/neo_hookean_3d.h"
// #include "constitutive_laws/hyperelastic_3d.h"
// #include "constitutive_laws/hyperelastic_2d.h"
// #include "constitutive_laws/viscoelastic_2d.h" // new VISCOELASTICITY
// #include "constitutive_laws/viscofibers_2d.h" // new VISCOELASTIC Fibers
// #include "constitutive_laws/viscofibers_hypermatrix_2d.h" // new VISCOELASTIC Fibers and Hyperelastic Matrix
// #include "constitutive_laws/von_mises_3d.h"
// #include "constitutive_laws/hypoelastic_2d.h"
#include "constitutive_laws/plane_strain.h"
#include "constitutive_laws/plane_stress.h"
// #include "constitutive_laws/fluid_2d.h"
// #include "constitutive_laws/external_isotropic_3d.h"
// #include "constitutive_laws/drucker_prager.h"
// #include "constitutive_laws/cam_clay_3d.h"
//#include "constitutive_laws/isotropic_elastic_large_strain.h"
// #include "constitutive_laws/hooks_law.h"
// #include "constitutive_laws/isotropic_planestress_wrinkling.h"
// #include "constitutive_laws/isotropic_damage_2d.h"
// #include "constitutive_laws/isotropic_rankine_damage_2d.h"
// #include "constitutive_laws/isotropic_rankine_damage_3d.h"
// #include "constitutive_laws/isotropic_damage_3d.h"
// #include "constitutive_laws/isotropic_damage_implex.h"
// #include "constitutive_laws/plasticity_2d.h"
// #include "constitutive_laws/plane_stress_J2.h"
// #include "constitutive_laws/brittle_material_2d.h"
// #include "constitutive_laws/orthotropic_3d.h"
// #include "includes/node.h"
// #include "includes/variables.h"
// #include "includes/mesh.h"
// #include "includes/element.h"
// #include "includes/condition.h"
// #include "includes/properties.h"
// #include "python/add_mesh_to_python.h"
// #include "python/pointer_vector_set_python_interface.h"
// #include "python/variable_indexing_python.h"
// #include "fluency_criteria/fluency_criteria.h"
// #include "constitutive_laws/von_mises_3d.h"
// #include "constitutive_laws/mohr_coulomb_plane_strain.h"
// #include "constitutive_laws/values_container_constitutive_law.h"
#include "constitutive_laws/hardening_law.h"
#include "constitutive_laws/linear_hardening_law.h"
#include "constitutive_laws/exponential_hardening_law.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

typedef ConstitutiveLaw ConstitutiveLawBaseType;
// typedef FluencyCriteria::Pointer FluencyCriteriaPointer;
// typedef SofteningHardeningCriteria::Pointer SofteningHardeningCriteriaPointer;
// typedef Properties::Pointer PropertiesPointer;

typedef std::vector<ConstitutiveLaw::Pointer> MaterialsContainer;
typedef ConstitutiveLaw::Pointer  ConstitutiveLawPointer;

void Push_Back_Constitutive_Laws( MaterialsContainer& ThisMaterialsContainer,
                                  ConstitutiveLawPointer ThisConstitutiveLaw )
{
    ThisMaterialsContainer.push_back( ThisConstitutiveLaw );
}

double HardeningLaw_GetValue(HardeningLaw& rDummy, const double& phi)
{
    return rDummy.GetValue(phi);
}

double HardeningLaw_GetDerivative(HardeningLaw& rDummy, const double& phi)
{
    return rDummy.GetDerivative(phi);
}

void StructuralApplication_AddConstitutiveLawsToPython(pybind11::module& m)
{
    // class_< MaterialsContainer >( m, "MaterialsContainer" )
    // .def( init<>() )
    // .def( "PushBack", Push_Back_Constitutive_Laws )
    // ;

    class_< DummyConstitutiveLaw, DummyConstitutiveLaw::Pointer, ConstitutiveLawBaseType >
    ( m, "DummyConstitutiveLaw" )
    .def( init<>() )
    ;

    // class_< TutorialDamageModel, ConstitutiveLawBaseType >
    // ( m, "TutorialDamageModel" )
    // .def( init<>() )
    // ;

    // class_< Isotropic2D, ConstitutiveLawBaseType >
    // ( m, "Isotropic2D" )
    // .def( init<>() )
    // // .def("Clone",              &Isotropic2D::Clone)
    // ;

    class_< PlaneStrain, PlaneStrain::Pointer, ConstitutiveLawBaseType >
    ( m, "PlaneStrain" )
    .def( init<>() )
    ;

    class_< PlaneStress, PlaneStress::Pointer, ConstitutiveLawBaseType >
    ( m, "PlaneStress" )
    .def( init<>() )
    ;

    // class_< NeoHookean2D, ConstitutiveLawBaseType >
    // ( m, "NeoHookean2D" )
    // init<>() );

    // class_< NeoHookean3D, ConstitutiveLawBaseType >
    // ( m, "NeoHookean3D" )
    // init<>() );

    // class_< MohrCoulombPlaneStrain, ConstitutiveLawBaseType >
    // ( m, "MohrCoulombPlaneStrain" )
    // .def( init<>() )
    // ;

    class_< Isotropic3D, Isotropic3D::Pointer, ConstitutiveLawBaseType >
    ( m, "Isotropic3D" )
    .def( init<>() )
    ;

    // class_< DruckerPrager, ConstitutiveLawBaseType >
    // ( m, "DruckerPrager" )
    // .def( init<>() )
    // ;

    // class_< Orthotropic3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
    // ( m, "Orthotropic3D" )
    // .def( init<>() )
    // ;

    // class_< Isotropic_Damage_2D, ConstitutiveLawBaseType >
    // ( m, "IsotropicDamage2D" )
    // .def( init<>() )
    // .def( init<FluencyCriteriaPointer, SofteningHardeningCriteriaPointer, PropertiesPointer>() )
    // ;

    // class_< Isotropic_Damage_3D, ConstitutiveLawBaseType >
    // ( m, "IsotropicDamage3D" )
    // .def( init<>() )
    // .def( init<FluencyCriteriaPointer, SofteningHardeningCriteriaPointer, PropertiesPointer>() )
    // ;

    // class_< IsotropicDamageIMPLEX, ConstitutiveLawBaseType >
    // ( m, "IsotropicDamageIMPLEX" )
    // .def( init<>() )
    // ;


    // class_<Plasticity2D, ConstitutiveLawBaseType >
    // ( m, "Plasticity2D" )
    // .def( init<>() )
    // .def( init<FluencyCriteriaPointer, PropertiesPointer>() )
    // ;

    // class_<PlaneStressJ2, ConstitutiveLawBaseType >
    // ( m, "PlaneStressJ2" )
    // .def( init<>() )
    ;

//             class_<Plasticity3D, ConstitutiveLawBaseType >
//             ("Plasticity3D")
//             .def( init<>() )
//             .def(init<FluencyCriteriaPointer,SofteningHardeningCriteriaPointer, PropertiesPointer>())
//             ;





    // class_<BrittleMaterial2D, ConstitutiveLawBaseType >
    // ( m, "BrittleMaterial2D" )
    // .def( init<>() )
    // .def( init<FluencyCriteriaPointer, PropertiesPointer>() )
    // ;


    // class_<IsotropicRankineDamage2D, ConstitutiveLawBaseType >
    // ( m, "IsotropicRankineDamage2D" )
    // .def( init<>() )
    // ;

    // class_<IsotropicRankineDamage3D, ConstitutiveLawBaseType >
    // ( m, "IsotropicRankineDamage3D" )
    // .def( init<>() )
    // ;

    // class_< VonMises3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
    // ( m, "VonMises3D" )
    // .def( init<>() )
    // ;

    // class_< Hypoelastic2D, ConstitutiveLawBaseType >
    // ( m, "Hypoelastic2D" )
    // .def( init<>() )
    // ;

    // class_< Fluid2D, ConstitutiveLawBaseType >
    // ( m, "Fluid2D" )
    // .def( init<>() )
    // ;

    // class_< ExternalIsotropic3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
    // ( m, "ExternalIsotropic3D" )
    // .def( init<>() )
    // ;

    // class_< HooksLaw, ConstitutiveLawBaseType >
    // ( m, "HooksLaw" )
    // .def( init<>() )
    // ;

    // class_< IsotropicPlaneStressWrinkling, ConstitutiveLawBaseType >
    // ( m, "IsotropicPlaneStressWrinkling" )
    // .def( init<>() )
    // ;

    // class_< Hyperelastic3D, ConstitutiveLawBaseType >
    // ( m, "Hyperelastic3D" )
    // .def( init<>() )
    // ;


    // class_< Hyperelastic2D, ConstitutiveLawBaseType >
    // ( m, "Hyperelastic2D" )
    // .def( init<>() )
    // ;

    // class_< CamClay3D, ConstitutiveLawBaseType >
    // ( m, "CamClay3D" )
    // .def( init<>() )
    // ;

    // class_< ValuesContainerConstitutiveLaw, ConstitutiveLawBaseType >
    // ( m, "ValuesContainerConstitutiveLaw")
    // .def( init<>() )
    // .def(init<ConstitutiveLawBaseType::Pointer>())
    // ;

    /*
           class_< Viscofibers2D, ConstitutiveLawBaseType >
                   (m, "Viscofibers2D")
                  .def( init<>() )
                   ;


           class_< Viscoelastic2D, ConstitutiveLawBaseType >
                   (m, "Viscoelastic2D")
                  .def( init<>() )
                   ;


    class_< Viscofibers_Hypermatrix2D, ConstitutiveLawBaseType >
                   (m, "Viscofibers_Hypermatrix2D")
                  .def( init<>() )
                   ;
     */
//    class_<Plane_Stress_Damage_Orthotropic_2D  , ConstitutiveLawBaseType >
//    (m, "PlaneStressDamageOrthotropic2D")
//    .def( init<>() )
//    //.def(init<FluencyCriteriaType const&>())
//                         .def(init<FluencyCriteriaPointer>())
//    ;
    /*
       class_<ComposeMaterial , ConstitutiveLawBaseType >
       (m, "ComposeMaterial")
       .def( init<>() )
                            .def(init<MaterialsContainer>())
       ;*/

    class_<Variable<HardeningLaw::Pointer>, VariableData>(m, "HardeningLawVariable");

    class_<HardeningLaw, HardeningLaw::Pointer, Flags>
    (m, "HardeningLaw")
    .def(init<>())
    .def("GetValue", &HardeningLaw_GetValue)
    .def("GetDerivative", &HardeningLaw_GetDerivative)
    ;

    class_<LinearHardeningLaw, LinearHardeningLaw::Pointer, HardeningLaw>
    (m, "LinearHardeningLaw")
    .def(init<>())
    .def(init<const double&, const double&>())
    ;

    class_<ExponentialHardeningLaw, ExponentialHardeningLaw::Pointer, HardeningLaw >
    (m, "ExponentialHardeningLaw")
    .def(init<>())
    .def(init<const double&, const double&, const double&>())
    ;

}
}  // namespace Python.
}  // namespace Kratos.
