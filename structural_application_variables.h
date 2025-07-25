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
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-03-20 08:55:34 $
//   Revision:            $Revision: 1.20 $
//
//


#if !defined(KRATOS_STRUCTURAL_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APPLICATION_VARIABLES_H_INCLUDED

#ifdef SD_APP_FORWARD_COMPATIBILITY
    #include "custom_python3/legacy_structural_app_vars.h"
#else
    #include "includes/legacy_structural_app_vars.h"
    #include "includes/deprecated_variables.h"
#endif

#include "includes/ublas_interface.h"
#include "includes/constitutive_law.h"

#include "phase_laws/hardening_law.h"

namespace Kratos
{

#ifdef SD_APP_FORWARD_COMPATIBILITY
KRATOS_DEFINE_APPLICATION_VARIABLE(STRUCTURAL_APPLICATION, double, ALPHA)
KRATOS_DEFINE_APPLICATION_VARIABLE(STRUCTURAL_APPLICATION, double, KAPPA)
#endif

// Variables definition
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, WRINKLING_APPROACH )
//  KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, GREEN_LAGRANGE_STRAIN_TENSOR )
//  KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, PK2_STRESS_TENSOR )
//  KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, AUXILIARY_MATRIX_1 )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, YOUNG_MODULUS )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, POISSON_RATIO )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, MU )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, RETRACTION_TIME )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, THICKNESS )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, NEGATIVE_FACE_PRESSURE )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, POSITIVE_FACE_PRESSURE )

typedef Matrix fix_matrix_33;
//typedef boost::numeric::ublas::bounded_matrix<double,3,3> fix_matrix_33;
typedef Vector array3;
//typedef array_1d<double,3> array3;

KRATOS_DEFINE_APPLICATION_CONSTITUTIVE_LAW_VARIABLE( STRUCTURAL_APPLICATION, CONSTITUTIVE_LAW_NO_INITIALIZE )

KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, fix_matrix_33 , MATRIX_A )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, fix_matrix_33 , MATRIX_B )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, fix_matrix_33 , MATRIX_D )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, array3, COMPOSITE_DIRECTION )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, array3, ORTHOTROPIC_YOUNG_MODULUS )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, array3, ORTHOTROPIC_SHEAR_MODULUS )
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, ORTHOTROPIC_POISSON_RATIO )
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION , GEOMETRIC_STIFFNESS )
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION , MATERIAL_DIRECTION )
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION , JOINT_STIFFNESS )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, LINING_JOINT_STIFFNESS )
//CONTACT_LINK_MASTER is defined in condition.h
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, Condition::Pointer, CONTACT_LINK_MASTER )
//CONTACT_LINK_SLAVE is defined in condition.h
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, Condition::Pointer, CONTACT_LINK_SLAVE )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, Node<3>::Pointer,  NEAR_NODE)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, Point<3>, MASTER_CONTACT_LOCAL_POINT )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, Point<3>, MASTER_CONTACT_CURRENT_LOCAL_POINT )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, Point<3>, MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, Point<3>, SLAVE_CONTACT_LOCAL_POINT )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, Point<3>, MASTER_CONTACT_GLOBAL_POINT )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, Point<3>, MASTER_CONTACT_CURRENT_GLOBAL_POINT )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, Point<3>, SLAVE_CONTACT_GLOBAL_POINT )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, INSITU_STRESS_SCALE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, LAGRANGE_SCALE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, REFERENCE_PRESSURE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, REFERENCE_WATER_PRESSURE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, WATER_PRESSURE_SCALE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, OVERCONSOLIDATION_RATIO )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, EXCESS_PORE_WATER_PRESSURE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PRESSURE_P)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PRESSURE_Q)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, COORDINATES)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, FLUID_FLOWS)

KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, CONTACT_PENETRATION )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DAMAGE_E0 )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DAMAGE_EF )

KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, BASE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, HEIGHT )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, CROSS_AREA)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AREA )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AREA_X )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AREA_Y )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AREA_Z )
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, LOCAL_INERTIA)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, INERTIA_X)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, INERTIA_Y)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, INERTIA_Z)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, FC)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, FT)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, CONCRETE_YOUNG_MODULUS_C)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, CONCRETE_YOUNG_MODULUS_T)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, FRACTURE_ENERGY)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, CRUSHING_ENERGY)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ELASTIC_ENERGY)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PLASTIC_ENERGY)

//       KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, YIELD_STRESS)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PLASTIC_MODULUS)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PLASTICITY_INDICATOR)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ISOTROPIC_HARDENING_MODULUS)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, KINEMATIC_HARDENING_MODULUS)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, ISOTROPIC_HARDENING_TYPE)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, KINEMATIC_HARDENING_TYPE)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, DRUCKER_PRAGER_MATCHING_TYPE)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, HARDENING_POINTS_ON_CURVE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, LAMNDA) // Load factor
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DAMAGE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DAMAGE_TENSION )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DAMAGE_COMPRESSION )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ORTHOTROPIC_ANGLE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, VOLUMEN_FRACTION)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, MAX_INTERNAL_FRICTION_ANGLE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DILATANCY_ANGLE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, MAX_DILATANCY_ANGLE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, COHESION)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ISOTROPIC_ELASTIC_LIMIT)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, ORTHOTROPIC_ELASTIC_LIMIT)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, VECTOR_DAMAGE)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, ORTHOTROPIC_YOUNG_MODULUS_2D) // [E1 E2 G12]
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, ORTHOTROPIC_POISSON_RATIO_2D) // [v12 v21]
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, ELASTIC_STRAIN_VECTOR)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, PLASTIC_STRAIN_VECTOR)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, CURRENT_STRAIN_VECTOR)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, CURRENT_DEFORMATION_GRADIENT)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, CURRENT_DEFORMATION_GRADIENT_DETERMINANT)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, INTEGRATION_POINT_STRAIN_VECTOR)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, EQUIVALENT_VOLUMETRIC_ELASTIC_STRAIN)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, EQUIVALENT_VOLUMETRIC_PLASTIC_STRAIN)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, EQUIVALENT_DEVIATORIC_ELASTIC_STRAIN)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, EQUIVALENT_DEVIATORIC_PLASTIC_STRAIN)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, EQUIVALENT_VOLUMETRIC_STRAIN)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, EQUIVALENT_DEVIATORIC_STRAIN)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, EQUIVALENT_STRAIN)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, ALMANSI_PLASTIC_STRAIN)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, ALMANSI_ELASTIC_STRAIN)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, NODAL_STRESS)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, NODAL_STRESS_VECTOR)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, NODAL_STRAIN)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, NODAL_STRAIN_VECTOR)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, CONSTRAINT_MATRIX)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, PRESTRAIN)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, PRESTRESS)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PRESTRESS_ZZ)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PRESTRESS_FACTOR )
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, INITIAL_STRESS)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, CONSTRAINT_VECTOR)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DISIPATION)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int,     NODAL_VALUES)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION,  NODAL_DAMAGE)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, IS_TARGET)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, IS_CONTACTOR)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DAMPING_RATIO)
//KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, KINETIC_ENERGY)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, POTENCIAL_ENERGY)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DEFORMATION_ENERGY)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, VON_MISES_STRESS)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, YIELD_STATE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, YIELD_SURFACE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, RHS_PRESSURE)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, COMPUTE_TANGENT_MATRIX)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, USE_NUMERICAL_TANGENT)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, IS_DISCRETE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, TENSILE_STRENGTH)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, SHEAR_STRENGTH)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, VISCOUS_DAMPING)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, MAX_FREQUENCY)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, JOINT_FORCE_REACTION)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, JOINT_MOMENT_REACTION)
//KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, INTERNAL_FORCE) //already put on variables.h (warning was appearing on Windows)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, ELASTIC_BEDDING_STIFFNESS)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, IS_BBAR)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, IS_FBAR)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, FBAR_MODE)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, NEIGHBOUR_EXPANSION_LEVEL)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, STRESS_RECOVERY_TYPE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, RAYLEIGH_DAMPING_ALPHA)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, RAYLEIGH_DAMPING_BETA)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, STABILISATION_FACTOR)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, SHEAR_MODULUS)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, SHEAR_MODULUS_EVOLUTION)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, RECOVERY_STRESSES)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, THREED_STRESSES)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, THREED_PRESTRESS)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, STRESSES_OLD)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, STRAIN_OLD)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, THREED_STRAIN)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, PRE_STRAIN_VECTOR)
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, POST_STRAIN_VECTOR)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, STRESS_TENSOR)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, STRAIN_TENSOR)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, ELASTIC_STRAIN_TENSOR)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, ELASTIC_STRAIN_TENSOR_OLD)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, PLASTIC_STRAIN_TENSOR)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, LEFT_STRETCH_TENSOR)
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, RIGHT_STRETCH_TENSOR)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, PRESCRIBED_DELTA_DISPLACEMENT)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, PRESCRIBED_DISPLACEMENT)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, INITIAL_DISPLACEMENT)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, LINE_LOAD)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool,   IS_CONTACT_NODE)
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, LAGRANGE_MULTIPLIER)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, LAGRANGE_MULTIPLIER_INDEX)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, HAS_STRAIN_AT_NODE)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, HAS_STRESSES_AT_NODE)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, HAS_NODAL_ERROR )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, FORCE_EQUAL_ORDER_INTERPOLATION )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, NODAL_ERROR_1 )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DUMMY_DOF )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PRECONSOLIDATION_PRESSURE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PRECONSOLIDATION_PRESSURE_MIN )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PRECONSOLIDATION_PRESSURE_DEF )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, CSL_SLOPE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, VIRGIN_COMPRESSION_INDEX )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, RECOMPRESSION_INDEX )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, SWELL_INDEX )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, VOID_RATIO )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, SPACING_RATIO )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ASSOCIATIVITY )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, SHAPE_PARAMETER )
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, NEIGHBOUR_WEIGHTS )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, MATERIAL_DENSITY )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, MATERIAL_DENSITY_NEW )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, MATERIAL_DENSITY_FILTERED )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ELEMENT_DC )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ELEMENT_DC_FILTERED )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ELEMENT_DV )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ELEMENT_DV_FILTERED )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, YOUNG_MODULUS_0 )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, YOUNG_MODULUS_MIN )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PENALIZATION_FACTOR )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, GEOMETRICAL_DOMAIN_SIZE )
// KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, JACOBIAN_0 )
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, PRINCIPAL_STRESS )
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, PRINCIPAL_STRAIN )
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, ALGORITHMIC_TANGENT )
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, THREED_ALGORITHMIC_TANGENT )
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, ELASTIC_TANGENT )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, IS_MARKED_FOR_REACTION )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, EXTRAPOLATION_FACTOR_1 )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, EXTRAPOLATION_FACTOR_2 )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, FIX_POROSITY )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, QUAD_POINT_STATUS )
KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, LOCAL_FRAME )
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, STRESS_LIKE_INTERNAL_VARIABLES )
KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, STRAIN_LIKE_INTERNAL_VARIABLES )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PRIMARY_HYDRATION_TIME )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PRIMARY_HYDRATION_TIME_GRADIENT )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, STIFFNESS_RATIO )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, L2_ERROR )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, H1_ERROR )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ROTATIONAL_STIFFNESS )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, RATE_SENSITIVITY )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, LOCAL_ERROR_TOLERANCE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, REFERENCE_STRAIN_RATE )

KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, PLASTIC_MODE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, ACCUMULATED_PLASTIC_STRAIN )

KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, DSTRESS_DTEMPERATURE )
KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, TAYLOR_QUINNEY_FACTOR )

KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, PERTURBATION_FACTOR )

KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, HardeningLaw::Pointer, HARDENING_LAW)
// KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, HardeningLaw::Pointer, HARDENING_LAW_1)
// KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, HardeningLaw::Pointer, HARDENING_LAW_2)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, HardeningLaw::Pointer, ISOTROPIC_HARDENING_LAW)
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, HardeningLaw::Pointer, KINEMATIC_HARDENING_LAW)

KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, INITIAL_DAMAGE )
KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, bool, FIX_DAMAGE )

//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DP_EPSILON )
//  KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, INSITU_STRESS )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DP_ALPHA1 )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, DP_K )
//KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION,TO_ERASE )
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, CALCULATE_INSITU_STRESS )
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, CONTACT_RAMP )
//  KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, PENALTY )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, INITIAL_PENALTY )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, MAXIMUM_PENALTY )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, RAMP_CRITERION )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, RAMP_FACTOR )
//  KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, PENALTY_T )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, INITIAL_PENALTY_T )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, MAXIMUM_PENALTY_T )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, RAMP_CRITERION_T )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, RAMP_FACTOR_T )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, FRICTION_COEFFICIENT )
//  KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, LAMBDAS )
//  KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, LAMBDAS_T )
//  KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, GAPS )
//  KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, DELTA_LAMBDAS )
//  KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, DELTA_LAMBDAS_T )
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, MAX_UZAWA_ITERATIONS)
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, CONTACT_SLAVE_INTEGRATION_POINT_INDEX )
//  KRATOS_DEFINE_APPLICATION_MATRIX_VARIABLE( STRUCTURAL_APPLICATION, CONTACT_LINK_M )
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, CONTACT_DOUBLE_CHECK )
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, IS_CONTACT_MASTER )
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, IS_CONTACT_SLAVE )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, K_CONTACT )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, K_CONTACT_T )
//  KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, STICK )
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, FIRST_TIME_STEP )
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, QUASI_STATIC_ANALYSIS )
//  KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, NORMAL_STRESS )
//  KRATOS_DEFINE_APPLICATION_VECTOR_VARIABLE( STRUCTURAL_APPLICATION, TANGENTIAL_STRESS )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, NORMAL_CONTACT_STRESS )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, TANGENTIAL_CONTACT_STRESS )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, CONTACT_STICK )

//
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, WATER_PRESSURE )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, WATER_PRESSURE_DT )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, WATER_PRESSURE_ACCELERATION )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, WATER_PRESSURE_NULL )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, WATER_PRESSURE_NULL_DT )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, WATER_PRESSURE_NULL_ACCELERATION )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, WATER_PRESSURE_EINS )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, WATER_PRESSURE_EINS_DT )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, WATER_PRESSURE_EINS_ACCELERATION )
//
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AIR_PRESSURE )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AIR_PRESSURE_DT )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AIR_PRESSURE_ACCELERATION )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AIR_PRESSURE_NULL )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AIR_PRESSURE_NULL_DT )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AIR_PRESSURE_NULL_ACCELERATION )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AIR_PRESSURE_EINS )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AIR_PRESSURE_EINS_DT )
//  KRATOS_DEFINE_APPLICATION_DOUBLE_VARIABLE( STRUCTURAL_APPLICATION, AIR_PRESSURE_EINS_ACCELERATION )
//
//  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, DISPLACEMENT_OLD)
//  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, DISPLACEMENT_DT)
//  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, ACCELERATION)
//  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, DISPLACEMENT_NULL)
//  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, DISPLACEMENT_NULL_DT)
//  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, ACCELERATION_NULL)
//  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, DISPLACEMENT_EINS)
//  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, DISPLACEMENT_EINS_DT)
//  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, ACCELERATION_EINS)
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION,  Matrix, ELASTIC_LEFT_CAUCHY_GREEN_OLD )
//
//  KRATOS_DEFINE_APPLICATION_VARIABLE( STRUCTURAL_APPLICATION, int, ACTIVATION_LEVEL)
//KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( STRUCTURAL_APPLICATION, VAUX);

} // namespace Kratos

#endif // KRATOS_STRUCTURAL_APPLICATION_VARIABLES_H_INCLUDED  defined
