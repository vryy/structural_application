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
//   Date:                $Date: 2009-03-20 08:55:31 $
//   Revision:            $Revision: 1.23 $
//
//



// System includes


// External includes


// Project includes
#include "structural_application_variables.h"

namespace Kratos
{

//Example
typedef Matrix fix_matrix_33;
//typedef boost::numeric::ublas::bounded_matrix<double,3,3> fix_matrix_33;
typedef Vector array3;
//typedef array_1d<double,3> array3;
KRATOS_CREATE_VARIABLE( fix_matrix_33 , MATRIX_A )
KRATOS_CREATE_VARIABLE( fix_matrix_33 , MATRIX_B )
KRATOS_CREATE_VARIABLE( fix_matrix_33 , MATRIX_D )
KRATOS_CREATE_VARIABLE( array3, COMPOSITE_DIRECTION )
KRATOS_CREATE_VARIABLE( array3, ORTHOTROPIC_YOUNG_MODULUS )
KRATOS_CREATE_VARIABLE( array3, ORTHOTROPIC_SHEAR_MODULUS )
KRATOS_CREATE_VARIABLE( Matrix, ORTHOTROPIC_POISSON_RATIO )
KRATOS_CREATE_VARIABLE( Matrix , GEOMETRIC_STIFFNESS )
KRATOS_CREATE_VARIABLE( Matrix , MATERIAL_DIRECTION )
KRATOS_CREATE_VARIABLE( Matrix , JOINT_STIFFNESS )
KRATOS_CREATE_VARIABLE( double , LINING_JOINT_STIFFNESS )

KRATOS_CREATE_VARIABLE( double, DAMAGE_E0 )
KRATOS_CREATE_VARIABLE( double, DAMAGE_EF )

KRATOS_CREATE_VARIABLE( ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_NO_INITIALIZE )

#ifdef SD_APP_FORWARD_COMPATIBILITY
KRATOS_CREATE_VARIABLE( double, ALPHA )
KRATOS_CREATE_VARIABLE( double, KAPPA )
#endif

//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VAUX);

// KRATOS_CREATE_VARIABLE(int, WRINKLING_APPROACH )
// KRATOS_CREATE_VARIABLE(Matrix, GREEN_LAGRANGE_STRAIN_TENSOR )
// KRATOS_CREATE_VARIABLE(Matrix, PK2_STRESS_TENSOR )
// KRATOS_CREATE_VARIABLE(Matrix, AUXILIARY_MATRIX_1 )
// KRATOS_CREATE_VARIABLE(double, YOUNG_MODULUS )
// KRATOS_CREATE_VARIABLE(double, POISSON_RATIO )
// KRATOS_CREATE_VARIABLE(double, MU )
KRATOS_CREATE_VARIABLE( double, RETRACTION_TIME )
// KRATOS_CREATE_VARIABLE(double, THICKNESS )
// KRATOS_CREATE_VARIABLE(double, NEGATIVE_FACE_PRESSURE )
// KRATOS_CREATE_VARIABLE(double, POSITIVE_FACE_PRESSURE )


//  KRATOS_CREATE_VARIABLE(double, DP_EPSILON)
//  KRATOS_CREATE_VARIABLE(Vector, INSITU_STRESS )
//  KRATOS_CREATE_VARIABLE(double, DP_ALPHA1 )
//  KRATOS_CREATE_VARIABLE(double, DP_K )
//KRATOS_CREATE_VARIABLE(double,TO_ERASE )
//  KRATOS_CREATE_VARIABLE(int, CALCULATE_INSITU_STRESS )
//CONTACT_LINK_MASTER is defined in condition.h
KRATOS_CREATE_VARIABLE( Condition::Pointer, CONTACT_LINK_MASTER )
//CONTACT_LINK_SLAVE is defined in condition.h
KRATOS_CREATE_VARIABLE( Condition::Pointer, CONTACT_LINK_SLAVE )
KRATOS_CREATE_VARIABLE( Node<3>::Pointer,   NEAR_NODE )
KRATOS_CREATE_VARIABLE( Point<3>, MASTER_CONTACT_LOCAL_POINT )
KRATOS_CREATE_VARIABLE( Point<3>, MASTER_CONTACT_CURRENT_LOCAL_POINT )
KRATOS_CREATE_VARIABLE( Point<3>, MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT )
KRATOS_CREATE_VARIABLE( Point<3>, SLAVE_CONTACT_LOCAL_POINT )
KRATOS_CREATE_VARIABLE( Point<3>, MASTER_CONTACT_GLOBAL_POINT )
KRATOS_CREATE_VARIABLE( Point<3>, MASTER_CONTACT_CURRENT_GLOBAL_POINT )
KRATOS_CREATE_VARIABLE( Point<3>, SLAVE_CONTACT_GLOBAL_POINT )
KRATOS_CREATE_VARIABLE( double, INSITU_STRESS_SCALE )
KRATOS_CREATE_VARIABLE( double, LAGRANGE_SCALE )
KRATOS_CREATE_VARIABLE( double, REFERENCE_PRESSURE )
KRATOS_CREATE_VARIABLE( double, REFERENCE_WATER_PRESSURE )
KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_SCALE )
KRATOS_CREATE_VARIABLE( double, OVERCONSOLIDATION_RATIO )
KRATOS_CREATE_VARIABLE( double, EXCESS_PORE_WATER_PRESSURE )
KRATOS_CREATE_VARIABLE( double, PRESSURE_P )
KRATOS_CREATE_VARIABLE( double, PRESSURE_Q )
KRATOS_CREATE_VARIABLE( Vector, COORDINATES )
KRATOS_CREATE_VARIABLE( Vector, FLUID_FLOWS )
KRATOS_CREATE_VARIABLE( double, CONTACT_PENETRATION )

KRATOS_CREATE_VARIABLE( double, BASE )
KRATOS_CREATE_VARIABLE( double, HEIGHT )
KRATOS_CREATE_VARIABLE( double, CROSS_AREA )
KRATOS_CREATE_VARIABLE( double, AREA )
KRATOS_CREATE_VARIABLE( double, AREA_X )
KRATOS_CREATE_VARIABLE( double, AREA_Y )
KRATOS_CREATE_VARIABLE( double, AREA_Z )
KRATOS_CREATE_VARIABLE( Matrix, LOCAL_INERTIA )
KRATOS_CREATE_VARIABLE( double, INERTIA_X )
KRATOS_CREATE_VARIABLE( double, INERTIA_Y )
KRATOS_CREATE_VARIABLE( double, INERTIA_Z )
KRATOS_CREATE_VARIABLE( double, FC )
KRATOS_CREATE_VARIABLE( double, FT )
KRATOS_CREATE_VARIABLE( double, CONCRETE_YOUNG_MODULUS_C )
KRATOS_CREATE_VARIABLE( double, CONCRETE_YOUNG_MODULUS_T )
KRATOS_CREATE_VARIABLE( double, FRACTURE_ENERGY )
KRATOS_CREATE_VARIABLE( double, CRUSHING_ENERGY )
KRATOS_CREATE_VARIABLE( double, ELASTIC_ENERGY )
KRATOS_CREATE_VARIABLE( double, PLASTIC_ENERGY )
//     KRATOS_CREATE_VARIABLE( double, YIELD_STRESS )
KRATOS_CREATE_VARIABLE( double, PLASTIC_MODULUS )
KRATOS_CREATE_VARIABLE( double, PLASTICITY_INDICATOR )
KRATOS_CREATE_VARIABLE( double, ISOTROPIC_HARDENING_MODULUS )
KRATOS_CREATE_VARIABLE( double, KINEMATIC_HARDENING_MODULUS )
KRATOS_CREATE_VARIABLE( int, ISOTROPIC_HARDENING_TYPE )
KRATOS_CREATE_VARIABLE( int, KINEMATIC_HARDENING_TYPE )
KRATOS_CREATE_VARIABLE( int, DRUCKER_PRAGER_MATCHING_TYPE )
KRATOS_CREATE_VARIABLE( Matrix, HARDENING_POINTS_ON_CURVE )
KRATOS_CREATE_VARIABLE( double, LAMNDA ) // Load factor
KRATOS_CREATE_VARIABLE( double, DAMAGE )
KRATOS_CREATE_VARIABLE( double, ORTHOTROPIC_ANGLE )
KRATOS_CREATE_VARIABLE( double, VOLUMEN_FRACTION )
KRATOS_CREATE_VARIABLE( double, MAX_INTERNAL_FRICTION_ANGLE )
KRATOS_CREATE_VARIABLE( double, DILATANCY_ANGLE )
KRATOS_CREATE_VARIABLE( double, MAX_DILATANCY_ANGLE )
KRATOS_CREATE_VARIABLE( double, COHESION )
KRATOS_CREATE_VARIABLE( double, DISIPATION )
KRATOS_CREATE_VARIABLE( double, ISOTROPIC_ELASTIC_LIMIT )
KRATOS_CREATE_VARIABLE( Vector, ORTHOTROPIC_ELASTIC_LIMIT )
KRATOS_CREATE_VARIABLE( Vector, VECTOR_DAMAGE )
KRATOS_CREATE_VARIABLE( Vector, ORTHOTROPIC_YOUNG_MODULUS_2D ) // [E1 E2 G12]
KRATOS_CREATE_VARIABLE( Vector, ORTHOTROPIC_POISSON_RATIO_2D ) // [v12 v21]
KRATOS_CREATE_VARIABLE( Matrix, GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
KRATOS_CREATE_VARIABLE( Vector, ELASTIC_STRAIN_VECTOR )
KRATOS_CREATE_VARIABLE( Vector, PLASTIC_STRAIN_VECTOR )
KRATOS_CREATE_VARIABLE( Vector, CURRENT_STRAIN_VECTOR )
KRATOS_CREATE_VARIABLE( Matrix, CURRENT_DEFORMATION_GRADIENT )
KRATOS_CREATE_VARIABLE( double, CURRENT_DEFORMATION_GRADIENT_DETERMINANT )
KRATOS_CREATE_VARIABLE( Vector, INTEGRATION_POINT_STRAIN_VECTOR )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_VOLUMETRIC_ELASTIC_STRAIN )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_VOLUMETRIC_PLASTIC_STRAIN )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_DEVIATORIC_ELASTIC_STRAIN )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_DEVIATORIC_PLASTIC_STRAIN )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_VOLUMETRIC_STRAIN )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_DEVIATORIC_STRAIN )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_STRAIN )
KRATOS_CREATE_VARIABLE( Vector, ALMANSI_PLASTIC_STRAIN )
KRATOS_CREATE_VARIABLE( Vector, ALMANSI_ELASTIC_STRAIN )
KRATOS_CREATE_VARIABLE( Matrix, NODAL_STRESS )
KRATOS_CREATE_VARIABLE( Vector, NODAL_STRESS_VECTOR )
KRATOS_CREATE_VARIABLE( Matrix, NODAL_STRAIN )
KRATOS_CREATE_VARIABLE( Vector, NODAL_STRAIN_VECTOR )
KRATOS_CREATE_VARIABLE( Matrix, CONSTRAINT_MATRIX )
KRATOS_CREATE_VARIABLE( Vector, PRESTRAIN )
KRATOS_CREATE_VARIABLE( Vector, PRESTRESS )
KRATOS_CREATE_VARIABLE( double, PRESTRESS_ZZ )
KRATOS_CREATE_VARIABLE( double, PRESTRESS_FACTOR )
KRATOS_CREATE_VARIABLE( Vector, INITIAL_STRESS )
KRATOS_CREATE_VARIABLE( Vector, CONSTRAINT_VECTOR )
KRATOS_CREATE_VARIABLE( int,    NODAL_VALUES )
KRATOS_CREATE_VARIABLE( double, NODAL_DAMAGE )
KRATOS_CREATE_VARIABLE( bool,   IS_TARGET )
KRATOS_CREATE_VARIABLE( bool,   IS_CONTACTOR )
KRATOS_CREATE_VARIABLE( bool,   COMPUTE_TANGENT_MATRIX )
KRATOS_CREATE_VARIABLE( bool,   USE_NUMERICAL_TANGENT )
KRATOS_CREATE_VARIABLE( double,   IS_DISCRETE )
KRATOS_CREATE_VARIABLE( double, DAMPING_RATIO )
//KRATOS_CREATE_VARIABLE( double, KINETIC_ENERGY )
KRATOS_CREATE_VARIABLE( double, POTENCIAL_ENERGY )
KRATOS_CREATE_VARIABLE( double, DEFORMATION_ENERGY )
KRATOS_CREATE_VARIABLE( double, VON_MISES_STRESS )
KRATOS_CREATE_VARIABLE( double, RHS_PRESSURE )
KRATOS_CREATE_VARIABLE( double,  TENSILE_STRENGTH )
KRATOS_CREATE_VARIABLE( double,  SHEAR_STRENGTH )
KRATOS_CREATE_VARIABLE( double,  VISCOUS_DAMPING )
KRATOS_CREATE_VARIABLE( int,  YIELD_STATE )
KRATOS_CREATE_VARIABLE( double,  YIELD_SURFACE )
KRATOS_CREATE_VARIABLE( double,  MAX_FRECUENCY )

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( JOINT_FORCE_REACTION )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( JOINT_MOMENT_REACTION )
//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( INTERNAL_FORCE ) //already put on variables.cpp (warning was appearing on Windows)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( ELASTIC_BEDDING_STIFFNESS )

KRATOS_CREATE_VARIABLE(bool, IS_BBAR )
KRATOS_CREATE_VARIABLE(bool, IS_FBAR )
KRATOS_CREATE_VARIABLE(int, FBAR_MODE )
KRATOS_CREATE_VARIABLE(int, NEIGHBOUR_EXPANSION_LEVEL )
KRATOS_CREATE_VARIABLE(int, STRESS_RECOVERY_TYPE )
KRATOS_CREATE_VARIABLE(double, RAYLEIGH_DAMPING_ALPHA )
KRATOS_CREATE_VARIABLE(double, RAYLEIGH_DAMPING_BETA )
KRATOS_CREATE_VARIABLE(double, STABILISATION_FACTOR )
KRATOS_CREATE_VARIABLE(double, SHEAR_MODULUS )
KRATOS_CREATE_VARIABLE(double, SHEAR_MODULUS_EVOLUTION )
KRATOS_CREATE_VARIABLE(Vector, RECOVERY_STRESSES )
KRATOS_CREATE_VARIABLE(Vector, THREED_STRESSES )
KRATOS_CREATE_VARIABLE(Vector, THREED_PRESTRESS )
KRATOS_CREATE_VARIABLE(Vector, STRESSES_OLD )
KRATOS_CREATE_VARIABLE(Vector, STRAIN_OLD )
KRATOS_CREATE_VARIABLE(Vector, THREED_STRAIN )
KRATOS_CREATE_VARIABLE(Vector, PRE_STRAIN_VECTOR )
KRATOS_CREATE_VARIABLE(Vector, POST_STRAIN_VECTOR )
KRATOS_CREATE_VARIABLE(Matrix, STRESS_TENSOR )
KRATOS_CREATE_VARIABLE(Matrix, STRAIN_TENSOR )
KRATOS_CREATE_VARIABLE(Matrix, ELASTIC_STRAIN_TENSOR )
KRATOS_CREATE_VARIABLE(Matrix, ELASTIC_STRAIN_TENSOR_OLD )
KRATOS_CREATE_VARIABLE(Matrix, PLASTIC_STRAIN_TENSOR )
KRATOS_CREATE_VARIABLE(Matrix, LEFT_STRETCH_TENSOR )
KRATOS_CREATE_VARIABLE(Matrix, RIGHT_STRETCH_TENSOR )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESCRIBED_DELTA_DISPLACEMENT)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESCRIBED_DISPLACEMENT)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(INITIAL_DISPLACEMENT)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
KRATOS_CREATE_VARIABLE(bool,   IS_CONTACT_NODE)
KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER)
KRATOS_CREATE_VARIABLE(int, LAGRANGE_MULTIPLIER_INDEX)
KRATOS_CREATE_VARIABLE(bool, HAS_STRAIN_AT_NODE)
KRATOS_CREATE_VARIABLE(bool, HAS_STRESSES_AT_NODE)
KRATOS_CREATE_VARIABLE( bool, HAS_NODAL_ERROR )
KRATOS_CREATE_VARIABLE( bool, FORCE_EQUAL_ORDER_INTERPOLATION )
KRATOS_CREATE_VARIABLE( double, NODAL_ERROR_1 )
KRATOS_CREATE_VARIABLE( double, DUMMY_DOF )

KRATOS_CREATE_VARIABLE( double, PRECONSOLIDATION_PRESSURE )
KRATOS_CREATE_VARIABLE( double, PRECONSOLIDATION_PRESSURE_MIN )
KRATOS_CREATE_VARIABLE( double, PRECONSOLIDATION_PRESSURE_DEF )
KRATOS_CREATE_VARIABLE( double, CSL_SLOPE )
KRATOS_CREATE_VARIABLE( double, VIRGIN_COMPRESSION_INDEX )
KRATOS_CREATE_VARIABLE( double, RECOMPRESSION_INDEX )
KRATOS_CREATE_VARIABLE( double, SWELL_INDEX )
KRATOS_CREATE_VARIABLE( double, VOID_RATIO )
KRATOS_CREATE_VARIABLE( double, SPACING_RATIO )
KRATOS_CREATE_VARIABLE( double, ASSOCIATIVITY )
KRATOS_CREATE_VARIABLE( double, SHAPE_PARAMETER )

KRATOS_CREATE_VARIABLE( Vector, NEIGHBOUR_WEIGHTS )
KRATOS_CREATE_VARIABLE( double, MATERIAL_DENSITY )
KRATOS_CREATE_VARIABLE( double, MATERIAL_DENSITY_NEW )
KRATOS_CREATE_VARIABLE( double, MATERIAL_DENSITY_FILTERED )
KRATOS_CREATE_VARIABLE( double, ELEMENT_DC )
KRATOS_CREATE_VARIABLE( double, ELEMENT_DC_FILTERED )
KRATOS_CREATE_VARIABLE( double, ELEMENT_DV )
KRATOS_CREATE_VARIABLE( double, ELEMENT_DV_FILTERED )
KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_0 )
KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_MIN )
KRATOS_CREATE_VARIABLE( double, PENALIZATION_FACTOR )
KRATOS_CREATE_VARIABLE( double, GEOMETRICAL_DOMAIN_SIZE )
// KRATOS_CREATE_VARIABLE( double, JACOBIAN_0 )
KRATOS_CREATE_VARIABLE( Vector, PRINCIPAL_STRESS )
KRATOS_CREATE_VARIABLE( Vector, PRINCIPAL_STRAIN )
KRATOS_CREATE_VARIABLE( Matrix, ALGORITHMIC_TANGENT )
KRATOS_CREATE_VARIABLE( Matrix, THREED_ALGORITHMIC_TANGENT )
KRATOS_CREATE_VARIABLE( Matrix, ELASTIC_TANGENT )
KRATOS_CREATE_VARIABLE( bool, IS_MARKED_FOR_REACTION )
KRATOS_CREATE_VARIABLE( double, EXTRAPOLATION_FACTOR_1 )
KRATOS_CREATE_VARIABLE( double, EXTRAPOLATION_FACTOR_2 )
KRATOS_CREATE_VARIABLE( bool, FIX_POROSITY )
KRATOS_CREATE_VARIABLE( int, QUAD_POINT_STATUS )
KRATOS_CREATE_VARIABLE( Matrix, LOCAL_FRAME )
KRATOS_CREATE_VARIABLE( Vector, STRESS_LIKE_INTERNAL_VARIABLES )
KRATOS_CREATE_VARIABLE( Vector, STRAIN_LIKE_INTERNAL_VARIABLES )
KRATOS_CREATE_VARIABLE ( double, PRIMARY_HYDRATION_TIME )
KRATOS_CREATE_VARIABLE ( double, PRIMARY_HYDRATION_TIME_GRADIENT )
KRATOS_CREATE_VARIABLE ( double, STIFFNESS_RATIO )
KRATOS_CREATE_VARIABLE ( double, L2_ERROR )
KRATOS_CREATE_VARIABLE ( double, H1_ERROR )
KRATOS_CREATE_VARIABLE ( double, ROTATIONAL_STIFFNESS )
KRATOS_CREATE_VARIABLE ( double, RATE_SENSITIVITY )
KRATOS_CREATE_VARIABLE ( double, LOCAL_ERROR_TOLERANCE )
KRATOS_CREATE_VARIABLE ( double, REFERENCE_STRAIN_RATE )

KRATOS_CREATE_VARIABLE ( int, PLASTIC_MODE )
KRATOS_CREATE_VARIABLE ( double, ACCUMULATED_PLASTIC_STRAIN )

KRATOS_CREATE_VARIABLE ( Vector, DSTRESS_DTEMPERATURE )
KRATOS_CREATE_VARIABLE ( double, TAYLOR_QUINNEY_FACTOR )

KRATOS_CREATE_VARIABLE ( double, PERTURBATION_FACTOR )

KRATOS_CREATE_VARIABLE( HardeningLaw::Pointer, HARDENING_LAW )
KRATOS_CREATE_VARIABLE( HardeningLaw::Pointer, ISOTROPIC_HARDENING_LAW )
KRATOS_CREATE_VARIABLE( HardeningLaw::Pointer, KINEMATIC_HARDENING_LAW )

KRATOS_CREATE_VARIABLE ( double, INITIAL_DAMAGE )
KRATOS_CREATE_VARIABLE ( bool, FIX_DAMAGE )

//  KRATOS_CREATE_VARIABLE(int, CONTACT_RAMP )
//  KRATOS_CREATE_VARIABLE(Vector, PENALTY )
// //  KRATOS_CREATE_VARIABLE(double, INITIAL_PENALTY )
//  KRATOS_CREATE_VARIABLE(double, MAXIMUM_PENALTY )
//  KRATOS_CREATE_VARIABLE(double, RAMP_CRITERION )
//  KRATOS_CREATE_VARIABLE(double, RAMP_FACTOR )
//  KRATOS_CREATE_VARIABLE(Vector, PENALTY_T )
//  KRATOS_CREATE_VARIABLE(double, INITIAL_PENALTY_T )
//  KRATOS_CREATE_VARIABLE(double, MAXIMUM_PENALTY_T )
//  KRATOS_CREATE_VARIABLE(double, RAMP_CRITERION_T )
//  KRATOS_CREATE_VARIABLE(double, RAMP_FACTOR_T )
//  KRATOS_CREATE_VARIABLE(double, FRICTION_COEFFICIENT )
//  KRATOS_CREATE_VARIABLE(Vector, LAMBDAS )
//  KRATOS_CREATE_VARIABLE(Matrix, LAMBDAS_T )
//  KRATOS_CREATE_VARIABLE(Vector, GAPS )
//  KRATOS_CREATE_VARIABLE(Vector, DELTA_LAMBDAS )
//  KRATOS_CREATE_VARIABLE(Matrix, DELTA_LAMBDAS_T )
//  KRATOS_CREATE_VARIABLE(int, MAX_UZAWA_ITERATIONS)
//  KRATOS_CREATE_VARIABLE(int, CONTACT_SLAVE_INTEGRATION_POINT_INDEX )
//  KRATOS_CREATE_VARIABLE( Matrix, CONTACT_LINK_M )
//  KRATOS_CREATE_VARIABLE( int, CONTACT_DOUBLE_CHECK )
//  KRATOS_CREATE_VARIABLE( int, IS_CONTACT_MASTER )
//  KRATOS_CREATE_VARIABLE( int, IS_CONTACT_SLAVE )
//  KRATOS_CREATE_VARIABLE( double, K_CONTACT )
//  KRATOS_CREATE_VARIABLE( double, K_CONTACT_T )
//  KRATOS_CREATE_VARIABLE( Vector, STICK )
//  KRATOS_CREATE_VARIABLE( int, FIRST_TIME_STEP )
//  KRATOS_CREATE_VARIABLE( int, QUASI_STATIC_ANALYSIS )
//  KRATOS_CREATE_VARIABLE( Vector, NORMAL_STRESS )
//  KRATOS_CREATE_VARIABLE( Vector, TANGENTIAL_STRESS )
//  KRATOS_CREATE_VARIABLE( double, NORMAL_CONTACT_STRESS )
//  KRATOS_CREATE_VARIABLE( double, TANGENTIAL_CONTACT_STRESS )
//  KRATOS_CREATE_VARIABLE( double, CONTACT_STICK )
//
//  KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE )
//  KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_DT )
//  KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_ACCELERATION )
//  KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_NULL )
//  KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_NULL_DT )
//  KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_NULL_ACCELERATION )
//  KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_EINS )
//  KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_EINS_DT )
//  KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_EINS_ACCELERATION )
//
//  KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE )
//  KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_DT )
//  KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_ACCELERATION )
//  KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_NULL )
//  KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_NULL_DT )
//  KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_NULL_ACCELERATION )
//  KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_EINS )
//  KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_EINS_DT )
//  KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_EINS_ACCELERATION )
//
//  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_OLD)
//  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_DT)
// //  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION)
//  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_NULL)
//  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_NULL_DT)
//  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION_NULL)
//  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_EINS)
//  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_EINS_DT)
//  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION_EINS)
//  KRATOS_CREATE_VARIABLE( Matrix, ELASTIC_LEFT_CAUCHY_GREEN_OLD )
//
//  KRATOS_CREATE_VARIABLE(int, ACTIVATION_LEVEL)

}  // namespace Kratos.
