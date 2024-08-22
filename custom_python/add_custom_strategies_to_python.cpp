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
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-15 18:49:07 $
//   Revision:            $Revision: 1.20 $
//
//


// System includes


// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "custom_python/add_custom_strategies_to_python.h"

//convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/convergencecriterias/multiphaseflow_criteria.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/residualbased_incrementalupdate_static_deactivation_scheme.h"
#include "custom_strategies/schemes/residualbased_newmark_scheme.h"
#include "custom_strategies/schemes/residualbased_theta_scheme.h"
#include "custom_strategies/schemes/residualbased_state_based_theta_scheme.h"
#include "custom_strategies/schemes/residualbased_central_difference_scheme.h"
#include "custom_strategies/schemes/residualbased_acc_based_forward_euler_scheme.h"
#include "custom_strategies/schemes/residualbased_acc_based_central_difference_scheme.h"
#include "custom_strategies/schemes/residualbased_mixed_forward_euler_scheme.h"
#include "custom_strategies/schemes/arc_length_displacement_control_support_scheme.h"
#include "custom_strategies/schemes/arc_length_displacement_control_energy_release_support_scheme.h"

//builder_and_solvers
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/modal_analysis_builder_and_solver.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

template<class TSchemeType>
void AddArcLengthDisplacementControlSupportScheme(const std::string& PostFix)
{
    typedef ArcLengthDisplacementControlSupportScheme<TSchemeType> ArcLengthDisplacementControlSupportSchemeType;

    std::stringstream Name1;
    Name1 << "ArcLengthDisplacementControl" << PostFix;
    class_< ArcLengthDisplacementControlSupportSchemeType, bases< TSchemeType >,  boost::noncopyable >
    ( Name1.str().c_str(), init<typename TSchemeType::Pointer>() )
    .def("GetForceVector", &ArcLengthDisplacementControlSupportSchemeType::GetForceVector, return_internal_reference<>())
    ;

    typedef ArcLengthDisplacementControlEnergyReleaseSupportScheme<TSchemeType> ArcLengthDisplacementControlEnergyReleaseSupportSchemeType;

    std::stringstream Name2;
    Name2 << "ArcLengthDisplacementControlEnergyRelease" << PostFix;
    class_< ArcLengthDisplacementControlEnergyReleaseSupportSchemeType, bases< ArcLengthDisplacementControlSupportSchemeType >,  boost::noncopyable >
    ( Name2.str().c_str(), init<typename TSchemeType::Pointer>() )
    .def("GetForceVector2", &ArcLengthDisplacementControlEnergyReleaseSupportSchemeType::GetForceVector2, return_internal_reference<>())
    ;
}

void  AddCustomStrategiesToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef ResidualBasedIncrementalUpdateStaticDeactivationScheme< SparseSpaceType, LocalSpaceType > ResidualBasedIncrementalUpdateStaticDeactivationSchemeType;
    typedef ResidualBasedNewmarkScheme< SparseSpaceType, LocalSpaceType > ResidualBasedNewmarkSchemeType;
    typedef ResidualBasedNewmarkScheme< SparseSpaceType, LocalSpaceType, 1 > ResidualBasedNonlinearMassDampingNewmarkSchemeType;
    typedef ResidualBasedThetaScheme< SparseSpaceType, LocalSpaceType > ResidualBasedThetaSchemeType;
    typedef ResidualBasedStateBasedThetaScheme< SparseSpaceType, LocalSpaceType > ResidualBasedStateBasedThetaSchemeType;
    typedef ResidualBasedCentralDifferenceScheme< SparseSpaceType, LocalSpaceType > ResidualBasedCentralDifferenceSchemeType;
    typedef ResidualBasedAccBasedForwardEulerScheme< SparseSpaceType, LocalSpaceType > ResidualBasedAccBasedForwardEulerSchemeType;
    typedef ResidualBasedAccBasedCentralDifferenceScheme< SparseSpaceType, LocalSpaceType > ResidualBasedAccBasedCentralDifferenceSchemeType;
    typedef ResidualBasedMixedForwardEulerScheme< SparseSpaceType, LocalSpaceType > ResidualBasedMixedForwardEulerSchemeType;

    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;

    typedef MultiPhaseFlowCriteria< SparseSpaceType,  LocalSpaceType > MultiPhaseFlowCriteriaType;

    typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;

    typedef ModalAnalysisBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ModalAnalysisBuilderAndSolverType;

    //********************************************************************
    //********************************************************************

    class_< ResidualBasedIncrementalUpdateStaticDeactivationSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
            "ResidualBasedIncrementalUpdateStaticDeactivationScheme", init< >()
            );

    class_< ResidualBasedNewmarkSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedNewmarkScheme", init< double >()
            )
            .def(init<>())
            .def(init<int, double>())
            .def("SetIntegrateRotation", &ResidualBasedNewmarkSchemeType::SetIntegrateRotation)
            .def("SetIntegrateMultiplier", &ResidualBasedNewmarkSchemeType::SetIntegrateMultiplier)
            .def("SetIntegrateLoad", &ResidualBasedNewmarkSchemeType::SetIntegrateLoad)
            .def("UpdateForces", &ResidualBasedNewmarkSchemeType::UpdateForces)
            ;

    class_< ResidualBasedNonlinearMassDampingNewmarkSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedNonlinearMassDampingNewmarkScheme", init< double >()
            )
            .def(init<>())
            .def(init<int, double>())
            .def("SetIntegrateRotation", &ResidualBasedNonlinearMassDampingNewmarkSchemeType::SetIntegrateRotation)
            .def("SetIntegrateMultiplier", &ResidualBasedNonlinearMassDampingNewmarkSchemeType::SetIntegrateMultiplier)
            .def("SetIntegrateLoad", &ResidualBasedNonlinearMassDampingNewmarkSchemeType::SetIntegrateLoad)
            .def("UpdateForces", &ResidualBasedNonlinearMassDampingNewmarkSchemeType::UpdateForces)
            ;

    class_< ResidualBasedThetaSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedThetaScheme", init< double >()
            )
            .def(init<>())
            .def("SetIntegrateRotation", &ResidualBasedThetaSchemeType::SetIntegrateRotation)
            .def("SetIntegrateMultiplier", &ResidualBasedThetaSchemeType::SetIntegrateMultiplier)
            .def("SetIntegrateLoad", &ResidualBasedThetaSchemeType::SetIntegrateLoad)
            .def("UpdateForces", &ResidualBasedThetaSchemeType::UpdateForces)
            ;

    class_< ResidualBasedStateBasedThetaSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedStateBasedThetaScheme", init< double >()
            )
            .def(init<>())
            .def("SetIntegrateRotation", &ResidualBasedStateBasedThetaSchemeType::SetIntegrateRotation)
            .def("SetIntegrateMultiplier", &ResidualBasedStateBasedThetaSchemeType::SetIntegrateMultiplier)
            .def("SetIntegrateLoad", &ResidualBasedStateBasedThetaSchemeType::SetIntegrateLoad)
            .def("UpdateForces", &ResidualBasedStateBasedThetaSchemeType::UpdateForces)
            ;

    class_< ResidualBasedCentralDifferenceSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedCentralDifferenceScheme", init<>()
            )
            ;

    class_< ResidualBasedAccBasedForwardEulerSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedAccBasedForwardEulerScheme", init<>()
            )
            .def(init<const bool>())
            ;

    class_< ResidualBasedAccBasedCentralDifferenceSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedAccBasedCentralDifferenceScheme", init<>()
            )
            .def(init<const bool>())
            ;

    class_< ResidualBasedMixedForwardEulerSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedMixedForwardEulerScheme", init<>()
            )
            .def(init<const bool>())
            ;

    AddArcLengthDisplacementControlSupportScheme<ResidualBasedIncrementalUpdateStaticDeactivationSchemeType>("ResidualBasedIncrementalUpdateStaticDeactivationScheme");

    class_< MultiPhaseFlowCriteriaType,
            bases< ConvergenceCriteriaBaseType >, boost::noncopyable >
            ("MultiPhaseFlowCriteria", init<double, double >() )
            .def("SetType", &MultiPhaseFlowCriteriaType::SetType)
            ;

    class_< ModalAnalysisBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable>
    ("ModalAnalysisBuilderAndSolver", init<LinearSolverType::Pointer>())
    .def("SetMaxEigenSolutions", &ModalAnalysisBuilderAndSolverType::SetMaxEigenSolutions)
    .def("SetMaxIterations", &ModalAnalysisBuilderAndSolverType::SetMaxIterations)
    .def("SetTolerance", &ModalAnalysisBuilderAndSolverType::SetTolerance)
    .def("ResizeAndInitializeEigenSystem", &ModalAnalysisBuilderAndSolverType::ResizeAndInitializeEigenSystem)
    .def("BuildEigenSystem", &ModalAnalysisBuilderAndSolverType::BuildEigenSystem)
    ;
}
}  // namespace Python.

} // Namespace Kratos

