// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Apr 1 2016 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "custom_processes/topology_update_process.h"
#include "custom_processes/calculate_reaction_process.h"
#include "custom_processes/calculate_strain_energy_process.h"
#include "custom_processes/arc_length_control_process.h"
#include "custom_utilities/arc_length_constraint.h"
#include "custom_utilities/arc_length_cylinder_constraint.h"
#include "custom_utilities/arc_length_cylinder_scalar_constraint.h"
#include "custom_utilities/arc_length_cylinder_ux_uy_uz_constraint.h"
#include "custom_utilities/arc_length_sphere_constraint.h"
#include "custom_utilities/arc_length_sphere_ux_uy_uz_constraint.h"
#include "custom_utilities/arc_length_load_control_energy_release_constraint.h"
#include "custom_utilities/arc_length_displacement_control_energy_release_constraint.h"
#include "add_custom_processes_to_python.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;

    typedef typename SparseSpaceType::MatrixType SparseMatrixType;
    typedef typename SparseSpaceType::VectorType SparseVectorType;

    class_<TopologyUpdateProcess, bases<Process>, boost::noncopyable>
    ("TopologyUpdateProcess", init<ModelPart&, double, double, double>())
    .def("SetBinSize", &TopologyUpdateProcess::SetBinSize)
    .def("GetTopologyChange", &TopologyUpdateProcess::GetTopologyChange)
    .def("GetObjective", &TopologyUpdateProcess::GetObjective)
    ;

    class_<CalculateReactionProcess, bases<Process>, boost::noncopyable>
    ("CalculateReactionProcess", init<ModelPart&, CalculateReactionProcess::SchemeType&>())
    ;

    class_<CalculateStrainEnergyProcess, bases<Process>, boost::noncopyable>
    ("CalculateStrainEnergyProcess", init<const ModelPart&>())
    .def("GetEnergy", &CalculateStrainEnergyProcess::GetEnergy)
    ;

    typedef ArcLengthConstraint<BuilderAndSolverType> ArcLengthConstraintType;
    class_<ArcLengthConstraintType, ArcLengthConstraintType::Pointer, boost::noncopyable>
    ( "ArcLengthConstraint", init<const double&>() )
    .def("SetRadius", &ArcLengthConstraintType::SetRadius)
    .def("NeedForceVector", &ArcLengthConstraintType::NeedForceVector)
    .def("SetForceVector", &ArcLengthConstraintType::SetForceVector)
    ;

    typedef ArcLengthCylinderConstraint<BuilderAndSolverType> ArcLengthCylinderConstraintType;
    class_<ArcLengthCylinderConstraintType, ArcLengthCylinderConstraintType::Pointer, bases<ArcLengthConstraintType>, boost::noncopyable>
    ( "ArcLengthCylinderConstraint", init<const double&>() )
    ;

    typedef ArcLengthCylinderScalarConstraint<BuilderAndSolverType> ArcLengthCylinderScalarConstraintType;
    class_<ArcLengthCylinderScalarConstraintType, ArcLengthCylinderScalarConstraintType::Pointer, bases<ArcLengthConstraintType>, boost::noncopyable>
    ( "ArcLengthCylinderScalarConstraint", init<const Variable<double>&, const double&>() )
    ;

    typedef ArcLengthCylinderUxUyUzConstraint<BuilderAndSolverType> ArcLengthCylinderUxUyUzConstraintType;
    class_<ArcLengthCylinderUxUyUzConstraintType, ArcLengthCylinderUxUyUzConstraintType::Pointer, bases<ArcLengthConstraintType>, boost::noncopyable>
    ( "ArcLengthCylinderUxUyUzConstraint", init<const double&>() )
    ;

    typedef ArcLengthSphereConstraint<BuilderAndSolverType> ArcLengthSphereConstraintType;
    class_<ArcLengthSphereConstraintType, ArcLengthSphereConstraintType::Pointer, bases<ArcLengthConstraintType>, boost::noncopyable>
    ( "ArcLengthSphereConstraint", init<const double&, const double&>() )
    .def("SetScale", &ArcLengthSphereConstraintType::SetScale)
    ;

    typedef ArcLengthSphereUxUyUzConstraint<BuilderAndSolverType> ArcLengthSphereUxUyUzConstraintType;
    class_<ArcLengthSphereUxUyUzConstraintType, ArcLengthSphereUxUyUzConstraintType::Pointer, bases<ArcLengthConstraintType>, boost::noncopyable>
    ( "ArcLengthSphereUxUyUzConstraint", init<const double&, const double&>() )
    .def("SetScale", &ArcLengthSphereUxUyUzConstraintType::SetScale)
    ;

    typedef ArcLengthLoadControlEnergyReleaseConstraint<BuilderAndSolverType> ArcLengthLoadControlEnergyReleaseConstraintType;
    class_<ArcLengthLoadControlEnergyReleaseConstraintType, ArcLengthLoadControlEnergyReleaseConstraintType::Pointer, bases<ArcLengthConstraintType>, boost::noncopyable>
    ( "ArcLengthLoadControlEnergyReleaseConstraint", init<const double&>() )
    ;

    typedef ArcLengthDisplacementControlEnergyReleaseConstraint<BuilderAndSolverType> ArcLengthDisplacementControlEnergyReleaseConstraintType;
    class_<ArcLengthDisplacementControlEnergyReleaseConstraintType, ArcLengthDisplacementControlEnergyReleaseConstraintType::Pointer, bases<ArcLengthConstraintType>, boost::noncopyable>
    ( "ArcLengthDisplacementControlEnergyReleaseConstraint", init<const double&>() )
    ;

    typedef ArcLengthControlProcess<BuilderAndSolverType> ArcLengthControlProcessType;
    void(ArcLengthControlProcessType::*ArcLengthControlProcess_Execute)(SparseVectorType&, const SparseVectorType&) = &ArcLengthControlProcessType::Execute;
    class_<ArcLengthControlProcessType, ArcLengthControlProcessType::Pointer, bases<Process>, boost::noncopyable>
    ( "ArcLengthControlProcess", init<ArcLengthConstraintType::Pointer>() )
    .def("SetPredictor", &ArcLengthControlProcessType::SetPredictor)
    .def("SetForcedMode", &ArcLengthControlProcessType::SetForcedMode)
    .def("SetSolveMode", &ArcLengthControlProcessType::SetSolveMode)
    .def("GetSolveMode", &ArcLengthControlProcessType::GetSolveMode)
    .def("SetModelPart", &ArcLengthControlProcessType::SetModelPart)
    .def("SetBuilderAndSolver", &ArcLengthControlProcessType::SetBuilderAndSolver)
    .def("Update", &ArcLengthControlProcessType::Update)
    .def("GetLambda", &ArcLengthControlProcessType::GetLambda)
    .def("GetLambdaOld", &ArcLengthControlProcessType::GetLambdaOld)
    .def("GetDeltaLambda", &ArcLengthControlProcessType::GetDeltaLambda)
    .def("GetDeltaLambdaOld", &ArcLengthControlProcessType::GetDeltaLambdaOld)
    .def("GetConstraint", &ArcLengthControlProcessType::pGetConstraint)
    .def("Reset", &ArcLengthControlProcessType::Reset)
    // .def("GetValue", &ArcLengthControlProcessType::GetValue)
    // .def("GetDerivativesDU", &ArcLengthControlProcessType::GetDerivativesDU)
    // .def("GetDerivativesDLambda", &ArcLengthControlProcessType::GetDerivativesDLambda)
    // .def("Predict", &ArcLengthControlProcessType::Predict)
    .def("Execute", ArcLengthControlProcess_Execute)
    .def("Copy", &ArcLengthControlProcessType::Copy)
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

} // Namespace Kratos
