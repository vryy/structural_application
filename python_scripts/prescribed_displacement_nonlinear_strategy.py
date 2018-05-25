from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.mpi import *
CheckForPreviousImport()

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(ELASTIC_BEDDING_STIFFNESS)
    model_part.AddNodalSolutionStepVariable(INTERNAL_FORCE)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT)
    model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL)
    model_part.AddNodalSolutionStepVariable(ACCELERATION_EINS)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(ELASTIC_LEFT_CAUCHY_GREEN_OLD)
    model_part.AddNodalSolutionStepVariable(INSITU_STRESS)
    model_part.AddNodalSolutionStepVariable(PRESTRESS)
    model_part.AddNodalSolutionStepVariable(STRESSES)
    model_part.AddNodalSolutionStepVariable(STRAIN)
    model_part.AddNodalSolutionStepVariable(FACE_LOAD)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_NULL)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_EINS)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_DT)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_NULL_DT)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_EINS_DT)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_NULL_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_EINS_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION_AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_DT)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL_DT)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS_DT)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(EXCESS_PORE_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(PRESCRIBED_DELTA_DISPLACEMENT)
    #auxiliary variables misused for mesh rezoning ;-)
    model_part.AddNodalSolutionStepVariable(IS_VISITED)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_MULTIPLIER)
#    model_part.AddNodalSolutionStepVariable(GAP)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_WATER_PRESSURE)
    #model_part.AddNodalSolutionStepVariable(INTERNAL_VARIABLES)
    model_part.AddNodalSolutionStepVariable(MOMENTUM)
    model_part.AddNodalSolutionStepVariable(ROTATION)
    model_part.AddNodalSolutionStepVariable(PRESSURE)        
    model_part.AddNodalSolutionStepVariable(ERROR_RATIO)
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(NODAL_ERROR_1)
    print("variables for the dynamic structural solution added correctly")

def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs 
#        node.AddDof(DISPLACEMENT_X)
#        node.AddDof(DISPLACEMENT_Y)
#        node.AddDof(DISPLACEMENT_Z)
        node.AddDof(DISPLACEMENT_X, REACTION_X)
        node.AddDof(DISPLACEMENT_Y, REACTION_Y)
        node.AddDof(DISPLACEMENT_Z, REACTION_Z)
        node.AddDof(WATER_PRESSURE)
        node.AddDof(AIR_PRESSURE)
#        node.AddDof(LAGRANGE_DISPLACEMENT_X, REACTION_LAGRANGE_DISPLACEMENT_X)
#        node.AddDof(LAGRANGE_DISPLACEMENT_Y, REACTION_LAGRANGE_DISPLACEMENT_Y)
#        node.AddDof(LAGRANGE_DISPLACEMENT_Z, REACTION_LAGRANGE_DISPLACEMENT_Z)
        node.AddDof(LAGRANGE_DISPLACEMENT_X)
        node.AddDof(LAGRANGE_DISPLACEMENT_Y)
        node.AddDof(LAGRANGE_DISPLACEMENT_Z)
        node.AddDof(ROTATION_X)
        node.AddDof(ROTATION_Y)
        node.AddDof(ROTATION_Z)
        #node.AddDof(LAGRANGE_AIR_PRESSURE)
        node.AddDof(LAGRANGE_WATER_PRESSURE)
    print("dofs for the dynamic structural solution added correctly")

class SampleSolver():

    def __init__(self, model_part, abs_tol, rel_tol):
        self.echo_level = 0
        self.toll = rel_tol
        self.absolute_tol = abs_tol
        self.structure_linear_solver = SkylineLUFactorizationSolver()
        self.conv_criteria = DisplacementCriteria(rel_tol, abs_tol)
        self.CalculateReactionFlag = False

    def Initialize(self):
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
#        self.conv_criteria = MultiPhaseFlowCriteria(self.toll, self.absolute_tol)
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = False
        self.MaxNewtonRapshonIterations = 30
        self.parallel_space = UblasSparseSpace()
        self.solver = NonlinearPrescribedDisplacementSolver(self.parallel_space, None, self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.builder_and_solver, self.MaxNewtonRapshonIterations, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag)

class NonlinearPrescribedDisplacementSolver:
    def __init__(self, parallel_space, comm, model_part, time_scheme, linear_solver, conv_criteria, builder_and_solver, MaxNewtonRapshonIterations, CalculateReactionsFlag, ReformDofSetAtEachStep, MoveMeshFlag):
        self.comm = comm
        self.parallel_space = parallel_space
        self.model_part = model_part
        self.time_scheme = time_scheme
        self.linear_solver = linear_solver
        self.conv_criteria = conv_criteria
        self.builder_and_solver = builder_and_solver
        self.max_iter = MaxNewtonRapshonIterations
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.MoveMeshFlag = MoveMeshFlag

        # default values for some variables
        self.echo_level = 1

        # local matrices and vectors
        if(self.comm != None):
            self.pA = self.parallel_space.CreateEmptyMatrixPointer(self.comm)
            self.pDx = self.parallel_space.CreateEmptyVectorPointer(self.comm)
            self.pb = self.parallel_space.CreateEmptyVectorPointer(self.comm)
        else:
            self.pA = self.parallel_space.CreateEmptyMatrixPointer()
            self.pDx = self.parallel_space.CreateEmptyVectorPointer()
            self.pb = self.parallel_space.CreateEmptyVectorPointer()

        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()

        # initialize flags
        self.SolutionStepIsInitialized = False
        self.InitializeWasPerformed = False
        self.StiffnessMatrixIsBuilt = False

        # provide settings to the builder and solver
        (self.builder_and_solver).SetCalculateReactionsFlag(self.CalculateReactionsFlag)
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)

        self.iterations_last_solution = 0

    #
    def Initialize(self):
        if(self.time_scheme.SchemeIsInitialized() == False):
            self.time_scheme.Initialize(self.model_part)

        if (self.time_scheme.ElementsAreInitialized() == False):
            self.time_scheme.InitializeElements(self.model_part)

        if (self.time_scheme.ConditionsAreInitialized() == False):
            self.time_scheme.InitializeConditions(self.model_part)

    def Solve(self):
        self.SolveOneStep()

    #
    def SolveOneStep(self):
        # perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            self.Initialize()
            self.InitializeWasPerformed = True

        # update the prescribed displacement
        for node in self.model_part.Nodes:
            if node.IsFixed(DISPLACEMENT_X) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_X):
                ux = node.GetSolutionStepValue(DISPLACEMENT_X)
                dux = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, ux + dux)
            if node.IsFixed(DISPLACEMENT_Y) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_Y):
                uy = node.GetSolutionStepValue(DISPLACEMENT_Y)
                duy = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y)
                node.SetSolutionStepValue(DISPLACEMENT_Y, uy + duy)
            if node.IsFixed(DISPLACEMENT_Z) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_Z):
                uz = node.GetSolutionStepValue(DISPLACEMENT_Z)
                duz = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z)
                node.SetSolutionStepValue(DISPLACEMENT_Z, uz + duz)

        # perform initializations for the current step
        # this operation implies:
        # identifying the set of DOFs that will be solved during this step
        # organizing the DOFs so to identify the dirichlet conditions
        # resizing the matrix preallocating the "structure"
        if (self.SolutionStepIsInitialized == False):
            if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
                reform_dofs = True
            else:
                reform_dofs = False
            self.InitializeSolutionStep(reform_dofs)
            self.SolutionStepIsInitialized = True

        # perform prediction
        self.Predict()

        # execute iteration - first iteration is ALWAYS executed
        converged = self.ExecuteIteration(self.echo_level, self.MoveMeshFlag)
        it = 1

        # set prescribed delta displacement to zero
        for node in self.model_part.Nodes:
            if node.IsFixed(DISPLACEMENT_X) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_X):
                node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, 0.0)
            if node.IsFixed(DISPLACEMENT_Y) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_Y):
                node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, 0.0)
            if node.IsFixed(DISPLACEMENT_Z) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_Z):
                node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z, 0.0)

        # non linear loop
        while(it < self.max_iter and converged == False):
            # calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            converged = self.ExecuteIteration(self.echo_level, self.MoveMeshFlag)

            # update iteration count
            it = it + 1

        # finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)
        self.SolutionStepIsInitialized = False

        self.iterations_last_solution = it - 1

        # clear if needed - deallocates memory
        mpi.world.barrier()
        if(self.ReformDofSetAtEachStep):
            self.Clear()
        if(mpi.rank == 0):
            print("SolveOneStep is Finished")

        if(it == self.max_iter or converged == False):
            print("Iteration does not converge in " + str(self.max_iter) + " steps")
            print("Solve for displacements failed")
            sys.exit(0)

    #
    def Predict(self):
        self.time_scheme.Predict(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

    #
    def InitializeSolutionStep(self, reform_dofs):
        if(reform_dofs):
            # initialize the list of degrees of freedom to be used
            self.builder_and_solver.SetUpDofSet(self.time_scheme, self.model_part)
            # reorder the list of degrees of freedom to identify fixity and system size
            self.builder_and_solver.SetUpSystem(self.model_part)
            # allocate memory for the system and preallocate the structure of the matrix
            self.builder_and_solver.ResizeAndInitializeVectors(self.pA, self.pDx, self.pb, self.model_part.Elements, self.model_part.Conditions, self.model_part.ProcessInfo)

            # updating references
            self.A = (self.pA).GetReference()
            self.Dx = (self.pDx).GetReference()
            self.b = (self.pb).GetReference()

            # clear scheme so dof update map is recomputed for the new Dof set
            self.time_scheme.Clear()
            
        self.builder_and_solver.InitializeSolutionStep(self.model_part, self.A, self.Dx, self.b)
        self.time_scheme.InitializeSolutionStep(self.model_part, self.A, self.Dx, self.b)

    #
    def ExecuteIteration(self, echo_level, MoveMeshFlag):
        # verify convergence
        converged = self.conv_criteria.PreCriteria(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

        # reset system matrices and vectors prior to rebuild
        self.parallel_space.SetToZeroMatrix(self.A)
        self.parallel_space.SetToZeroVector(self.Dx)
        self.parallel_space.SetToZeroVector(self.b)

        self.time_scheme.InitializeNonLinIteration(self.model_part, self.A, self.Dx, self.b)

        # build and solve the problem
#        self.builder_and_solver.Build(self.time_scheme, self.model_part, self.A, self.b)
#        self.builder_and_solver.ApplyDirichletConditions(self.time_scheme, self.model_part, self.A, self.Dx, self.b)
# #        self.linear_solver.ProvideAdditionalData(self.A, self.Dx, self.b, self.builder_and_solver.GetDofSet(), self.model_part)
#        self.builder_and_solver.SystemSolve(self.A, self.Dx, self.b)
        self.builder_and_solver.BuildAndSolve(self.time_scheme, self.model_part, self.A, self.Dx, self.b)

        # perform update
        self.time_scheme.Update(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b);

        # move the mesh as needed
        if(MoveMeshFlag):
            self.time_scheme.MoveMesh(self.model_part.Nodes);

        self.time_scheme.FinalizeNonLinIteration(self.model_part, self.A, self.Dx, self.b)

        # verify convergence
        converged = self.conv_criteria.PostCriteria(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)
        return converged

    #
    def FinalizeSolutionStep(self, CalculateReactionsFlag):
        if(CalculateReactionsFlag):
            self.builder_and_solver.CalculateReactions(self.time_scheme, self.model_part, self.A, self.Dx, self.b)
        # Finalisation of the solution step,
        self.time_scheme.FinalizeSolutionStep(self.model_part, self.A, self.Dx, self.b)
        self.builder_and_solver.FinalizeSolutionStep(self.model_part, self.A, self.Dx, self.b)
        self.time_scheme.Clean()
        # reset flags for the next step
        self.mSolutionStepIsInitialized = False

    #
    def Clear(self):
        mpi.world.barrier()
        if(mpi.rank == 0):
            print("Entered in Clear")
        self.parallel_space.ClearMatrix(self.pA)
        self.parallel_space.ClearVector(self.pDx)
        self.parallel_space.ClearVector(self.pb)

        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()

        self.builder_and_solver.SetDofSetIsInitializedFlag(False)

        self.builder_and_solver.Clear()

        self.time_scheme.Clear()

        if(mpi.rank == 0):
            print("Clear is completed")

    #
    def SetEchoLevel(self, level):
        self.echo_level = level
        self.builder_and_solver.SetEchoLevel(level)

