#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MultithreadedSolversApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

import sys

import structural_solver_static

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
    print("variables for the modal structural solution added correctly")
    
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs 
        node.AddDof(DISPLACEMENT_X, REACTION_X)
        node.AddDof(DISPLACEMENT_Y, REACTION_Y)
        node.AddDof(DISPLACEMENT_Z, REACTION_Z)
    print("dofs for the modal structural solution added correctly")
        
#######################################################################
class ModalSolver():
    def __init__( self, model_part, time_scheme, structure_eigen_solver, analysis_parameters, parallel_space, builder_and_solver):
        self.echo_level = 0
        self.model_part = model_part
        self.structure_eigen_solver = structure_eigen_solver
        self.time_scheme = time_scheme
        self.parallel_space = parallel_space
        self.builder_and_solver = builder_and_solver
        self.analysis_parameters = analysis_parameters

        if 'number_of_required_eigenvalues' not in self.analysis_parameters:
            self.analysis_parameters['number_of_required_eigenvalues'] = 1
        if 'maximum_number_of_iterations' not in self.analysis_parameters:
            self.analysis_parameters['maximum_number_of_iterations'] = 1000
        if 'tolerance' not in self.analysis_parameters:
            self.analysis_parameters['tolerance'] = 1.0e-10
        if 'export_eigensystem' not in self.analysis_parameters:
            self.analysis_parameters['export_eigensystem'] = False
        if 'eigensystem_solver' not in self.analysis_parameters:
            self.analysis_parameters['eigensystem_solver'] = 'None'

        if 'feast_number_of_requested_eigenvalues' not in self.analysis_parameters:
            self.analysis_parameters['feast_number_of_requested_eigenvalues'] = 10
        if 'feast_lower_bound' not in self.analysis_parameters:
            self.analysis_parameters['feast_lower_bound'] = 0.0
        if 'feast_upper_bound' not in self.analysis_parameters:
            self.analysis_parameters['feast_upper_bound'] = 1.0e2

    #######################################################################
    def Initialize(self):
        #definition of time integration scheme
        self.dof_util = DofUtility()

    #######################################################################   
    def SetEchoLevel(self, level):
        self.echo_level = level
        self.builder_and_solver.SetEchoLevel(level)

    #######################################################################   
    def Solve(self):
        #local matrices and vectors
        self.pK = self.parallel_space.CreateEmptyMatrixPointer()
        self.pM = self.parallel_space.CreateEmptyMatrixPointer()

        self.K = (self.pK).GetReference()
        self.M = (self.pM).GetReference()

        #initialize the list of degrees of freedom to be used 
        self.builder_and_solver.SetUpDofSet(self.time_scheme,self.model_part);

        #reorder the list of degrees of freedom to identify fixity and system size	  			
        self.builder_and_solver.SetUpSystem(self.model_part)

        #allocate memory for the system and preallocate the structure of the matrix
        self.builder_and_solver.ResizeAndInitializeEigenSystem(self.pK,self.pM,self.model_part.Elements,self.model_part.Conditions,self.model_part.ProcessInfo)

        #build the eigensystem
        self.builder_and_solver.BuildEigenSystem(self.time_scheme, self.model_part, self.K, self.M)
        self.dof_util.ListDofs(self.builder_and_solver.GetDofSet(),self.builder_and_solver.GetEquationSystemSize())

        # export the matrices
        if self.analysis_parameters['export_eigensystem'] == True:
            self.parallel_space.WriteMatrixMarketMatrix("K.mm", self.K, False)
            self.parallel_space.WriteMatrixMarketMatrix("M.mm", self.M, False)
            print("Export Eigensystem completed")

        #solve for eigenvalues
        if self.analysis_parameters['eigensystem_solver'] == "Feast":
            eigen_solver = FeastSolver()
            self.E = Vector()
            self.V = []
#            eigen_solver.Solve(self.K, 100, V, 0.0, 1.0e2)
            ne = self.analysis_parameters['feast_number_of_requested_eigenvalues']
            lb = self.analysis_parameters['feast_lower_bound']
            ub = self.analysis_parameters['feast_upper_bound']
            eigen_solver.SolveGeneralized(self.K, self.M, ne, self.E, self.V, lb, ub)
 #            return [E, V] # this will not work

        elif self.analysis_parameters['eigensystem_solver'] == "Arpack":
#            for now it doesn't work
            eigen_solver = ArpackSolver()
            V = Vector();
            eigen_solver.Solve(self.K, 5, V)
            print V
        else:
            pass

    #######################################################################   
    def Update(self, V):
        # set the displacement in the system to zero
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)
            node.SetSolutionStepValue(DISPLACEMENT_Z, 0.0)

        # update the solution vector
        Adummy = CompressedMatrix()
        bdummy = Vector()
        self.time_scheme.Update(self.model_part,self.builder_and_solver.GetDofSet(),Adummy,V,bdummy);

    #######################################################################   
    def Clear(self):
        #clear the system
        self.parallel_space.ClearMatrix(self.pK)
        self.parallel_space.ClearMatrix(self.pM)

        #updating references
        self.K = (self.pK).GetReference()
        self.M = (self.pM).GetReference()

        self.builder_and_solver.Clear()
        self.time_scheme.Clear()

