import KratosMultiphysics
import KratosMultiphysics.StructuralApplication as StructuralApplication

def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_LOAD)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.PRESCRIBED_DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LAGRANGE_DISPLACEMENT)
    print("variables for the static structural solution added correctly")

def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
        node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
        node.AddDof(KratosMultiphysics.LAGRANGE_DISPLACEMENT_X)
        node.AddDof(KratosMultiphysics.LAGRANGE_DISPLACEMENT_Y)
        node.AddDof(KratosMultiphysics.LAGRANGE_DISPLACEMENT_Z)
    print("dofs for the static structural solution added correctly")

class StaticStructuralSolver:
    #

    def __init__(self, model_part, domain_size, abs_tol, rel_tol):

        self.model_part = model_part
        self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        # if called here, Check may be called before the system is completely set up!!!
        self.time_scheme.Check(self.model_part)

        # definition of the solvers
        self.structure_linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()

        # definition of the convergence criteria
        self.conv_criteria = KratosMultiphysics.DisplacementCriteria(rel_tol, abs_tol)
        self.conv_criteria.Check(self.model_part)
        self.MaxNewtonRapshonIterations = 100

        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = False

    #
    def Initialize(self):

        # creating the solution strategy
        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.MaxNewtonRapshonIterations, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag)

    #
    def Check(self):
        self.solver.Check()

    #
    def Solve(self):
        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

# def CreateSolver(model_part, config):
#     solver = StaticStructuralSolver(model_part, config.domain_size)

#     import linear_solver_factory
#     if(hasattr(config, "linear_solver_config")):
#         solver.structure_linear_solver = linear_solver_factory.ConstructSolver(config.linear_solver_config)

#     return solver
