from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
CheckForPreviousImport()


def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL)
    model_part.AddNodalSolutionStepVariable(ACCELERATION_EINS)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(ROTATION)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    # model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(INSITU_STRESS)
    model_part.AddNodalSolutionStepVariable(FACE_LOAD)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(PRESCRIBED_DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_DISPLACEMENT)

    print("variables for the dynamic structural solution added correctly")


def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(DISPLACEMENT_X, REACTION_X)
        node.AddDof(DISPLACEMENT_Y, REACTION_Y)
        node.AddDof(DISPLACEMENT_Z, REACTION_Z)
        node.AddDof(LAGRANGE_DISPLACEMENT_X)
        node.AddDof(LAGRANGE_DISPLACEMENT_Y)
        node.AddDof(LAGRANGE_DISPLACEMENT_Z)
    print("dofs for the dynamic structural solution added correctly")


class StaticStructuralSolver:
    #

    def __init__(self, model_part, domain_size, abs_tol=1e-9, rel_tol=1e-6):

        self.domain_size = domain_size
        self.model_part = model_part

        if self.model_part.Type == "ModelPart":
            self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
            # if called here, Check may be called before the system is completely set up!!!
            # self.time_scheme.Check(self.model_part) # Check shall be called in Initialize

            # definition of the solvers
            self.structure_linear_solver = SkylineLUFactorizationSolver()

            # definition of the convergence criteria
            self.conv_criteria = DisplacementCriteria(rel_tol, abs_tol)
            # self.conv_criteria.Check(self.model_part)

        elif self.model_part.Type == "ComplexModelPart":
            self.time_scheme = ComplexResidualBasedIncrementalUpdateStaticScheme()
            # if called here, Check may be called before the system is completely set up!!!
            # self.time_scheme.Check(self.model_part) # Check shall be called in Initialize

            # definition of the solvers
            self.structure_linear_solver = ComplexSkylineLUFactorizationSolver()

            # definition of the convergence criteria
            self.conv_criteria = ComplexDisplacementCriteria(rel_tol, abs_tol)
            # self.conv_criteria.Check(self.model_part)

        elif self.model_part.Type == "GComplexModelPart":
            self.time_scheme = GComplexResidualBasedIncrementalUpdateStaticScheme()
            # if called here, Check may be called before the system is completely set up!!!
            # self.time_scheme.Check(self.model_part) # Check shall be called in Initialize

            # definition of the solvers
            self.structure_linear_solver = GComplexSkylineLUFactorizationSolver()

            # definition of the convergence criteria
            self.conv_criteria = GComplexDisplacementCriteria(rel_tol, abs_tol)
            # self.conv_criteria.Check(self.model_part)

        else:

            raise Exception("Unknown " + self.model_part.Type)

        self.MaxNewtonRapshonIterations = 100

        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = False
        self.MoveMeshFlag = True

    #
    def Initialize(self):
        # creating the solution strategy

       # import strategy_python
       # self.solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag)
        # (self.solver).SetEchoLevel(2)

        # creating the solution strategy
        if self.model_part.Type == "ModelPart":
            self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.MaxNewtonRapshonIterations, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag)
        elif self.model_part.Type == "ComplexModelPart":
            self.solver = ComplexResidualBasedNewtonRaphsonStrategy(self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.MaxNewtonRapshonIterations, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag)
        elif self.model_part.Type == "GComplexModelPart":
            self.solver = GComplexResidualBasedNewtonRaphsonStrategy(self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.MaxNewtonRapshonIterations, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag)
        # self.solver.Check()

        #(self.solver).SetReformDofSetAtEachStepFlag(True)
        #(self.solver).SetMoveMeshFlag(True)

    #
    def Solve(self):
        return (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

def CreateSolver(model_part, config):
    solver = StaticStructuralSolver(model_part, config.domain_size)

    import linear_solver_factory
    if(hasattr(config, "linear_solver_config")):
        solver.structure_linear_solver = linear_solver_factory.ConstructSolver(config.linear_solver_config)

    return solver
