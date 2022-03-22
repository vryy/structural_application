#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
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
    model_part.AddNodalSolutionStepVariable(FACE_LOAD_NULL)
    model_part.AddNodalSolutionStepVariable(FACE_LOAD_EINS)
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
    model_part.AddNodalSolutionStepVariable(PRESCRIBED_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(ELASTIC_STRAIN_VECTOR)
    #auxiliary variables misused for mesh rezoning ;-)
    model_part.AddNodalSolutionStepVariable(IS_VISITED)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_MULTIPLIER)
#    model_part.AddNodalSolutionStepVariable(GAP)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(REACTION_LAGRANGE_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_MULTIPLIER_CONSTRAINT)
    #model_part.AddNodalSolutionStepVariable(INTERNAL_VARIABLES)
    model_part.AddNodalSolutionStepVariable(MOMENTUM)
    model_part.AddNodalSolutionStepVariable(MOMENT)
    model_part.AddNodalSolutionStepVariable(ROTATION)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(ERROR_RATIO)
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(NODAL_ERROR_1)
    model_part.AddNodalSolutionStepVariable(VON_MISES_STRESS)
    model_part.AddNodalSolutionStepVariable(PLASTICITY_INDICATOR)
    model_part.AddNodalSolutionStepVariable(THREED_STRESSES)
    model_part.AddNodalSolutionStepVariable(INTEGRATION_POINT_STRAIN_VECTOR)
    model_part.AddNodalSolutionStepVariable(ROTATION_OLD)
    model_part.AddNodalSolutionStepVariable(ROTATION_NULL)
    model_part.AddNodalSolutionStepVariable(ROTATION_EINS)
    model_part.AddNodalSolutionStepVariable(ROTATION_DT)
    model_part.AddNodalSolutionStepVariable(ROTATION_NULL_DT)
    model_part.AddNodalSolutionStepVariable(ROTATION_EINS_DT)
    model_part.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION_NULL)
    model_part.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION_EINS)
    print("variables for the dynamic structural solution added correctly")

def AddDofsForNode(node):
    #adding dofs
#        node.AddDof(DISPLACEMENT_X)
#        node.AddDof(DISPLACEMENT_Y)
#        node.AddDof(DISPLACEMENT_Z)
    node.AddDof(DISPLACEMENT_X, REACTION_X)
    node.AddDof(DISPLACEMENT_Y, REACTION_Y)
    node.AddDof(DISPLACEMENT_Z, REACTION_Z)
    node.AddDof(WATER_PRESSURE, REACTION_WATER_PRESSURE)
    node.AddDof(AIR_PRESSURE, REACTION_AIR_PRESSURE)
#        node.AddDof(LAGRANGE_DISPLACEMENT_X, REACTION_LAGRANGE_DISPLACEMENT_X)
#        node.AddDof(LAGRANGE_DISPLACEMENT_Y, REACTION_LAGRANGE_DISPLACEMENT_Y)
#        node.AddDof(LAGRANGE_DISPLACEMENT_Z, REACTION_LAGRANGE_DISPLACEMENT_Z)
    node.AddDof(LAGRANGE_DISPLACEMENT_X)
    node.AddDof(LAGRANGE_DISPLACEMENT_Y)
    node.AddDof(LAGRANGE_DISPLACEMENT_Z)
    node.AddDof(ROTATION_X)
    node.AddDof(ROTATION_Y)
    node.AddDof(ROTATION_Z)
    # node.AddDof(LAGRANGE_AIR_PRESSURE, REACTION_LAGRANGE_AIR_PRESSURE)
    # node.AddDof(LAGRANGE_WATER_PRESSURE, REACTION_LAGRANGE_WATER_PRESSURE) # do not add here, add in the include file instead

def AddDofsForNodes(nodes):
    for node in nodes:
        AddDofsForNode(node)
    print("dofs for the dynamic structural solution added correctly")

def AddDofs(model_part):
    AddDofsForNodes(model_part.Nodes)

#######################################################################
class SolverAdvanced(structural_solver_static.StaticStructuralSolver):
    def __init__( self, model_part, domain_size, time_steps, analysis_parameters, abs_tol, rel_tol ):
        structural_solver_static.StaticStructuralSolver.__init__( self, model_part, domain_size )
        self.time_steps = time_steps
        self.analysis_parameters = self.CheckAndConvertParameters(analysis_parameters)
        self.echo_level = 0
        self.dissipation_radius = self.analysis_parameters['dissipation_radius']
        self.toll = rel_tol
        self.absolute_tol = abs_tol
        #definition of the solvers
        self.structure_linear_solver =  SkylineLUFactorizationSolver()
        #pDiagPrecond = ParallelDiagonalPreconditioner()
        #self.structure_linear_solver =  ParallelCGSolver(1e-8, 5000,pDiagPrecond)
        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(0.000001,1e-9)
        #self.conv_criteria = ParallelDisplacementCriteria(0.000001,1e-9)
        self.CalculateReactionFlag = False

    #######################################################################
    def CheckAndConvertParameters(self, analysis_parameters):
        if( type( analysis_parameters ) == dict ):
            if 'builder_and_solver_type' not in analysis_parameters:
                analysis_parameters['builder_and_solver_type'] = "residual-based elimination deactivation"
            return analysis_parameters
        elif( type( analysis_parameters ) == list ):
            new_analysis_parameters = {}
            new_analysis_parameters['perform_contact_analysis_flag'] = analysis_parameters[0]
            new_analysis_parameters['penalty'] = analysis_parameters[1]
            new_analysis_parameters['maxuzawa'] = analysis_parameters[2]
            new_analysis_parameters['friction'] = analysis_parameters[3]
            new_analysis_parameters['frictionpenalty'] = analysis_parameters[4]
            new_analysis_parameters['contact_double_check_flag'] = analysis_parameters[5]
            new_analysis_parameters['contact_ramp_penalties_flag'] = analysis_parameters[6]
            new_analysis_parameters['maxpenalty'] = analysis_parameters[7]
            new_analysis_parameters['rampcriterion'] = analysis_parameters[8]
            new_analysis_parameters['rampfactor'] = analysis_parameters[9]
            new_analysis_parameters['fricmaxpenalty'] = analysis_parameters[10]
            new_analysis_parameters['fricrampcriterion'] = analysis_parameters[11]
            new_analysis_parameters['fricrampfactor'] = analysis_parameters[12]
            new_analysis_parameters['print_sparsity_info_flag'] = analysis_parameters[13]
            new_analysis_parameters['analysis_type'] = analysis_parameters[14]
            if(len(analysis_parameters) > 15):
                new_analysis_parameters['dissipation_radius'] = analysis_parameters[15]
            else:
                if new_analysis_parameters['analysis_type'] == 2:
                    new_analysis_parameters['dissipation_radius'] = 0.1
                else:
                    new_analysis_parameters['dissipation_radius'] = 1.0
            new_analysis_parameters['decouple_build_and_solve'] = False
            new_analysis_parameters['builder_and_solver_type'] = "residual-based elimination deactivation"
            return new_analysis_parameters
        else:
            print 'unsupported type of analysis parameters'
            sys.exit(0)


    #######################################################################
    def Initialize(self):
        #definition of time integration scheme
        if( self.analysis_parameters['analysis_type'] == 0 ):
            print("using static scheme")
            # self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
            self.time_scheme = ResidualBasedIncrementalUpdateStaticDeactivationScheme()
            #self.time_scheme = ParallelResidualBasedIncrementalUpdateStaticScheme()
            self.MoveMeshFlag = True
        elif( self.analysis_parameters['analysis_type'] == 1 ):
            print("using newmark quasi-static scheme, dissipation_radius=" + str(self.dissipation_radius))
            self.model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, True )
            if(self.dissipation_radius >= 0.0 and self.dissipation_radius <= 1.0): # generalized Newmark-alpha
                self.time_scheme = ResidualBasedNewmarkScheme(self.dissipation_radius)
            elif(self.dissipation_radius >= 2.0 and self.dissipation_radius <= 3.0): # Bossak-alpha
                self.time_scheme = ResidualBasedNewmarkScheme(2, self.dissipation_radius-2.0)
            elif(self.dissipation_radius >= 4.0 and self.dissipation_radius <= 5.0): # Hilber-alpha
                self.time_scheme = ResidualBasedNewmarkScheme(3, self.dissipation_radius-4.0)
            else:
                self.time_scheme = ResidualBasedNewmarkScheme() # pure Newmarkscheme
            self.MoveMeshFlag = True
        elif( self.analysis_parameters['analysis_type'] == 2 ):
            print("using newmark dynamic scheme, dissipation_radius=" + str(self.dissipation_radius))
            self.model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, False )
            if(self.dissipation_radius >= 0.0 and self.dissipation_radius <= 1.0): # generalized Newmark-alpha
                self.time_scheme = ResidualBasedNewmarkScheme(self.dissipation_radius)
            elif(self.dissipation_radius >= 2.0 and self.dissipation_radius <= 3.0): # Bossak-alpha
                self.time_scheme = ResidualBasedNewmarkScheme(2, self.dissipation_radius-2.0)
            elif(self.dissipation_radius >= 4.0 and self.dissipation_radius <= 5.0): # Hilber-alpha
                self.time_scheme = ResidualBasedNewmarkScheme(3, self.dissipation_radius-4.0)
            else:
                self.time_scheme = ResidualBasedNewmarkScheme() # pure Newmarkscheme
            #self.time_scheme.Check(self.model_part)
        elif( self.analysis_parameters['analysis_type'] == 3 ):
            self.model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, True )
            print("using theta quasi-static scheme, theta=" + str(self.dissipation_radius))
            self.time_scheme = ResidualBasedThetaScheme(self.dissipation_radius)
        elif( self.analysis_parameters['analysis_type'] == 4 ):
            self.model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, False )
            print("using theta dynamics scheme, theta=" + str(self.dissipation_radius))
            self.time_scheme = ResidualBasedThetaScheme(self.dissipation_radius)
        else:
            print("analysis type is not defined! Define in analysis_parameters['analysis_type']:")
            print("   'using static scheme': static analysis")
            print("   'using newmark quasi-static scheme': quasi-static analysis")
            print("   'using newmark dynamic scheme': dynamic analysis")
            sys.exit(0)
        #definition of the convergence criteria
        self.conv_criteria = MultiPhaseFlowCriteria(self.toll,self.absolute_tol)
        #self.conv_criteria = MultiPhaseFlowCriteria(1.0e-13,1.0e-13)
        #self.conv_criteria = ResidualBasedMultiPhaseCriteria(self.toll,self.absolute_tol)
        #self.conv_criteria = ResidualCriteria(1.0e-9,1.0e-9)
#        self.conv_criteria = DisplacementCriteria(self.toll,self.absolute_tol)
        if(self.analysis_parameters['decouple_build_and_solve'] == False):
            if(self.analysis_parameters['builder_and_solver_type'] == "residual-based elimination deactivation"):
                builder_and_solver = ResidualBasedEliminationBuilderAndSolverDeactivation(self.structure_linear_solver)
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block"):
                builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.structure_linear_solver)
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block with constraints"):
                builder_and_solver = ResidualBasedBlockBuilderAndSolverWithConstraints(self.structure_linear_solver)
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block with constraints element-wise"):
                builder_and_solver = ResidualBasedBlockBuilderAndSolverWithConstraintsElementWise(self.structure_linear_solver)
            #builder_and_solver = MultiPhaseBuilderAndSolver(self.structure_linear_solver)
            #builder_and_solver = BuiMultiPhaseBuilderAndSolver(self.structure_linear_solver)
            #builder_and_solver = ParallelResidualBasedEliminationBuilderAndSolverDeactivation(self.structure_linear_solver)
            #builder_and_solver = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block with constraints deactivation"):
                builder_and_solver = ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivation(self.structure_linear_solver)
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block with constraints deactivation element-wise"):
                builder_and_solver = ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivationElementWise(self.structure_linear_solver)
        else:
            if(self.analysis_parameters['builder_and_solver_type'] == "residual-based elimination deactivation"):
                builder_and_solver = ResidualBasedEliminationBuilderAndSolverDeactivation(LinearSolver())
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block"):
                builder_and_solver = ResidualBasedBlockBuilderAndSolver(LinearSolver())
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block with constraints"):
                builder_and_solver = ResidualBasedBlockBuilderAndSolverWithConstraints(LinearSolver())
        print("builder_and_solver type: " + str(self.analysis_parameters['builder_and_solver_type']))

        #creating the solution strategy
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = True
        self.space_utils = UblasSparseSpace()
        self.model_part.ProcessInfo[RESET_CONFIGURATION] = 0

        if self.analysis_parameters['perform_contact_analysis_flag'] == True:
            import uzawa_contact_strategy
            self.solver = uzawa_contact_strategy.SolvingStrategyPython( self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, builder_and_solver )
        else:
            if 'solution_strategy' in self.analysis_parameters:
                if self.analysis_parameters['solution_strategy'] == "implicit_Newton_Raphson":
                    import newton_raphson_strategy
                    self.solver = newton_raphson_strategy.SolvingStrategyPython( self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, builder_and_solver )
                elif self.analysis_parameters['solution_strategy'] == "initial_stiffness":
                    import initial_stiffness_strategy
                    self.solver = initial_stiffness_strategy.SolvingStrategyPython( self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, builder_and_solver )
                elif self.analysis_parameters['solution_strategy'] == "modified_thomas":
                    import modified_thomas_strategy
                    self.solver = modified_thomas_strategy.SolvingStrategyPython( self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, builder_and_solver )
                elif self.analysis_parameters['solution_strategy'] == "modified_sloan":
                    import modified_sloan_strategy
                    self.solver = modified_sloan_strategy.SolvingStrategyPython( self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, builder_and_solver )
                elif self.analysis_parameters['solution_strategy'] == "arc_length_load_control":
                    import arc_length_load_control_strategy
                    self.solver = arc_length_load_control_strategy.SolvingStrategyPython( self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, builder_and_solver )
                elif self.analysis_parameters['solution_strategy'] == "arc_length_displacement_control":
                    if( self.analysis_parameters['analysis_type'] == 0 ):
                        self.time_scheme = ArcLengthDisplacementControlResidualBasedIncrementalUpdateStaticDeactivationScheme(self.time_scheme)
                    else:
                        raise Exception("analysis_type > 0 is not yet supported for arc-length displacement control")
                    import arc_length_displacement_control_strategy
                    self.solver = arc_length_displacement_control_strategy.SolvingStrategyPython( self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, builder_and_solver )
                else:
                    raise Exception("Unknown solution_strategy " + str(self.analysis_parameters['solution_strategy']))
            else:
                import newton_raphson_strategy
                self.solver = newton_raphson_strategy.SolvingStrategyPython( self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, builder_and_solver )

    #######################################################################
    def InitializeSolver(self):
        self.solver.Initialize()

    #######################################################################
    def SolveLagrange(self):
        (self.solver).SolveLagrange()

    #######################################################################
    def SolveLagrangeLocal(self):
        (self.solver).SolveLagrangeLocal()

    #######################################################################
    def SolveOneStep(self):
        (self.solver).SolveOneStep()

