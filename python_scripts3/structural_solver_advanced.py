#importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralApplication as StructuralApplication

import KratosMultiphysics.StructuralApplication.structural_solver_static as structural_solver_static
import KratosMultiphysics.StructuralApplication.newton_raphson_strategy as newton_raphson_strategy

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.ELASTIC_BEDDING_STIFFNESS)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FORCE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.DISPLACEMENT_OLD)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.DISPLACEMENT_NULL)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.DISPLACEMENT_EINS)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.DISPLACEMENT_DT)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.DISPLACEMENT_NULL_DT)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.DISPLACEMENT_EINS_DT)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.ACCELERATION_NULL)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.ACCELERATION_EINS)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.ELASTIC_LEFT_CAUCHY_GREEN_OLD)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.INSITU_STRESS)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.PRESTRESS)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRESSES)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRAIN)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_LOAD)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.AIR_PRESSURE_NULL)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.AIR_PRESSURE_EINS)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.AIR_PRESSURE_DT)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.AIR_PRESSURE_NULL_DT)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.AIR_PRESSURE_EINS_DT)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.AIR_PRESSURE_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.AIR_PRESSURE_NULL_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.AIR_PRESSURE_EINS_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.WATER_PRESSURE_NULL)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.WATER_PRESSURE_EINS)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.WATER_PRESSURE_DT)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.WATER_PRESSURE_NULL_DT)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.WATER_PRESSURE_EINS_DT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.WATER_PRESSURE_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.WATER_PRESSURE_NULL_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.WATER_PRESSURE_EINS_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.EXCESS_PORE_WATER_PRESSURE)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.VISCOSITY)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.PRESCRIBED_DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.ELASTIC_STRAIN_VECTOR)
    #auxiliary variables misused for mesh rezoning ;-)
    # model_part.AddNodalSolutionStepVariable(IS_VISITED)
    # model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.LAGRANGE_MULTIPLIER)
#    model_part.AddNodalSolutionStepVariable(GAP)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LAGRANGE_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LAGRANGE_AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LAGRANGE_WATER_PRESSURE)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.LAGRANGE_MULTIPLIER_CONSTRAINT)
    #model_part.AddNodalSolutionStepVariable(INTERNAL_VARIABLES)
    # model_part.AddNodalSolutionStepVariable(MOMENTUM)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ERROR_RATIO)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.NODAL_ERROR_1)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.VON_MISES_STRESS)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.PLASTICITY_INDICATOR)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.THREED_STRESSES)
    model_part.AddNodalSolutionStepVariable(StructuralApplication.INTEGRATION_POINT_STRAIN_VECTOR)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.ROTATION_OLD)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.ROTATION_NULL)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.ROTATION_EINS)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.ROTATION_DT)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.ROTATION_NULL_DT)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.ROTATION_EINS_DT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.ANGULAR_ACCELERATION_NULL)
    # model_part.AddNodalSolutionStepVariable(StructuralApplication.ANGULAR_ACCELERATION_EINS)
    print("variables for the dynamic structural solution added correctly")

def AddDofsForNode(node):
    #adding dofs
#        node.AddDof(DISPLACEMENT_X)
#        node.AddDof(DISPLACEMENT_Y)
#        node.AddDof(DISPLACEMENT_Z)
    node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
    node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
    node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
    node.AddDof(KratosMultiphysics.WATER_PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE)
    node.AddDof(KratosMultiphysics.AIR_PRESSURE, KratosMultiphysics.REACTION_AIR_PRESSURE)
#        node.AddDof(LAGRANGE_DISPLACEMENT_X, REACTION_LAGRANGE_DISPLACEMENT_X)
#        node.AddDof(LAGRANGE_DISPLACEMENT_Y, REACTION_LAGRANGE_DISPLACEMENT_Y)
#        node.AddDof(LAGRANGE_DISPLACEMENT_Z, REACTION_LAGRANGE_DISPLACEMENT_Z)
    node.AddDof(KratosMultiphysics.LAGRANGE_DISPLACEMENT_X)
    node.AddDof(KratosMultiphysics.LAGRANGE_DISPLACEMENT_Y)
    node.AddDof(KratosMultiphysics.LAGRANGE_DISPLACEMENT_Z)
    node.AddDof(KratosMultiphysics.ROTATION_X)
    node.AddDof(KratosMultiphysics.ROTATION_Y)
    node.AddDof(KratosMultiphysics.ROTATION_Z)
    #node.AddDof(LAGRANGE_AIR_PRESSURE)
    node.AddDof(KratosMultiphysics.LAGRANGE_WATER_PRESSURE)

def AddDofsForNodes(nodes):
    for node in nodes:
        AddDofsForNode(node)
    print("dofs for the dynamic structural solution added correctly")

def AddDofs(model_part):
    AddDofsForNodes(model_part.Nodes)

#######################################################################
class SolverAdvanced(structural_solver_static.StaticStructuralSolver):
    def __init__( self, model_part, domain_size, time_steps, analysis_parameters, abs_tol, rel_tol ):
        structural_solver_static.StaticStructuralSolver.__init__( self, model_part, domain_size, abs_tol, rel_tol )
        self.time_steps = time_steps
        self.analysis_parameters = self.CheckAndConvertParameters(analysis_parameters)
        self.echo_level = 0
        self.dissipation_radius = self.analysis_parameters['dissipation_radius']
        self.toll = rel_tol
        self.absolute_tol = abs_tol
        #definition of the solvers
        self.structure_linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        #pDiagPrecond = ParallelDiagonalPreconditioner()
        #self.structure_linear_solver =  ParallelCGSolver(1e-8, 5000,pDiagPrecond)
        #definition of the convergence criteria
        self.conv_criteria = KratosMultiphysics.DisplacementCriteria(0.000001,1e-9)
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
            print('unsupported type of analysis parameters')
            sys.exit(0)


    #######################################################################
    def Initialize(self):
        #definition of time integration scheme
        if( self.analysis_parameters['analysis_type'] == 0 ):
            print("using static scheme")
            self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
            #self.time_scheme = ParallelResidualBasedIncrementalUpdateStaticScheme()
            self.MoveMeshFlag = True
        elif( self.analysis_parameters['analysis_type'] == 1 ):
            print("using newmark quasi-static scheme")
            self.model_part.ProcessInfo.SetValue( StructuralApplication.QUASI_STATIC_ANALYSIS, True )
            if(self.dissipation_radius >= 0.0 and self.dissipation_radius <= 1.0):
                self.time_scheme = StructuralApplication.ResidualBasedNewmarkScheme(self.dissipation_radius)
#                self.time_scheme = ResidualBasedStaggeredNewmarkScheme(self.dissipation_radius)
            else:
                self.time_scheme = StructuralApplication.ResidualBasedNewmarkScheme() #pure Newmarkscheme
#                self.time_scheme = ResidualBasedStaggeredNewmarkScheme() #pure Newmarkscheme
            self.MoveMeshFlag = True
        elif( self.analysis_parameters['analysis_type'] == 2 ):
            print("using newmark dynamic scheme")
            self.model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, False )
            if(self.dissipation_radius >= 0.0 and self.dissipation_radius <= 1.0):
                self.time_scheme = StructuralApplication.ResidualBasedNewmarkScheme(self.dissipation_radius)
            else:
                self.time_scheme = StructuralApplication.ResidualBasedNewmarkScheme() #pure Newmarkscheme
            #self.time_scheme = ResidualBasedPredictorCorrectorBossakScheme(self.dissipation_radius)
            #self.time_scheme.Check(self.model_part)
        else:
            print("analysis type is not defined! Define in analysis_parameters['analysis_type']:")
            print("   'using static scheme': static analysis")
            print("   'using newmark quasi-static scheme': quasi-static analysis")
            print("   'using newmark dynamic scheme': dynamic analysis")
            sys.exit(0)
        #definition of the convergence criteria
        self.conv_criteria = StructuralApplication.MultiPhaseFlowCriteria(self.toll,self.absolute_tol)
        #self.conv_criteria = MultiPhaseFlowCriteria(1.0e-13,1.0e-13)
        #self.conv_criteria = ResidualBasedMultiPhaseCriteria(self.toll,self.absolute_tol)
        #self.conv_criteria = ResidualCriteria(1.0e-9,1.0e-9)
#        self.conv_criteria = DisplacementCriteria(self.toll,self.absolute_tol)
        if(self.analysis_parameters['decouple_build_and_solve'] == False):
            if(self.analysis_parameters['builder_and_solver_type'] == "residual-based elimination deactivation"):
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block"):
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.structure_linear_solver)
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block with constraints"):
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolverWithConstraints(self.structure_linear_solver)
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block with constraints element-wise"):
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolverWithConstraintsElementWise(self.structure_linear_solver)
            #builder_and_solver = MultiPhaseBuilderAndSolver(self.structure_linear_solver)
            #builder_and_solver = BuiMultiPhaseBuilderAndSolver(self.structure_linear_solver)
            #builder_and_solver = ParallelResidualBasedEliminationBuilderAndSolverDeactivation(self.structure_linear_solver)
            #builder_and_solver = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block with constraints deactivation"):
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivation(self.structure_linear_solver)
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block with constraints deactivation element-wise"):
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivationElementWise(self.structure_linear_solver)
        else:
            if(self.analysis_parameters['builder_and_solver_type'] == "residual-based elimination deactivation"):
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(KratosMultiphysics.LinearSolver())
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block"):
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(KratosMultiphysics.LinearSolver())
            elif(self.analysis_parameters['builder_and_solver_type'] == "residual-based block with constraints"):
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolverWithConstraints(KratosMultiphysics.LinearSolver())
        print("builder_and_solver type: " + str(self.analysis_parameters['builder_and_solver_type']))

        #creating the solution strategy
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = True
        self.space_utils = KratosMultiphysics.UblasSparseSpace()
        self.model_part.ProcessInfo[StructuralApplication.RESET_CONFIGURATION] = 0
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
