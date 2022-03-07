#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()
import os,time,math

# Reference:
# Sloan et al, Accelerated initial sti!ness schemes for elastoplasticity

class SolvingStrategyPython:
    #######################################################################
    def __init__( self, model_part, time_scheme, linear_solver, convergence_criteria, CalculateReactionsFlag, ReformDofSetAtEachStep, MoveMeshFlag, Parameters, space_utils, builder_and_solver ):
        #save the input parameters
        self.model_part = model_part
        self.time_scheme = time_scheme
        self.linear_solver = linear_solver
        self.convergence_criteria = convergence_criteria
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.MoveMeshFlag = MoveMeshFlag
        self.Parameters = Parameters
        self.PrintSparsity = self.Parameters['print_sparsity_info_flag']
        self.space_utils = space_utils

        #default values for some variables
        self.max_iter = 30
        self.echo_level = 1
        self.builder_and_solver = builder_and_solver
        self.dof_util = DofUtility()
        self.output_util = OutputUtility()
        self.calculate_reaction_process = None

        #local matrices and vectors
        self.pA = self.space_utils.CreateEmptyMatrixPointer()
        self.pDx = self.space_utils.CreateEmptyVectorPointer()
        self.pDxe = self.space_utils.CreateEmptyVectorPointer()
        self.pb = self.space_utils.CreateEmptyVectorPointer()
        self.pbold = self.space_utils.CreateEmptyVectorPointer()

        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.Dxe = (self.pDxe).GetReference()
        self.b = (self.pb).GetReference()
        ##local matrices and vectors
        #self.A = CompressedMatrix()
        #self.Dx = Vector()
        #self.b = Vector()

        #initialize flags
        self.SolutionStepIsInitialized = False
        self.InitializeWasPerformed = False
        self.StiffnessMatrixIsBuilt = False
        #provide settings to the builder and solver
        (self.builder_and_solver).SetCalculateReactionsFlag(self.CalculateReactionsFlag)
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)

        self.solveCounter = 0 #hbui added this variable

        self.system_reorderer = Process()
#        self.system_reorderer = SystemRCMReordererProcess(self.model_part)
#        self.system_reorderer = SystemBoostRCMReordererProcess(self.model_part)
#        self.system_reorderer = SystemMetisReordererProcess(self.model_part)
#        self.system_reorderer = SystemAMDReordererProcess(self.model_part)

        self.attached_processes = []

        # turn off the calculate reaction flag
        self.model_part.ProcessInfo[SET_CALCULATE_REACTION] = False

        if 'stop_Newton_Raphson_if_not_converge' in self.Parameters:
            self.Parameters['stop_Newton_Raphson_if_not_converged'] = self.Parameters['stop_Newton_Raphson_if_not_converge']

        if 'list_plastic_points' not in self.Parameters:
            self.Parameters['list_plastic_points'] = False

        if not('residuum_tolerance' in self.Parameters):
            self.erbar = 1.0e-4
        else:
            self.erbar = self.Parameters['residuum_tolerance']

        if not('log_residuum_name' in self.Parameters):
            self.log_residuum = open('residuum_modified_thomas.log', 'w')
        else:
            self.log_residuum = open(self.Parameters['log_residuum_name'], 'w')

    def __del__(self):
        self.log_residuum.close()

    #######################################################################
    def Initialize(self):
        if(self.time_scheme.SchemeIsInitialized() == False):
            self.time_scheme.Initialize(self.model_part)
        if (self.time_scheme.ElementsAreInitialized() == False):
            self.time_scheme.InitializeElements(self.model_part)
        if (self.time_scheme.ConditionsAreInitialized() == False):
            self.time_scheme.InitializeConditions(self.model_part)
        for proc in self.attached_processes:
            proc.ExecuteInitialize()
        self.InitializeWasPerformed = True
        print("modified_thomas_strategy.Initialize is called")

    #######################################################################
    def Solve(self):
        #print self.model_part
        self.solveCounter = self.solveCounter + 1
        #solve the nonlinear equilibrium
        self.PerformNewtonRaphsonIteration()
        #finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)
        #clear if needed - deallocates memory
        if(self.ReformDofSetAtEachStep == True):
            self.Clear()

    #######################################################################
    def PerformNewtonRaphsonIteration( self ):
        print("time = " + str(self.model_part.ProcessInfo[TIME]))
        #perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            self.Initialize()
        #perform initializations for the current step
        #this operation implies:
        #identifying the set of DOFs that will be solved during this step
        #organizing the DOFs so to identify the dirichlet conditions
        #resizing the matrix preallocating the "structure"
        if (self.SolutionStepIsInitialized == False):
            self.InitializeSolutionStep()
            self.SolutionStepIsInitialized = True
        #perform prediction
        self.Predict()

        #execute iteration - first iteration is ALWAYS executed, with alpha = 1.0
        calculate_norm = False
        self.iterationCounter = 0 #hbui added this variable
        self.iterationCounter = self.iterationCounter + 1
        alpha = 1.0
        normDx = self.ExecuteIteration(self.echo_level,calculate_norm,0,alpha)
        self.FinalizeNonLinIteration(False,self.MoveMeshFlag)
        er_0 = self.space_utils.TwoNorm(self.b)
        self.space_utils.Copy(self.Dx, self.Dxe)
        normDxe = self.space_utils.TwoNorm(self.Dxe)
        if normDxe > 0.0:
            normDx = self.ExecuteIteration(self.echo_level,calculate_norm,0,1.0)
            self.FinalizeNonLinIteration(False,self.MoveMeshFlag)

        er_n = er_0
        self.log_residuum.write('time: ' + str(self.model_part.ProcessInfo[TIME]) + '\n')
        self.log_residuum.write('it\talpha\tresidual\tratio\treduction\tremaining\n')
        self.log_residuum.write('0\t' + str(alpha) + '\t' + str(er_0) + '\n')
        self.log_residuum.flush()

        #non linear loop
        converged = False
        it = 0
        while(it < self.max_iter and converged == False):
            #verify convergence
            converged = self.convergence_criteria.PreCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

            # compute new alpha parameter
            if normDxe > 0.0:
                alpha += self.space_utils.Dot(self.Dxe, self.Dx) / (normDxe**2)

            #calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            self.iterationCounter = self.iterationCounter + 1
            normDx = self.ExecuteIteration(self.echo_level,calculate_norm,it+1,alpha)
            self.FinalizeNonLinIteration(False,self.MoveMeshFlag)
            self.space_utils.Copy(self.Dx, self.Dxe)
            normDxe = self.space_utils.TwoNorm(self.Dxe)
            if normDxe > 0.0:
                normDx = self.ExecuteIteration(self.echo_level,calculate_norm,it+1,1.0)
                # print("normDx at iteration " + str(it+1) + ": " + str(normDx))

                #verify convergence
                converged = self.convergence_criteria.PostCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

                # finalize the iteration
                self.FinalizeNonLinIteration(converged,self.MoveMeshFlag)
            else:
                print("modified_thomas_strategy.PerformNewtonRaphsonIteration converged since the increment is zero")
                converged = True

            #update iteration count
            it = it + 1

            # estimate the number of remaining iterations
            er = self.space_utils.TwoNorm(self.b)
            # if not(self.alpha_m == 1.0):
            #     n = (it*math.log(er) - math.log(erbar)) / math.log(abs(self.alpha_m))
            # else:
            #     n = self.max_iter - it
            # print("Residual: " + str(er))
            # print("Estimated number of remaining iterations: " + str(n))
            if er_0 > 0.0:
                er_ratio = er/er_0
            else:
                er_ratio = 1.0
            if er > 0.0:
                er_reduction = er_n/er
            else:
                if er_n == 0.0:
                    er_reduction = 1.0
                else:
                    er_reduction = 1.0e99
            er_n = er

            if not(er_reduction == 1.0):
                # n = (it*math.log(er_ratio) - math.log(erbar)) / math.log(er_reduction)
                n = -(math.log(self.erbar) - math.log(er_ratio)) / math.log(er_reduction)
            else:
                n = self.max_iter - it
            # print("Estimated number of remaining iterations: " + str(n))

            # override the convergence
            if (converged == False)  and (er_ratio < self.erbar):
                converged = True
                print("modified_thomas_strategy.PerformNewtonRaphsonIteration converged when residuum ratio (" + str(er_ratio) + ") reached tolerance " + str(self.erbar))

            self.log_residuum.write(str(it) + '\t' + str(alpha) + '\t' + str(er) + '\t' + str(er_ratio) + '\t' + str(er_reduction) + '\t' + str(n) + '\n')
            self.log_residuum.flush()

        self.log_residuum.write("------------------------------------------\n")

        if( it == self.max_iter and converged == False):
            print("Iteration did not converge at time step " + str(self.model_part.ProcessInfo[TIME]))
            if('stop_Newton_Raphson_if_not_converged' in self.Parameters):
                if(self.Parameters['stop_Newton_Raphson_if_not_converged'] == True):
                    sys.exit("Sorry, my boss does not allow me to continue. The time step did not converge at time step " + str(self.model_part.ProcessInfo[TIME]) + ", it = " + str(it) + ", max_iter = " + str(self.max_iter))
                else:
                    print('However, the iteration will still be proceeded' + ", it = " + str(it) + ", max_iter = " + str(self.max_iter))
            else:
                sys.exit("Sorry, my boss does not allow me to continue. The time step did not converge at time step " + str(self.model_part.ProcessInfo[TIME]) + ", it = " + str(it) + ", max_iter = " + str(self.max_iter))
        print("modified_thomas_strategy.PerformNewtonRaphsonIteration converged after " + str(it) + " steps")

    #######################################################################
    def Predict(self):
        self.time_scheme.Predict(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

        if self.model_part.NumberOfMasterSlaveConstraints() > 0:

            for constraint in self.model_part.MasterSlaveConstraints:
                constraint.ResetSlaveDofs(self.model_part.ProcessInfo)
                constraint.Apply(self.model_part.ProcessInfo)

            self.space_utils.SetToZeroVector(self.Dx)
            self.time_scheme.Update(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

    #######################################################################
    def InitializeSolutionStep(self):
        if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
            #initialize the list of degrees of freedom to be used
            self.builder_and_solver.SetUpDofSet(self.time_scheme,self.model_part)
            #reorder the list of degrees of freedom to identify fixity and system size
            self.builder_and_solver.SetUpSystem(self.model_part)
            #reorder the system dof id
            self.system_reorderer.Execute()
            #allocate memory for the system and preallocate the structure of the matrix
            self.builder_and_solver.ResizeAndInitializeVectors(self.pA,self.pDx,self.pb,self.model_part)
            #updating references
            self.A = (self.pA).GetReference()
            self.Dx = (self.pDx).GetReference()
            self.b = (self.pb).GetReference()
        if(self.SolutionStepIsInitialized == False):
            self.builder_and_solver.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)
            self.time_scheme.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        for proc in self.attached_processes:
            proc.ExecuteInitializeSolutionStep()

    #######################################################################
    def ExecuteIteration(self,echo_level,MoveMeshFlag,CalculateNormDxFlag,it,alpha):
        #reset system matrices and vectors prior to rebuild
        if (it == 0):
            self.space_utils.SetToZeroMatrix(self.A)
        self.space_utils.SetToZeroVector(self.Dx)
        self.space_utils.SetToZeroVector(self.b)

        self.time_scheme.InitializeNonLinIteration(self.model_part,self.A,self.Dx,self.b)
        print("modified_thomas_strategy.ExecuteIteration:InitializeNonLinIteration is called")

        #build and solve the problem
        if(self.Parameters['decouple_build_and_solve'] == False):
            if it == 0:
                self.linear_solver.SetAdditionalPhysicalData(True)
                self.builder_and_solver.BuildAndSolve(self.time_scheme, self.model_part, self.A, self.Dx, self.b)
            else:
                self.linear_solver.SetAdditionalPhysicalData(False)
                self.builder_and_solver.BuildRHSAndSolve(self.time_scheme, self.model_part, self.A, self.Dx, self.b)
        else:
            if it == 0:
                self.builder_and_solver.Build(self.time_scheme,self.model_part,self.A,self.b)
                self.builder_and_solver.ApplyDirichletConditions(self.time_scheme,self.model_part,self.A,self.Dx,self.b)
                self.dof_util.ListDofs(self.builder_and_solver.GetDofSet(),self.builder_and_solver.GetEquationSystemSize())
                #provide data for the preconditioner and linear solver
                if self.linear_solver.AdditionalPhysicalDataIsNeeded():
                    self.linear_solver.ProvideAdditionalData(self.A,self.Dx,self.b,self.builder_and_solver.GetDofSet(),self.model_part)
                self.linear_solver.Solve(self.A,self.Dx,self.b)
            else:
                self.builder_and_solver.BuildRHS(self.time_scheme,self.model_part,self.b)
                # print("normb at BuildRHS: " + str(self.space_utils.TwoNorm(self.b)))
                self.builder_and_solver.ApplyDirichletConditions(self.time_scheme,self.model_part,self.A,self.Dx,self.b)
                self.linear_solver.Solve(self.A,self.Dx,self.b)

        # print("At solve, A = " + str(self.A))
        # print("At solve, rhs = " + str(self.b))
        # print("At solve, dx = " + str(self.Dx))

        #calculate reaction process
        if not(self.calculate_reaction_process == None):
            self.calculate_reaction_process.Execute()

        if(self.Parameters['list_plastic_points'] == True):
            self.output_util.ListPlasticPoints(self.model_part)

#        diagAstr = ""
#        for i in range(0, self.A.Size1()):
#            diagAstr = diagAstr + ", " + str(self.A[(i, i)])
#        print("diagonal A:" + diagAstr)

        #full output if needed
        if( self.PrintSparsity ):
            self.PlotSparsityScheme( self.A )
#            wr = UblasMatrixIO()
#            wr.WriteHB(self.A, self.b, "matrix" + str(self.solveCounter) + "." + str(self.iterationCounter) + ".hb.dat")
#            self.space_utils.WriteMatrixMarketMatrix("matrix" + str(self.solveCounter) + "." + str(self.iterationCounter) + ".mm",self.A,False)
            #petsc_utils.DumpUblasCompressedMatrixVector("tempAb", self.A, self.b, False)

        if(echo_level >= 3):
            print("SystemMatrix = " + str(self.A))
        #printA = []
        #printdx = []
        #printb = []
        #for i in range(0,len(self.Dx)):
        #    if( abs(self.Dx[i]) < 1.0e-10 ):
         #       printdx.append(0.0)
         #   else:
         #       printdx.append(self.Dx[i])
         #   if( abs(self.b[i]) < 1.0e-6 ):
         #       printb.append(0.0)
         #   else:
         #       printb.append(self.b[i])
         #   row = []
         #   for j in range(0,len(self.Dx)):
         #       if( abs(self.A[(i,j)]) < 1.0 ):
         #           row.append( 0.0 )
         #       else:
         #           row.append(self.A[(i,j)])
         #   printA.append(row)
            print("solution obtained = " + str(self.Dx))
            #formatted_printdx = [ '%.6f' % elem for elem in printdx ]
            #print formatted_printdx
            #formatted_printb = [ '%.4f' % elem for elem in printb ]
            print("RHS = " + str(self.b))
        #print formatted_printb
        #print "Matrix: "
        #for i in range(0,len(self.Dx)):
        #    formatted_printA = [ '%.1f' % elem for elem in printA[i] ]
        #    print(formatted_printA)
        self.AnalyseSystemMatrix(self.A)

        #calculate the norm of the "correction" Dx
        if(CalculateNormDxFlag == True):
            normDx = self.space_utils.TwoNorm(self.Dx)
        else:
            normDx = 0.0

        return normDx

    #######################################################################
    def FinalizeNonLinIteration(self,ConvergedFlag,MoveMeshFlag):
        if not ConvergedFlag:
            #perform update
            self.time_scheme.Update(self.model_part,self.builder_and_solver.GetDofSet(),self.A,alpha*self.Dx,self.b)
            print("modified_thomas_strategy.FinalizeNonLinIteration:Update is called")

            #move the mesh as needed
            if(MoveMeshFlag == True):
                self.time_scheme.MoveMesh(self.model_part.Nodes)
#        print("b:" + str(self.b))
#        print("Dx:" + str(self.Dx))
#        print("A:" + str(self.A))

        self.time_scheme.FinalizeNonLinIteration(self.model_part,self.A,self.Dx,self.b)
        print("modified_thomas_strategy.FinalizeNonLinIteration:FinalizeNonLinIteration is called")

    #######################################################################
    def FinalizeSolutionStep(self,CalculateReactionsFlag):
        if(CalculateReactionsFlag == True):
            self.model_part.ProcessInfo[SET_CALCULATE_REACTION] = True
            self.builder_and_solver.CalculateReactions(self.time_scheme,self.model_part,self.A,self.Dx,self.b)
            self.model_part.ProcessInfo[SET_CALCULATE_REACTION] = False

        #Finalisation of the solution step
        self.time_scheme.FinalizeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        self.builder_and_solver.FinalizeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        self.time_scheme.Clean()
        self.linear_solver.Clear()
        #reset flags for the next step
        self.SolutionStepIsInitialized = False

        for proc in self.attached_processes:
            proc.ExecuteFinalizeSolutionStep()
        print("modified_thomas_strategy.FinalizeSolutionStep is called")

    #######################################################################
    def Clear(self):
        self.space_utils.ClearMatrix(self.pA)
        self.space_utils.ResizeMatrix(self.A,0,0)

        self.space_utils.ClearVector(self.pDx)
        self.space_utils.ResizeVector(self.Dx,0)

        self.space_utils.ClearVector(self.pb)
        self.space_utils.ResizeVector(self.b,0)

        #updating references
        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()

        self.builder_and_solver.SetDofSetIsInitializedFlag(False)
        self.builder_and_solver.Clear()
        self.time_scheme.Clear()

        for proc in self.attached_processes:
            proc.ExecuteFinalize()

    #######################################################################
    def SetEchoLevel(self,level):
        self.echo_level = level
        self.builder_and_solver.SetEchoLevel(level)

    #######################################################################
    def AnalyseSystemMatrix(self,  A):
        max = 0.0
        for i in range(0,  A.Size1()):
           if( abs(A[(i, i)]) > max ):
               max = A[(i, i)]

#        nonzeros = 0
#        for i in range(0,  A.Size1()):
#            for j in range(0,  A.Size2()):
#                if( abs(A[(i, j)]) > 1e-16 ):
#                    nonzeros = nonzeros + 1

        print("#############################")
        print("Number of rows: " +str(A.Size1()) )
        print("Number of columns: " +str(A.Size2()) )
#        print("Number of entries: " +str(nonzeros) )
        print("Max in Diagonal: " +str(max) )
        print("#############################")

    #######################################################################
    def PlotSparsityScheme(self, A):
        try:
            import Gnuplot, Gnuplot.PlotItems, Gnuplot.funcutils
        except ImportError:
            # kludge in case Gnuplot hasn't been installed as a module yet:
            import __init__
            Gnuplot = __init__
            import PlotItems
            Gnuplot.PlotItems = PlotItems
            import funcutils
            Gnuplot.funcutils = funcutils
        print("gnuplot-python imported")
        g = Gnuplot.Gnuplot(debug=1)
        g.clear()
        #g.plot(Gnuplot.Func('sin(x)'))
        #self.wait('hit Return to continue')
        file = open("matrix.dat",'w')
        for i in range(0, A.Size1()):
            for j in range(0, A.Size2()):
                tmp = A[(i,j)]
                if( (tmp > 1.0e-9) or (tmp < -1.0e-9) ):
                   #file.write( str(tmp) +"\t" )
                   file.write( "1.0 " )
                else:
                   file.write("0.0 ")
            file.write("\n")
        file.close()
        g("set term postscript")
        g("set size square")
        g("set output 'matrix.ps'")
        g("set zrange [0.5:1.5]")
        g("set pm3d map")
        g("splot 'matrix.dat' matrix with dots")


    def wait(self,str=None, prompt='Press return to show results...\n'):
        if str is not None:
            print(str)
        raw_input(prompt)


