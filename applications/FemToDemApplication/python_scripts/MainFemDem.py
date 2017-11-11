from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import system python modules
import time as timer
import os

# Import kratos core and applications
import KratosMultiphysics.SolidMechanicsApplication	 as KratosSolid
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers
import KratosMultiphysics.FemToDemApplication   as KratosFemDem
import MainSolidFEM
import main_script as MainDEM

import adaptive_mesh_refinement_utility
import gid_output_utility
import cleaning_utility

def Wait():
	input("Press Something")


class FEM_Solution(MainSolidFEM.Solution):

	def Info(self):
		print("FEM part of the FEMDEM application")

#============================================================================================================================					
	def __init__(self):

		#### TIME MONITORING START ####
		# Time control starts		
		print(timer.ctime())
		# Measure process time
		self.t0p = timer.clock()
		# Measure wall time
		self.t0w = timer.time()
		#### TIME MONITORING END ####


		#### PARSING THE PARAMETERS ####

		# Import input
		parameter_file = open("ProjectParameters.json",'r')
		self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

		# set echo level
		self.echo_level = self.ProjectParameters["problem_data"]["echo_level"].GetInt()

		print(" ")

		# defining the number of threads:
		num_threads =  self.GetParallelSize()
		print("::[KSM Simulation]:: [OMP USING",num_threads,"THREADS ]")
		#parallel.PrintOMPInfo()


		print(" ")
		print("::[KSM Simulation]:: [Time Step:", self.ProjectParameters["problem_data"]["time_step"].GetDouble()," echo:", self.echo_level,"]")

		#### Model_part settings start ####

		# Defining the model_part
		self.main_model_part = KratosMultiphysics.ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())

		if (self.ProjectParameters["solver_settings"]["solution_type"].GetString() == "Dynamic"):
			self.main_model_part.ProcessInfo.SetValue(KratosFemDem.IS_DYNAMIC, 1)
		else:
			self.main_model_part.ProcessInfo.SetValue(KratosFemDem.IS_DYNAMIC, 0)

		self.time_step  = self.ProjectParameters["problem_data"]["time_step" ].GetDouble()
		self.start_time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
		self.end_time   = self.ProjectParameters["problem_data"]["end_time"  ].GetDouble()


		#self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DIMENSION, self.ProjectParameters["problem_data"]["dimension"].GetInt())
		self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())
		self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.time_step)
		self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.start_time)


		###TODO replace this "model" for real one once available in kratos core
		self.Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

		#construct the solver (main setting methods are located in the solver_module)
		solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
		self.solver   = solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])

		#### Output settings start ####
		self.problem_path = os.getcwd()
		self.problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()

		
#============================================================================================================================
	def AddMaterials(self):

		# Assign material to model_parts (if Materials.json exists)
		import process_factory

		if os.path.isfile("Materials.json"):
			materials_file = open("Materials.json",'r')
			MaterialParameters = KratosMultiphysics.Parameters(materials_file.read())

			if(MaterialParameters.Has("material_models_list")):

				## Get the list of the model_part's in the object Model
				for i in range(self.ProjectParameters["solver_settings"]["problem_domain_sub_model_part_list"].size()):
					part_name = self.ProjectParameters["solver_settings"]["problem_domain_sub_model_part_list"][i].GetString()
					if( self.main_model_part.HasSubModelPart(part_name) ):
						self.Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})


				assign_materials_processes = process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses( MaterialParameters["material_models_list"] )

			for process in assign_materials_processes:
				process.Execute()
		else:
			print(" No Materials.json found ")
				
#============================================================================================================================		  
	def AddProcesses(self):

		# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
		## Get the list of the submodel part in the object Model
		for i in range(self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
			part_name = self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
			if( self.main_model_part.HasSubModelPart(part_name) ):
				self.Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})
		
		# Obtain the list of the processes to be applied
		import process_handler

		process_parameters = KratosMultiphysics.Parameters("{}") 
		process_parameters.AddValue("echo_level", self.ProjectParameters["problem_data"]["echo_level"])
		process_parameters.AddValue("constraints_process_list", self.ProjectParameters["constraints_process_list"])
		process_parameters.AddValue("loads_process_list", self.ProjectParameters["loads_process_list"])
		if( self.ProjectParameters.Has("problem_process_list") ):
			process_parameters.AddValue("problem_process_list", self.ProjectParameters["problem_process_list"])
		if( self.ProjectParameters.Has("output_process_list") ):
			process_parameters.AddValue("output_process_list", self.ProjectParameters["output_process_list"])

		return (process_handler.ProcessHandler(self.Model, process_parameters))

#============================================================================================================================	
	def Run(self):

		self.Initialize()

		self.RunMainTemporalLoop()

		self.Finalize()
		
#============================================================================================================================		
	def Initialize(self):


		#### INITIALIZE ####
		
		# Add variables (always before importing the model part)
		self.solver.AddVariables()
		
		# Read model_part (note: the buffer_size is set here) (restart is read here)
		self.solver.ImportModelPart()

		# Add dofs (always after importing the model part)
		if((self.main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
			if(self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
				self.solver.AddDofs()
		else:
			self.solver.AddDofs()


		# Add materials (assign material to model_parts if Materials.json exists)
		self.AddMaterials()
		

		# Add processes
		self.model_processes = self.AddProcesses()
		self.model_processes.ExecuteInitialize()


		# Print model_part and properties
		if(self.echo_level > 1):
			print("")
			print(self.main_model_part)
			for properties in self.main_model_part.Properties:
				print(properties)


		#### START SOLUTION ####

		self.computing_model_part = self.solver.GetComputingModelPart()

		## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
		self.solver.Initialize()
		#self.solver.InitializeStrategy()
		self.solver.SetEchoLevel(self.echo_level)

		
		# Initialize GiD  I/O (gid outputs, file_lists)
		self.SetGraphicalOutput()
		
		self.GraphicalOutputExecuteInitialize()


		print(" ")
		print("::[KSM Simulation]:: Analysis -START- ")

		self.model_processes.ExecuteBeforeSolutionLoop()

		self.GraphicalOutputExecuteBeforeSolutionLoop()		

		# Set time settings
		self.step	   = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
		self.time	   = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

		self.end_time   = self.ProjectParameters["problem_data"]["end_time"].GetDouble()
		self.delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()


		###   ------  Initializing Adaptive Mesh Refinement  ----------####
		self.cleaning_util = cleaning_utility.CleaningUtility(self.problem_path)

		self.gid_output_util = gid_output_utility.GidOutputUtility(self.ProjectParameters,
			                                                       self.problem_name,
			                                                       self.start_time,
			                                                       self.end_time,
			                                                       self.delta_time)

		self.constitutive_law_utility = []  # must be changed->provisional TODO
		self.conditions_util = []



		self.activate_AMR = self.ProjectParameters["AMR_data"]["activate_AMR"].GetBool()
		self.current_id = 1


		# Initialize the AMR_util
		if(self.activate_AMR):
			self.AMR_util = adaptive_mesh_refinement_utility.AdaptiveMeshRefinementUtility(self.ProjectParameters,
				                                                                           self.start_time,
				                                                                           self.solver,
				                                                                           self.constitutive_law_utility,
				                                                                           gid_output_utility,
				                                                                           self.conditions_util,
				                                                                           self.problem_path)
			self.activate_AMR = self.AMR_util.Initialize() # check the amr criteria
			#print("point alex 223333", self.activate_AMR)
			#Wait()



#============================================================================================================================
	def RunMainTemporalLoop(self):
		
		# Solving the problem (time integration)
		while(self.time < self.end_time):
			
			self.InitializeSolutionStep()
			self.SolveSolutionStep()
			self.FinalizeSolutionStep()
	  


#============================================================================================================================
	def InitializeSolutionStep(self):

		neighbour_elemental_finder =  KratosMultiphysics.FindElementalNeighboursProcess(self.main_model_part, 2, 5)
		neighbour_elemental_finder.Execute()

		self.delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

		self.time = self.time + self.delta_time
		self.step = self.step + 1

		self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.step
		self.main_model_part.CloneTimeStep(self.time) 


		print(" [STEP:",self.step," TIME:", self.time,"]")

		#print("priemro")
		#Wait()

		# processes to be executed at the begining of the solution step
		self.model_processes.ExecuteInitializeSolutionStep()
		#print("segundo")
		#Wait()
		self.GraphicalOutputExecuteInitializeSolutionStep()
		#print("tercero")
		#Wait()
		self.solver.InitializeSolutionStep()
		#print("cuarto")
		#Wait()

#============================================================================================================================
	def SolveSolutionStep(self):

		self.clock_time = self.StartTimeMeasuring();
		#print("quinto")
		#Wait()
		self.solver.Solve()
		#print("sexto")
		#Wait()

		# ************* PRINTS DE PRUEBA ******************************************* 
		#print("*************************************")
		#print(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
		#new_list = []
		#for node in self.main_model_part.Nodes:
		#	pass
			#if node.Id == 3:
				#print("Id: ",node.Id)
			#	print("X: ",node.X)
				#print("Y: ",node.Y)
			#	print("Displ: ", node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT))


		#damaged_elements = []
		#for elem in self.main_model_part.Elements:
		#	damage = elem.GetValuesOnIntegrationPoints(KratosFemDem.DAMAGE_ELEMENT, self.main_model_part.ProcessInfo)
		#	if damage[0][0] > 0.0:
		#		print("Id elemento dañado: ", elem.Id)
		#		Wait()
				#print(damaged_elements)

			#print("damage : ", damage)
			#print("\n")

		#print("*************************************")
		
		#print(activate_AMR)
		#print("fjfjfjjfjfj  ",self.activate_AMR)
		#Wait()
		# Adaptive Mesh Refinement
		#print("solve sol step",self.activate_AMR )
		#Wait()


		#print("despues de solve sol step  / execute")
		#Wait()

		self.StopTimeMeasuring(self.clock_time,"Solving", False);

#============================================================================================================================
	def FinalizeSolutionStep(self):


		self.GraphicalOutputExecuteFinalizeSolutionStep()			

		# processes to be executed at the end of the solution step
		self.model_processes.ExecuteFinalizeSolutionStep()

		# processes to be executed before witting the output	  
		self.model_processes.ExecuteBeforeOutputStep()

		# write output results GiD: (frequency writing is controlled internally)
		self.GraphicalOutputPrintOutput()			

		# processes to be executed after witting the output
		self.model_processes.ExecuteAfterOutputStep()

		# Eliminates elements from the mesh with damage > 0.98
		#self.main_model_part.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)

		# let's check the damage threshold
		#if self.step > 14:
		#	for elem in self.main_model_part.Elements:
		#		if elem.Id	== 10:			
		#			thres = elem.GetValuesOnIntegrationPoints(KratosFemDem.STRESS_THRESHOLD, self.main_model_part.ProcessInfo)
		#			print("en python thres", thres)
		#			Wait()










		if(self.activate_AMR):
			self.refine, self.last_mesh = self.AMR_util.CheckAMR(self.time)
			if(self.refine):
				self.main_model_part, self.solver = self.AMR_util.Execute(self.main_model_part,
					                                         self.solver,
					                                         self.gid_output_util,
					                                         self.time,
					                                         self.current_id)

				# save print parameters
				printed_step_count = self.graphical_output.printed_step_count
				step_count = self.graphical_output.step_count

				#construct the new solver (main setting methods are located in the solver_module)
				#self.ProjectParameters["solver_settings"].RemoveValue("damp_factor_m")
				#self.ProjectParameters["solver_settings"].RemoveValue("dynamic_factor")

				#solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
				#self.solver   = solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])
				self.InitializeAfterAMR()

				# assign print parameters
				self.graphical_output.printed_step_count = printed_step_count
				self.graphical_output.step_count = step_count

				# test
				
				print("antes de imprimir")
				#Wait()
				self.GraphicalOutputPrintOutput()	


			
			elif(self.last_mesh):
				self.AMR_util.Finalize(self.main_model_part, self.current_id)





		#print("NODOS: ", cont) # es correcto el nuevo mdpa
		#print("computing model part",self.computing_model_part)
		#print(" *********model part*******",self.main_model_part)
		#print("tiempooo :", self.time)
		#print(self.graphical_output.printed_step_count)
		#print(self.graphical_output.step_count)
		#print("************************")
		#Wait()



#============================================================================================================================
	def Finalize(self):
		
		# Ending the problem (time integration finished)
		self.GraphicalOutputExecuteFinalize()		

		self.model_processes.ExecuteFinalize()

		print("::[KSM Simulation]:: Analysis -END- ")
		print(" ")

		# Check solving information for any problem
		#~ self.solver.InfoCheck() # InfoCheck not implemented yet.

		#### END SOLUTION ####

		# Measure process time
		tfp = timer.clock()
		# Measure wall time
		tfw = timer.time()

		print("::[KSM Simulation]:: [Elapsed Time = %.2f" % (tfw - self.t0w),"seconds] (%.2f" % (tfp - self.t0p),"seconds of cpu/s time)")

		print(timer.ctime())


 #============================================================================================================================	   
	def SetGraphicalOutput(self):
		from gid_output_process import GiDOutputProcess
		self.output_settings = self.ProjectParameters["output_configuration"]
		self.graphical_output = GiDOutputProcess(self.computing_model_part,
									  self.problem_name,
									  self.output_settings)		
	#============================================================================================================================
	def GraphicalOutputExecuteInitialize(self):
		self.graphical_output.ExecuteInitialize() 
	#============================================================================================================================		
	def GraphicalOutputExecuteBeforeSolutionLoop(self):
		# writing a initial state results file or single file (if no restart)
		if((self.main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
			if(self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
				self.graphical_output.ExecuteBeforeSolutionLoop()
	#============================================================================================================================			   
	def GraphicalOutputExecuteInitializeSolutionStep(self):
		self.graphical_output.ExecuteInitializeSolutionStep()
	#============================================================================================================================		
	def GraphicalOutputExecuteFinalizeSolutionStep(self):
		self.graphical_output.ExecuteFinalizeSolutionStep() 
	#============================================================================================================================		
	def GraphicalOutputPrintOutput(self):
		if(self.graphical_output.IsOutputStep()):
				self.graphical_output.PrintOutput()
	#============================================================================================================================
	def GraphicalOutputExecuteFinalize(self):
		self.graphical_output.ExecuteFinalize()
				

	#============================================================================================================================   
	def SetParallelSize(self, num_threads):
		parallel = KratosMultiphysics.OpenMPUtils()
		parallel.SetNumThreads(int(num_threads))
	#============================================================================================================================
	def GetParallelSize(self):
		parallel = KratosMultiphysics.OpenMPUtils()
		return parallel.GetNumThreads()	
	#============================================================================================================================	
	def StartTimeMeasuring(self):
		# Measure process time
		time_ip = timer.clock()
		return time_ip
	#============================================================================================================================
	def StopTimeMeasuring(self, time_ip, process, report):
		# Measure process time
		time_fp = timer.clock()
		if( report ):
			used_time = time_fp - time_ip
			print("::[KSM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

	#============================================================================================================================

	def InitializeAfterAMR(self):

		#### INITIALIZE ####
		"""
		# Add variables (always before importing the model part)
		self.solver.AddVariables()
		print("en main antes de import:", self.main_model_part.ProcessInfo[KratosMultiphysics.STEP])
		# Read model_part (note: the buffer_size is set here) (restart is read here)
		self.solver.ImportModelPart()
		print("en main despeus de import:", self.main_model_part.ProcessInfo[KratosMultiphysics.STEP])

		# Add dofs (always after importing the model part)
		if((self.main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
			if(self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
				self.solver.AddDofs()
		else:
			self.solver.AddDofs()

		"""
		# Add materials (assign material to model_parts if Materials.json exists)
		#self.AddMaterials()
		

		# Add processes
		self.model_processes = self.AddProcesses()
		self.model_processes.ExecuteInitialize()



		#### START SOLUTION ####

		self.computing_model_part = self.solver.GetComputingModelPart()

		## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
		self.solver.Initialize()
		#self.solver.InitializeStrategy()
		self.solver.SetEchoLevel(self.echo_level)

		
		# Initialize GiD  I/O (gid outputs, file_lists)
		self.SetGraphicalOutput()
		#print("stepppp2 :", self.main_model_part.ProcessInfo[KratosMultiphysics.STEP])
		self.GraphicalOutputExecuteInitialize()
		#print("stepppp3 :", self.main_model_part.ProcessInfo[KratosMultiphysics.STEP])


		#print(" ")
		#print("::[KSM Simulation]:: Analysis -START- ")

		#self.model_processes.ExecuteBeforeSolutionLoop()

		self.GraphicalOutputExecuteBeforeSolutionLoop()		

		# Set time settings
		#self.step	   = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
		#self.time	   = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

		#self.end_time   = self.ProjectParameters["problem_data"]["end_time"].GetDouble()
		#self.delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()


		###   ------  Initializing Adaptive Mesh Refinement  ----------####
		#self.cleaning_util = cleaning_utility.CleaningUtility(self.problem_path)

		#self.gid_output_util = gid_output_utility.GidOutputUtility(self.ProjectParameters,
			                                                       #self.problem_name,
			                                                       #self.start_time,
			                                                      # self.end_time,
			                                                       #self.delta_time)

		#self.constitutive_law_utility = []  # must be changed->provisional TODO
		#self.conditions_util = []



		#self.activate_AMR = self.ProjectParameters["AMR_data"]["activate_AMR"].GetBool()
		#self.current_id = 1


		# Initialize the AMR_util
		#if(self.activate_AMR):
			#self.AMR_util = adaptive_mesh_refinement_utility.AdaptiveMeshRefinementUtility(self.ProjectParameters,
				                                                                          # self.start_time,
				                                                                           #self.solver,
				                                                                          # self.constitutive_law_utility,
				                                                                           #gid_output_utility,
				                                                                           #self.conditions_util,
				                                                                           #self.problem_path)
			#self.activate_AMR = self.AMR_util.Initialize() # check the amr criteria
			#print("point alex 223333", self.activate_AMR)
			#Wait()



if __name__ == "__main__": 
	Solution().Run()


#============================================================================================================================
#============================================================================================================================
# Main de DEM application

class DEM_Solution(MainDEM.Solution):

    def Info(self):
        print("DEM part of the FEM-DEM application")



