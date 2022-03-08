from ipm.general_methods.ParametersDefault import Parameters
from ipm.ipm_methods._inexact_infeasible_IPM import inexact_infeasible_IPM
from ipm.ipm_methods._inexact_feasible_IPM import inexact_feasible_IPM
from ipm.general_methods._iterative_refinement import iterative_refinement

import numpy as np
import sys


#===========================================================================
# Linear optimization problem class 
#===========================================================================
class Model(Parameters):

	#-----------------------------------------------------------------------
	# Initialize the problem
	#-----------------------------------------------------------------------
	def __init__(self, A=None, b=None, c=None):


		#-------------------------------------------------------------------
		# Check the data type
		#-------------------------------------------------------------------
		if isinstance(A, np.ndarray) == True or isinstance(A, list) == True :
			self.A = np.array(A)
		else:
			sys.exit("A is not a 2D list or array.")


		if isinstance(b, np.ndarray) == True or isinstance(b, list) == True :
			self.b = np.array(b)
			self.m = len(b)
		else:
			sys.exit("b is not a list or array.")


		if isinstance(c, np.ndarray) == True or isinstance(c, list) == True :
			self.c = np.array(c)
			self.n = len(c)
		else:
			sys.exit("c is not a list or array.")


		#-------------------------------------------------------------------
		# Check the compatibility of dimensions 
		#-------------------------------------------------------------------
		if self.A.shape[0] != self.m:
			sys.exit("The dimension of matrix A and vector b does not match.")

		if self.A.shape[1] != self.n:
			sys.exit("The dimension of matrix A and vector c does not match.")



		#-------------------------------------------------------------------
		# Set the intial value for the solution of the problem 
		#-------------------------------------------------------------------	
		self.x 	= np.zeros(self.n)
		self.s 	= np.zeros(self.n)
		self.y	= np.zeros(self.m)


		self.Params = Parameters()


	#-----------------------------------------------------------------------
	# Import methods
	#-----------------------------------------------------------------------
	


	#-----------------------------------------------------------------------
	# Calculate the primal problem objective value (c^T.x)
	#-----------------------------------------------------------------------
	def primal_objective(self):
		return np.dot(self.c, self.x) 	


	#-----------------------------------------------------------------------
	# Calculate the dual problem objective value (b^T.y)
	#-----------------------------------------------------------------------
	def dual_objective(self):
		return np.dot(self.b, self.y) 


	#-----------------------------------------------------------------------
	# Calculate the norm of the primal problem's residual (||A.x - b||)
	#-----------------------------------------------------------------------
	def primal_residual(self):
		return np.linalg.norm(np.dot(self.A, self.x) - self.b)


	#-----------------------------------------------------------------------
	# Calculate the norm of the dual problem's residual (||c - A^T.y - s||)
	#-----------------------------------------------------------------------
	def dual_residual(self):
		return np.linalg.norm(self.c - np.dot(self.A.T, self.y) - self.s)


	#-----------------------------------------------------------------------
	# Calculate the complementarity (s^T.x)
	#-----------------------------------------------------------------------
	def complementarity(self):
		return np.dot(self.x, self.s)

	#-----------------------------------------------------------------------
	# Solve the LO problem
	#-----------------------------------------------------------------------
	def solve(self):

		self.print_paramters()
		
		if self.Params.Method == "II-IPM":
			self.Params.Is_Quantum 		= False
			solver 						= inexact_infeasible_IPM


		elif self.Params.Method == "II-QIPM":
			self.Params.Is_Quantum 		= True
			solver 						= inexact_infeasible_IPM

		elif self.Params.Method == "IF-IPM":
			self.Params.Is_Quantum 		= False
			solver 						= inexact_feasible_IPM 

		elif self.Params.Method == "IF-QIPM":
			self.Params.Is_Quantum 		= True
			solver 						= inexact_feasible_IPM


		elif "IR-" in self.Params.Method:
			solver 						= iterative_refinement

			if "II-QIPM" in self.Params.Method:
				self.Params.Is_Quantum 	= True
				self.linear_optimizer 	= inexact_infeasible_IPM

			elif "II-IPM" in self.Params.Method:
				self.Params.Is_Quantum 	= False
				self.linear_optimizer 	= inexact_infeasible_IPM

			elif "IF-QIPM" in self.Params.Method:
				self.Params.Is_Quantum 	= True
				self.linear_optimizer 	= inexact_feasible_IPM

			elif "IF-IPM" in self.Params.Method:
				self.Params.Is_Quantum 	= False
				self.linear_optimizer 	= inexact_feasible_IPM


			else:
				print("Please choose one of the following methods:")
				print("  1- \"II-QIPM\": Inexact Infeasible Quantum Interior-Point method")
				print("  2- \"II-IPM\": Inexact Infeasible Interior-Point method")
				print("  3- \"IF-QIPM\": Inexact Feasible Quantum Interior-Point method")
				print("  4- \"IF-IPM\": Inexact Feasible Interior-Point method")
				sys.exit("Please select the method correctly!")


		return solver(self)


	#-----------------------------------------------------------------------
	# Print updated paramters
	#-----------------------------------------------------------------------
	def print_paramters(self):
		default_paramters 	= Parameters()

		print("{:25}{}".format("Method:", self.Params.Method) )

		if "IR-" in self.Params.Method:

			if (default_paramters.IR_Precision != self.Params.IR_Precision):
				print("{:25}{:>8.2e}".format("IR_Precision:", self.Params.IR_Precision) )

			print("{:25}{:>8.2e}".format("IR_LS_Precision:", self.Params.IR_LS_Precision) )

			if (default_paramters.IR_Verbosity != self.Params.IR_Verbosity):
				print("{:25}{}".format("IR_Verbosity:", self.Params.IR_Verbosity) )

			if (default_paramters.ScalFact != self.Params.ScalFact):
				print("{:25}{:>8.2e}".format("ScalFact:", self.Params.ScalFact) )

			if (default_paramters.IncScalLim != self.Params.IncScalLim):
				print("{:25}{:>8.2e}".format("IncScalLim:", self.Params.IncScalLim) )



		if (default_paramters.LO_Precision != self.Params.LO_Precision):
			print("{:25}{:>8.2e}".format("LO_Precision:", self.Params.LO_Precision) )


		if (default_paramters.Stop_Precision != self.Params.Stop_Precision):
			print("{:25}{:>8.2e}".format("Stop_Precision:", self.Params.Stop_Precision) )
			

		if (default_paramters.LO_Verbosity != self.Params.LO_Verbosity):
			print("{:25}{}".format("LO_Verbosity:", self.Params.LO_Verbosity) )


		if (default_paramters.Beta_1 != self.Params.Beta_1):
			print("{:25}{:>8.2e}".format("Beta_1:", self.Params.Beta_1) )

		if (default_paramters.Beta_2 != self.Params.Beta_2):
			print("{:25}{:>8.2e}".format("Beta_2:", self.Params.Beta_2) )

		if (default_paramters.Omega != self.Params.Omega):
			print("{:25}{:>8.2e}".format("Omega:", self.Params.Omega) )

		if (default_paramters.Gamma != self.Params.Gamma):
			print("{:25}{:>8.2e}".format("Gamma:", self.Params.Gamma) )

		if (default_paramters.AlphaHatDec != self.Params.AlphaHatDec):
			print("{:25}{:>8.2e}".format("AlphaHatDec:", self.Params.AlphaHatDec) )
			


		if "Q" in self.Params.Method:

			if (default_paramters.Is_Simulator != self.Params.Is_Simulator):
				print("{:25}{}".format("Is_Simulator:", self.Params.Is_Simulator) )

			if (default_paramters.num_ancillae != self.Params.num_ancillae):
				print("{:25}{}".format("num_ancillae:", self.Params.num_ancillae) )

			if (default_paramters.num_time_slices != self.Params.num_time_slices):
				print("{:25}{}".format("num_time_slices:", self.Params.num_time_slices) )

			if (default_paramters.expansion_order != self.Params.expansion_order):
				print("{:25}{}".format("expansion_order:", self.Params.expansion_order) )




