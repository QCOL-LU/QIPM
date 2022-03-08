import numpy as np
import sys


#===========================================================================
# Linear optimization problem class 
#===========================================================================
class Parameters():

	#-----------------------------------------------------------------------
	# Initialize the problem
	#-----------------------------------------------------------------------
	def __init__(self):

		#-------------------------------------------------------------------
		# Default values of problem generator
		#-------------------------------------------------------------------
		self.Problem_Type 		= "LO"
		self.norm_x 			= 5
		self.norm_s 			= 5
		self.norm_y 			= 10 
		self.norm_b 			= -1
		self.norm_c 			= -1
		self.norm_A				= 1 

		self.decimals			= 3 
		self.condition_number	= 2

		self.seed 				= 6261997 

		self.has_interior		= True
		self.has_optimal		= True

		self.make_psd			= False 
		self.do_print			= False 
		self.qlsa_print			= False 
		self.symmetry			= False

		self.time_limit 		= 100


		#-------------------------------------------------------------------
		# Default values of paramters
		#-------------------------------------------------------------------
		self.Method 			= "II-IPM"
		self.LS_Method 			= "LS"
		self.Is_Quantum 		= False
		self.Is_Simulator 		= True
		self.Is_Noisy			= False

		self.LS_Precision 		= 1e-1
		self.IR_LS_Precision 	= 1e-8

		self.IR_Precision 		= 1e-10
		self.LO_Precision 		= 1e-8
		self.Stop_Precision 	= 1e-16					# If the step-length becomes less than the algorithm stops
		self.Stop_Cond_Num		= 1e3

		self.IR_Verbosity 		= 2
		self.LO_Verbosity 		= 1

		self.num_ancillae 		= 3
		self.num_time_slices	= 1
		self.expansion_order	= 2
		self.HHL_Method 		= 1
		self.qlsa_precision  	= 1e0




		#-------------------------------------------------------------------
		# Default values of inexact_infeasible_IPM paramters
		#-------------------------------------------------------------------
		self.Beta_1 			= 0.1
		self.Beta_2 			= 1 - 5e-4

		self.Omega 				= 1e8
		self.Gamma 				= 0.5

		self.AlphaHatDec 		= 1 - 1e-3


		#-------------------------------------------------------------------
		# Default values of Iterative Refinement algorithm paramters
		#------------------------------------------------------------------- 
		self.ScalFact 			= 1e2							# Scaling factor
		self.IncScalLim 		= 10  							# Incremental scaling limit


		self.LS_ScalFact 		= 1								# Scaling factor
		self.LS_IncScalLim 		= 2  							# Incremental scaling limit

