from ipm.print_methods.print_iterative_refinement import *
from copy import deepcopy


from time import time
import numpy as np
import sys


#===========================================================================
# Iterative Refinement Algorithm
#===========================================================================
def iterative_refinement(LO):

	#-----------------------------------------------------------------------
	# Initialize the paramters of the algorithm
	#-----------------------------------------------------------------------
	nabla 			= LO.Params.ScalFact					# Scaling factor (1e2)
	rho 			= LO.Params.IncScalLim  				# Incremental scaling limit (10)
	precision 		= LO.Params.IR_Precision				# Desired precision (1e-10)
	
	iteration 		= 0
	start_time 		= time()

	old_r 			= 1e8
	old_y 			= LO.y
	old_x 			= LO.x
	old_s 			= LO.c - np.dot(LO.A.T, LO.y)


	if LO.Params.IR_Verbosity > -1:
		print("\nThe Iterative Refinement algorithm starts running ...\n")


	#-----------------------------------------------------------------------
	# Solve the modified LO problem approximately
	#----------------------------------------------------------------------- 
	LO  			= LO.linear_optimizer(LO)
	LO_residual 	= deepcopy(LO)
	

	#-----------------------------------------------------------------------
	# Infinite Loop
	#-----------------------------------------------------------------------
	while (True):

		#-------------------------------------------------------------------
		# Calculate b_bar and c_bar 
		#-------------------------------------------------------------------
		b_hat 		= LO.b - np.dot(LO.A,  LO.x)
		c_hat 		= LO.c - np.dot(LO.A.T, LO.y) 
		
		LO.s		= c_hat


		#-------------------------------------------------------------------
		# Calculate the max of residual (r)
		#-------------------------------------------------------------------
		# print(min(np.multiply(c_hat, LO.x)) )
		# r			= max(max(abs(b_hat)), max(-c_hat), sum(abs(np.multiply(c_hat, LO.x))) )
		r  			= max(np.linalg.norm(b_hat), max(-c_hat), sum(abs(np.multiply(c_hat, LO.x)))) 
		

		if r > old_r: 
			LO.y 	= old_y
			LO.x 	= old_x
			LO.s 	= old_s
			print("The iterative refinement algorithm stopped because of inexact solution returned by IPM.")
			break

		old_x 		= LO.x 
		old_y 		= LO.y
		old_s		= LO.s

		old_r 		= r

		#-------------------------------------------------------------------
		# Calculate the step_size using r_bar and nabla
		#-------------------------------------------------------------------
		r_bar 		= max(r, 1/(nabla * rho))
		nabla 		= 2**(np.ceil(- np.log2(r_bar)))


		#-------------------------------------------------------------------
		# Calculate the run-time
		#-------------------------------------------------------------------
		end_time 	= time()
		run_time 	= end_time - start_time


		#-------------------------------------------------------------------
		# Print the summary of an iteration
		#-------------------------------------------------------------------
		print_iteration_iterative_refinement(LO, r, r_bar, nabla, iteration, run_time)


		if (1/nabla < LO.Params.Stop_Precision):
			print()
			print("The solution quality is limited by the precision of the linear system solver.")
			break

		#-------------------------------------------------------------------
		# Stopping condition
		#-------------------------------------------------------------------
		if (r <= precision): break


		#-------------------------------------------------------------------
		# Calculate the b_bar and c_bar
		#-------------------------------------------------------------------
		b_bar 		= nabla * b_hat
		c_bar 		= nabla * c_hat
		l_bar 		= - nabla * LO.x

		#-------------------------------------------------------------------
		# Solve the residual LO problem approximately
		#-------------------------------------------------------------------
		LO_residual.b 	= b_bar - np.dot(LO.A, l_bar)
		LO_residual.c 	= c_bar
		
		LO_residual 	= LO_residual.linear_optimizer(LO_residual)	

		#-------------------------------------------------------------------
		# Update the solution 
		#-------------------------------------------------------------------
		LO.x 		= LO.x + (LO_residual.x + l_bar) * (1 / nabla)
		LO.y 		= LO.y + LO_residual.y * (1 / nabla)


		iteration 	+= 1



	#-----------------------------------------------------------------------
	# Update the iterates
	#-----------------------------------------------------------------------
	print_final_iterative_refinement(LO, iteration, run_time)

	
	return LO



