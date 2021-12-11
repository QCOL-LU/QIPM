from ipm.general_methods.qlsa import QLSA
from ipm.print_methods.print_IPM import *

import numpy as np
from time import time
import sys

#===========================================================================
# Calculate the ratio such that s and x remains non-negative
#===========================================================================
def ratio(x_vec, delta_x_vec):
	ratio 	= 1
	for (x, delta_x) in zip(x_vec, delta_x_vec):
		if delta_x < 0:
			ratio 	= min(- x / delta_x , ratio)
	return ratio



#===========================================================================
# Calculate alpha-star using backtracking
#===========================================================================
def calculate_alpha_star(LO, alpha_star, delta_x, delta_y, delta_s):

	beta_1 			= LO.Params.Beta_1					# Default value: 0.1
	beta_2 			= LO.Params.Beta_2 					# Default value: 1 - 5e-4

	omega 			= (LO.Params.Omega)					# Default value: 1e8
	gamma 			= LO.Params.Gamma					# Default value: 0.5

	alpha_hat_dec 	= LO.Params.AlphaHatDec				# Default value: 1 - 1e-3
	precision 		= LO.Params.LO_Precision

	is_neighbor 	= False
	alpha_hat 		= alpha_star
	compl 			= LO.complementarity()

	while not is_neighbor: 

		x_temp 			= LO.x + alpha_hat * delta_x
		y_temp 			= LO.y + alpha_hat * delta_y	
		s_temp 			= LO.s + alpha_hat * delta_s

		compl_temp 		= np.dot(x_temp, s_temp)

		is_neighbor		= True



		for (xi, si) in zip(x_temp, s_temp):
			if (xi*si < gamma * compl_temp / LO.n):
				alpha_hat 		*= alpha_hat_dec
				is_neighbor 	= False
				break

		if not is_neighbor: 
			continue


		epsilon_primal 	= np.linalg.norm(np.dot(LO.A, x_temp) - LO.b)

		if (epsilon_primal > max(compl_temp/gamma, precision)):
			alpha_hat 		*= alpha_hat_dec
			is_neighbor 	= False
			continue


		epsilon_dual 	= np.linalg.norm(np.dot(LO.A.T, y_temp) + s_temp - LO.c)

		if (epsilon_dual > max(compl_temp/gamma, precision)):
			alpha_hat 		*= alpha_hat_dec
			is_neighbor 	= False
			continue


		if (compl_temp > (1 - alpha_hat * (1 - beta_2) )* compl) :
			alpha_hat 		*= alpha_hat_dec
			is_neighbor 	= False
			continue

	#-------------------------------------------------------------------
	# Update iterates
	#-------------------------------------------------------------------
	LO.x 		= x_temp
	LO.y 		= y_temp
	LO.s 		= s_temp

	return (LO, alpha_hat)


#===========================================================================
# Inexact Infeasible Interior Point Method
#===========================================================================
def inexact_infeasible_IPM(LO):
	#-----------------------------------------------------------------------
	# Initialize the paramters of the algorithm
	#-----------------------------------------------------------------------
	beta_1 			= LO.Params.Beta_1					# Default value: 0.1
	omega 			= LO.Params.Omega					# Default value: 1e8
	precision 		= LO.Params.LO_Precision

	iteration 		= 0
	
	start_time 		= time()
	

	if LO.Params.LO_Verbosity > -1:
		if LO.Params.Is_Quantum == False:
			print("\nThe Inexact Infeasible Interior-Point Method starts running ...\n")
		else:
			print("\nThe Quantum Inexact Infeasible Interior-Point Method starts running ...\n")



	#-----------------------------------------------------------------------
	# Initialize the paramters of the problem
	#-----------------------------------------------------------------------
	LO.x 		= np.ones(LO.n) * omega  				# primal solution
	LO.y 		= np.zeros(LO.m)						# dual solution
	LO.s 		= np.ones(LO.n) * omega 	 			# dual slack
	


	#-----------------------------------------------------------------------
	# Loop with 
	#-----------------------------------------------------------------------
	while (True):

		#-------------------------------------------------------------------
		# Calculate complementarity 
		#-------------------------------------------------------------------
		mu 		= LO.complementarity() * beta_1 / LO.n


		#-------------------------------------------------------------------
		# Calculate the intermediate terms
		#-------------------------------------------------------------------
		X		= np.diag(LO.x)
		S 		= np.diag(LO.s)

		c_ATy 	= LO.c - np.dot(LO.A.T, LO.y) 							# c - A^T.y
		XS_1 	= np.dot(X, np.linalg.inv(S))							# X.S^(-1)
		SX_1 	= np.dot(S, np.linalg.inv(X))							# S.X^(-1)
		s_1 	= np.reciprocal(LO.s)									# S^(-1).e
		x_1 	= np.reciprocal(LO.x)									# X^(-1).e

		M 		= LO.A.dot(XS_1).dot(LO.A.T) 							# A.X.S^(-1).A^T
		
		r 		= LO.b - np.dot(LO.A, LO.x) - mu * np.dot(LO.A, s_1) + \
					LO.A.dot(XS_1).dot(c_ATy) 							# b - A.x - mu.As^(-1) + A.X.S^(-1).(c - A^T.y)

		

		#-------------------------------------------------------------------
		# Liner system solver solution
		#-------------------------------------------------------------------
		if LO.Params.Is_Quantum == True:
			delta_y, norm_of_residual, is_sign_changed 	= QLSA(M, r, LO.Params) 

		else:
			is_sign_changed 	= False
			norm_of_residual 	= 0
			delta_y 			= np.linalg.solve(M, r)		


		#-------------------------------------------------------------------
		# Calculate delta_x and delta_s
		#-------------------------------------------------------------------
		delta_s 		= c_ATy - LO.s - np.dot(LO.A.T, delta_y)
		delta_x 		= mu * s_1 - LO.x - np.dot(XS_1, delta_s)

		 
		#-------------------------------------------------------------------
		# Calculate the step-size for updating the iterates
		#-------------------------------------------------------------------
		alpha_star_x 	= ratio(LO.x, delta_x)
		alpha_star_s 	= ratio(LO.s, delta_s)
		alpha_star 		= min(alpha_star_x, alpha_star_s, 1.0)

		
		#-------------------------------------------------------------------
		# Backtracking to calculate alpha_hat
		#-------------------------------------------------------------------
		LO, alpha_hat 	= calculate_alpha_star(LO, alpha_star, delta_x, delta_y, delta_s)
		
		#-------------------------------------------------------------------
		# Calculate the run-time
		#-------------------------------------------------------------------
		end_time 	= time()
		run_time 	= end_time - start_time


		#-------------------------------------------------------------------
		# Check infeasibility of the problem
		#-------------------------------------------------------------------
		if (max(abs(entry) for entry in np.concatenate((LO.x, LO.s)) ) > 2 * LO.n * omega):
			print("The problem is infeasible.")
			break
			


		#-------------------------------------------------------------------
		# Print the results
		#-------------------------------------------------------------------
		print_iteration_IPM(LO, M, r, iteration , norm_of_residual, is_sign_changed, alpha_star, alpha_hat, run_time)


		iteration 	+= 1


		#-------------------------------------------------------------------
		# Stopping condition
		#-------------------------------------------------------------------
		if (LO.complementarity() <= precision  ): break

		if (alpha_hat < LO.Params.Stop_Precision):
			print()
			print("The solution quality is limited by the precision of the linear system solver.")
			break


		condition_number 	= np.linalg.cond(M)

		if (condition_number > LO.Params.Stop_Cond_Num):
			print()
			print("The solution quality is limited by the condition number that the linear system solver can handel.")
			break
	
	print_final_IPM(LO, iteration, run_time)

	return LO



