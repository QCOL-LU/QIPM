import numpy as np
import sys
import time






#===========================================================================
# Eigenvalues using Quantum Phase Estimation
#===========================================================================
class Output:
	def __init__(self, cofficent_matrix, rhs_vector):
		self.cofficent_matrix 		= cofficent_matrix
		self.rhs_vector 			= rhs_vector

	def __str__(self): 
		np.set_printoptions(precision=3, formatter={'float': '{:+.3f}'.format})
		result 		= "\n"
		result 		+= (12*"-" + "Results of the linear system solver" + 12*"-") + "\n"
		if self.Is_Quantum == True:
			result 	+= "{:<25}{:}".format("Quantum solution:", self.solution) + "\n"
		result 		+= "{:<25}{:}".format("Exact solution:", self.exact_solution) + 2 * "\n"
		result 		+= "{:<25}{:.2e}".format("Norm of difference:", self.norm_of_difference) +"\n"
		result 		+= "{:<25}{:.2e}".format("Norm of residual:", self.norm_of_residual) + 2*"\n"
		if self.Is_Quantum == True:
			result 	+= "{:<25}{:.2e}".format("QLSA relative error:", self.LS_Precision)  + "\n" 
			result 	+= "{:<25}{:}".format("Is sign changed:", self.is_sign_changed)  + "\n"

		if self.LS_Method == "IR-LS":
			result 	+= "{:<25}{:}".format("iteration:", self.iteration)  + "\n"

		result 		+= "\n" + "{:<25}{:.2e}".format("Time (s):", self.end_time - self.start_time) + "\n"
		result 		+= (61*"=")

		

		return result

	


