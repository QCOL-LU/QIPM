import numpy as np
from copy import deepcopy



#===========================================================================
# Linear optimization problem with a known interior solution 
#===========================================================================
def generate_lo_problem_with_int(m, n, parameters):

	np.random.seed(parameters.seed)	

	int_x 		= np.random.rand(n)
	int_s 		= np.random.rand(n)
	int_y 		= np.random.rand(m) - 0.5

	int_x 		= (parameters.norm_x/np.linalg.norm(int_x))*int_x
	int_s 		= (parameters.norm_s/np.linalg.norm(int_s))*int_s
	int_y 		= (parameters.norm_y/np.linalg.norm(int_y))*int_y

	A 			= np.random.rand(m, n) - .5 
	A 			= np.matmul(A, A.T) if parameters.make_psd == True else A

	u, s, v 	= np.linalg.svd(A, full_matrices=False)
	s 			= np.linspace(parameters.norm_A, parameters.norm_A/parameters.condition_number, min(m, n))

	A 			= np.dot(u * s, v) if parameters.symmetry == False else np.dot(u * s, u.T)

	b 			= np.matmul(A, int_x)
	c 			= np.matmul(A.T, int_y) + int_s

	if parameters.norm_b != -1:
		temp_norm_b = np.linalg.norm(b)
		coef 		= parameters.norm_b / temp_norm_b
		int_x 		= coef * int_x
		b 			= coef * b

	results 	= (A, b, c, int_x, int_y, int_s)

		
	return results




