import numpy as np
from copy import deepcopy



#===========================================================================
# Linear optimization problem with a known optimal solution 
#===========================================================================
def generate_lo_problem_with_opt(m, n, parameters):

	np.random.seed(parameters.seed)	

	mask        = [1 if ind < m else 0 for ind in range(n)]
    np.random.shuffle(mask)

    opt_x 		= np.random.rand(n)
    opt_x 		= np.multiply(opt_x,mask)
	opt_s 		= np.random.rand(n)
    opt_s 		= opt_s-np.multiply(opt_s,mask)
	opt_y 		= np.random.rand(m) - 0.5
	opt_x 		= (parameters.norm_x/np.linalg.norm(opt_x))*opt_x
	opt_s 		= (parameters.norm_s/np.linalg.norm(opt_s))*opt_s
	opt_y 		= (parameters.norm_y/np.linalg.norm(opt_y))*opt_y

	A 			= np.random.rand(m, n) - 0.5
	A 			= np.matmul(A, A.T) if parameters.make_psd == True else A

	u, s, v 	= np.linalg.svd(A, full_matrices=False)
	s 			= np.linspace(parameters.norm_A, parameters.norm_A/parameters.condition_number, min(m, n))

	A 			= np.dot(u * s, v) if parameters.symmetry == False else np.dot(u * s, u.T)

	b 			= np.matmul(A, opt_x)
	c 			= np.matmul(A.T, opt_y) + opt_s

	if parameters.norm_b != -1:
		temp_norm_b = np.linalg.norm(b)
		coef 		= parameters.norm_b / temp_norm_b
		opt_x 		= coef * opt_x
		b 			= coef * b

	results 	= (A, b, c, opt_x, opt_y, opt_s)

		
	return results




