import numpy as np
from copy import deepcopy



#===========================================================================
# Linear optimization problem with known interior and optimal solutions
#===========================================================================
def generate_lo_problem_with_int_and_opt(m, n, parameters):

	np.random.seed(parameters.seed)	

	mask 			= [1 if ind < m else 0 for ind in range(n)]
	np.random.shuffle(mask)

	opt_x 			= np.random.rand(n)
	opt_x 			= np.multiply(opt_x,mask)
	opt_s 			= np.random.rand(n)
	opt_s 			= opt_s-np.multiply(opt_s,mask)
	opt_y 			= np.random.rand(m) - 0.5
	A 				= np.random.rand(m, n) - 0.5

	int_x 			= np.random.rand(n+1)
	int_s 			= np.random.rand(n+1)
	int_y 			= np.random.rand(m+1) - 0.5
	int_y[m]		= 0.5

	delta			= np.dot((int_x[:n]-opt_x).T,(int_s[:n]-opt_s))
	de 				= max(0,-delta/int_x[n])
	int_s[n] 		= de if int_s[n] <= de else int_s[n]

	opt_x 			= np.append(opt_x,0)
	opt_s 			= np.append(opt_s,(delta/int_x[n])+int_s[n])
	opt_y 			= np.append(opt_y,0)

	opt_x 			= (parameters.norm_x/np.linalg.norm(opt_x))*opt_x
	opt_s 			= (parameters.norm_s/np.linalg.norm(opt_s))*opt_s
	opt_y 			= (parameters.norm_y/np.linalg.norm(opt_y))*opt_y

	int_x 			= (parameters.norm_x/np.linalg.norm(int_x))*int_x
	int_s 			= (parameters.norm_s/np.linalg.norm(int_s))*int_s
	int_y 			= (parameters.norm_y/np.linalg.norm(int_y))*int_y

	ahat 			= (1/int_x[n])*(np.dot(A,(opt_x[:n]-int_x[:n])))
	dhat 			= (1/int_y[m])*(np.dot(A.T,(opt_y[:m]-int_y[:m]))+opt_s[:n]-int_s[:n])
	dlast 			= (1/int_x[n])*(np.dot(dhat,(opt_x[:n]-int_x[:n])))
	AA 				= np.blocks([[A,np.matrix(ahat).T],[np.matrix(dhat),np.matrix(dlast)]])

	AA 				= np.matmul(AA, AA.T) if parameters.make_psd == True else AA

	u, s, v 		= np.linalg.svd(A, full_matrices=False)
	s 				= np.linspace(parameters.norm_A, parameters.norm_A/parameters.condition_number, min(m, n))

	AA 				= np.dot(u * s, v) if parameters.symmetry == False else np.dot(u * s, u.T)

	b 				= np.matmul(AA, opt_x)
	c 				= np.matmul(AA.T, opt_y) + opt_s

	if parameters.norm_b != -1:
		temp_norm_b = np.linalg.norm(b)
		coef 		= parameters.norm_b / temp_norm_b
		opt_x 		= coef * opt_x
		int_x 		= coef * int_x
		b 			= coef * b
	    
	results     	= (A, b, c, opt_x, opt_y, opt_s, int_x, int_y, int_s)

		
	return results




