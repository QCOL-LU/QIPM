import numpy as np
from copy import deepcopy
import sys

from ipm.general_methods.ParametersDefault import Parameters


#===========================================================================
# Linear optimization problem class 
#===========================================================================
def generate_problem(m, n, parameters):

	np.random.seed(parameters.seed)
	n 				= m if n == -1 else n

	#-----------------------------------------------------------------------
	# Initialize the problem
	#-----------------------------------------------------------------------
	random_vec		= np.random.rand(n)
	optimal_x 		= deepcopy(random_vec) 
	optimal_s 		= deepcopy(random_vec)
	optimal_y 		= np.random.rand(m) - .5


	#-----------------------------------------------------------------------
	# Construct an optimal solution (x, y, s)
	#-----------------------------------------------------------------------
	mask 			= [True if ind < m else False for ind in range(n)]
	np.random.shuffle(mask)

	if m != n:
		if mask[-1] == False and parameters.have_interior == True:
			true_ind 		= n - 1 - mask[::-1].index(True)
			mask[-1] 		= True
			mask[true_ind] 	= False		
		
		optimal_x[mask]						= 0
		optimal_s[np.logical_not(mask)]		= 0

	optimal_x 			= parameters.norm_x * optimal_x / np.linalg.norm(optimal_x)
	optimal_s 			= parameters.norm_s * optimal_s / np.linalg.norm(optimal_s)

	if  parameters.have_interior == True:
		norm_yy 		= np.linalg.norm(optimal_y[1:] )
		norm_yy  		= 1 if norm_yy == 0 else norm_yy
		optimal_y 		= np.sqrt(parameters.norm_y**2 - 1) * optimal_y/ norm_yy
		optimal_y[0] 	= 1
	else:
		optimal_y 		= parameters.norm_y * optimal_y/ np.linalg.norm(optimal_y)

	

	#-----------------------------------------------------------------------
	# Construct coefficient matrix 
	#-----------------------------------------------------------------------
	A 			= np.random.rand(m, n) - .5 
	A 			= np.matmul(A, A.T) if parameters.make_psd == True else A

	u, s, v 	= np.linalg.svd(A, full_matrices=False)
	s 			= np.linspace(parameters.norm_A, parameters.norm_A/parameters.condition_number, min(m, n))
	
	A 			= np.dot(u * s, v) if parameters.symmetry == False else np.dot(u * s, u.T)


	#-----------------------------------------------------------------------
	# Modify coefficient matrix to have an interior point
	#-----------------------------------------------------------------------
	if  parameters.have_interior == True:
		A[0, mask] 			= 0 
		A[0, np.logical_not(mask)] = 1

		A[:,-1] 			= 0

		last_col 			= A[:,mask]

		for (ind, elm) in enumerate(last_col):
			A[ind, -1] 		= A[ind, -1] - sum(elm)
			
	

	#-----------------------------------------------------------------------
	# Calculate the right hand side vector and objective function coefficient
	#-----------------------------------------------------------------------		
	b 			= np.matmul(A, optimal_x)
	c 			= np.matmul(A.T, optimal_y) + optimal_s


	#-----------------------------------------------------------------------
	# Determine an interior point
	#-----------------------------------------------------------------------
	if parameters.have_interior == True:
		int_x 				= deepcopy(optimal_x)
		int_y 				= deepcopy(optimal_y)
		int_y[0] 			= 0.0
		int_s 				= deepcopy(optimal_s)
		int_s[np.logical_not(mask)] 	= 1
		
		int_x[mask] 		= 1
		primal_res 			= b - np.matmul(A, int_x)

		dual_res   			= c - np.matmul(A.T, int_y) - int_s

		results 			= (A, b, c, int_x, int_y, int_s)

	else:
	
		if parameters.norm_b != -1:
			temp_norm_b		= np.linalg.norm(b)
			coef 			= parameters.norm_b / temp_norm_b
			optimal_x 		= coef * optimal_x
			b 				= coef * b

		results 			= (A, b, c)



	#-----------------------------------------------------------------------
	# Print the problem data
	#-----------------------------------------------------------------------
	if parameters.do_print == True:
		np.set_printoptions(precision=parameters.decimals, formatter={'float': ('{:+.'+str(parameters.decimals)+'f}').format})
		print()

		print(17*"=", "Linear optimization problem Information", 17*"=")
		print("{:<25}{:d}".format("Number of variables:",n) )
		print("{:<25}{:d}".format("Number of constraints:",m) )
		print()
		print("{:<25}{:.2e}".format("Norm of primal sol:",np.linalg.norm(optimal_x)) )
		print("{:<25}{:.2e}".format("Norm of dual sol:",np.linalg.norm(optimal_y)) )
		print("{:<25}{:.2e}".format("Norm of dual slacks:",np.linalg.norm(optimal_s)) )
		print()
		print("{:<25}{:.2e}".format("Norm of vector b:",np.linalg.norm(b)) )
		print("{:<25}{:.2e}".format("Norm of vector c:",np.linalg.norm(c)) )
		print()
		print("{:<25}{:.2e}".format("Norm of matrix A:",np.linalg.norm(A, 2)) )
		print("{:<25}{:.2e}".format("Condition number:",np.linalg.cond(A)) )
		print("{:<25}{:d}".format("Seed:",parameters.seed) )
		print()

		print(23*"=", "Linear optimization problem", 23*"=")
		print("{:<25}[{:}".format("coeff matrix A:", A[0]))
		for row in A[1:-1]:
			print("{:<25} {:}".format(" ", row))
		print("{:<25} {:}]".format(" ", A[-1]))
		print()
		print("{:<25}{:}".format("RHS vector b:", b))
		print("{:<25}{:}".format("Obj vector c:", c))
		print()
		
		print(28*"=", "Optimal solution", 29*"=")
		print("{:<25}".format("Primal variables:"), optimal_x)
		print("{:<25}".format("Dual slacks:"), optimal_s)
		print("{:<25}".format("Dual variables:"), optimal_y)
		print()
		print("{:<25}{:+.8e}".format("Optimal objective:", np.dot(c, optimal_x)) )
		print(75*"=")
		print()


	return results




