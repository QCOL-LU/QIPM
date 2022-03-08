import numpy as np
from copy import deepcopy
import sys
import mosek.fusion as mosek
import time
import threading


#===========================================================================
# Linear optimization problem class 
#===========================================================================
def generate_problem(m, n, parameters):

	np.random.seed(parameters.seed)

	#-----------------------------------------------------------------------
	# Generate a Linear Optimization problem
	#-----------------------------------------------------------------------
	if parameters.Problem_Type == "LO":
		#-------------------------------------------------------------------
		# Initialize the problem
		#-------------------------------------------------------------------
		random_vec		= np.random.rand(n)
		optimal_x 		= deepcopy(random_vec) 
		optimal_s 		= deepcopy(random_vec)
		optimal_y 		= np.random.rand(m) - .5


		#-------------------------------------------------------------------
		# Construct an optimal solution (x, y, s)
		#-------------------------------------------------------------------
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

		if parameters.have_interior == True:
			norm_yy 		= np.linalg.norm(optimal_y[1:] )
			norm_yy  		= 1 if norm_yy == 0 else norm_yy
			optimal_y 		= np.sqrt(parameters.norm_y**2 - 1) * optimal_y/ norm_yy
			optimal_y[0] 	= 1
		else:
			optimal_y 		= parameters.norm_y * optimal_y/ np.linalg.norm(optimal_y)

		

		#-------------------------------------------------------------------
		# Construct coefficient matrix 
		#-------------------------------------------------------------------
		A 			= np.random.rand(m, n) - .5 
		A 			= np.matmul(A, A.T) if parameters.make_psd == True else A

		u, s, v 	= np.linalg.svd(A, full_matrices=False)

		s 			= np.linspace(parameters.norm_A, parameters.norm_A/parameters.condition_number, min(m, n))

		A 			= np.dot(u * s, v) if parameters.symmetry == False else np.dot(u * s, u.T)
		

		#-------------------------------------------------------------------
		# Modify coefficient matrix to have an interior point
		#-------------------------------------------------------------------
		if  parameters.have_interior == True:
			A[0, mask] 			= 0 
			A[0, np.logical_not(mask)] = 1

			A[:,-1] 			= 0

			last_col 			= A[:,mask]

			for (ind, elm) in enumerate(last_col):
				A[ind, -1] 		= A[ind, -1] - sum(elm)
				
		

		#-------------------------------------------------------------------
		# Calculate the right hand side vector and objective function coefficient
		#-------------------------------------------------------------------		
		b 			= np.matmul(A, optimal_x)
		c 			= np.matmul(A.T, optimal_y) + optimal_s


		#-------------------------------------------------------------------
		# Determine an interior point
		#-------------------------------------------------------------------
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


			if parameters.norm_c != -1:
				temp_norm_c		= np.linalg.norm(c)
				coef 			= parameters.norm_c / temp_norm_c
				optimal_y 		= coef * optimal_y
				optimal_s 		= coef * optimal_s
				c 				= coef * c

			results 			= (A, b, c)


	#-----------------------------------------------------------------------
	# Generate a Semi-Definite Optimization problem 
	#-----------------------------------------------------------------------
	elif parameters.Problem_Type == "SDO":

		#-------------------------------------------------------------------
		# Construct an interior point x and s are 
		#-------------------------------------------------------------------
		int_x 					= np.random.rand(n, n) - 0.5
		int_x 					= (int_x + int_x.T) 
		int_x 					= int_x * parameters.norm_x / np.linalg.norm(int_x, 2)

		int_s 					= np.random.rand(n, n) - 0.5
		int_s 					= (int_s + int_s.T)
		int_s 					= int_s * parameters.norm_s / np.linalg.norm(int_s, 2)

		int_y 					= np.random.rand(m) - .5
		int_y 					= int_y * parameters.norm_y / np.linalg.norm(int_y)


		#-------------------------------------------------------------------
		# Construct coefficient matrices
		#-------------------------------------------------------------------
		coeff_matrices_A 		= [np.random.rand(n, n) - 0.5 for _ in range(m)]
		coeff_matrices_A 		= [matrix + matrix.T for matrix in coeff_matrices_A]
		coeff_matrices_A 		= [matrix * parameters.norm_A / np.linalg.norm(matrix, 2) for matrix in coeff_matrices_A]


		#-------------------------------------------------------------------
		# Calculate the right hand side matrix and objective function coefficient
		#-------------------------------------------------------------------
		b 						= np.array([np.trace(np.dot(matrix, int_x)) for matrix in coeff_matrices_A])
		c 						= sum([y * matrix for (matrix, y) in zip(coeff_matrices_A, int_y)]) + int_s

		results 				= (coeff_matrices_A, b, c, int_x, int_y, int_s)


		#-------------------------------------------------------------------
		# Construct the primal-dual model to obtain an optimal solution using MOSEK solver
		#-------------------------------------------------------------------
		model 		= mosek.Model("SDO problem")
		X 			= model.variable(mosek.Domain.isTrilPSD(n))
		S 			= model.variable(mosek.Domain.isTrilPSD(n))
		y 			= model.variable(m)




		model.objective(mosek.ObjectiveSense.Minimize,mosek.Expr.sub(mosek.Expr.dot(c, X), mosek.Expr.dot(b, y)) )

		for (matrix, rhs) in zip(coeff_matrices_A, b):
			model.constraint(mosek.Expr.dot( matrix, X), mosek.Domain.equalsTo(rhs))

		exper 		= 0
		for (ind, matrix) in enumerate(coeff_matrices_A):
			exper 	= mosek.Expr.add(exper, mosek.Expr.mul(matrix, y.index(ind))) 

		model.constraint(mosek.Expr.add(exper, S), mosek.Domain.equalsTo(c))



		#-------------------------------------------------------------------
		# Solve the problem using MOSEK solver
		#-------------------------------------------------------------------
		T 			= threading.Thread(target=model.solve)

		start_time 	= time.time()

		try:
			T.start() 			# optimization now running in background

			# Loop until we get a solution or you run out of patience and press Ctrl-C
			while True:
				if not T.is_alive():
					break
				elif time.time() - start_time > parameters.time_limit:
					print("Mosek solver terminated due to time_limit!")
					model.breakSolver()
					break
		except KeyboardInterrupt:
			print("Signalling the solver that it can give up now!")
			model.breakSolver()
		finally:
			try: T.join() 		# wait for the solver to return
			except: pass


		optimal_x 	= np.array([[X.index(i, j).level()[0] for j in range(n)] for i in range(n)])
		optimal_s	= np.array([[S.index(i, j).level()[0] for j in range(n)] for i in range(n)])
		optimal_y 	= np.array([y.index(i).level()[0] for i in range(m)])


 
	#-----------------------------------------------------------------------
	# Print the problem data
	#-----------------------------------------------------------------------
	if parameters.do_print == True:
		np.set_printoptions(precision=parameters.decimals, formatter={'float': ('{:+.'+str(parameters.decimals)+'f}').format})
		print()

		if parameters.Problem_Type == "LO":
			print(15*"-", "Linear optimization problem characteristics", 15*"-")
		elif parameters.Problem_Type == "SDO":
			print(12*"-", "Semi-Definite optimization problem characteristics", 11*"-")


		print("{:<25}{:d}".format("Seed:", parameters.seed) )
		print("{:<25}{:d}".format("Number of variables:",n) )
		print("{:<25}{:d}".format("Number of constraints:",m) )
		print()
		print("{:<25}{:.2e}".format("Norm of primal sol:",np.linalg.norm(optimal_x)) )
		print("{:<25}{:.2e}".format("Norm of dual sol:",np.linalg.norm(optimal_y)) )
		print("{:<25}{:.2e}".format("Norm of dual slacks:",np.linalg.norm(optimal_s)) )
		print()
		print("{:<25}{:.2e}".format("Norm of vector b:",np.linalg.norm(b)) )
		
		
		if parameters.Problem_Type == "LO":
			print("{:<25}{:.2e}".format("Norm of vector c:",np.linalg.norm(c)) )
			print()
			print("{:<25}{:.2e}".format("Norm of matrix A:",np.linalg.norm(A, 2)) )
			print("{:<25}{:.2e}".format("Condition number:",np.linalg.cond(A)) )
			print()

			print(23*"-", "Linear optimization problem", 23*"-")
			print("{:<25}[{:}".format("Coeff matrix A:", A[0]))
			for row in A[1:-1]:
				print("{:<25} {:}".format(" ", row))
			print("{:<25} {:}]".format(" ", A[-1]))
			print()

			print("{:<25}{:}".format("RHS vector b:", b))
			print("{:<25}{:}".format("Obj vector c:", c))
			print()
		
			print(28*"-", "Optimal solution", 29*"-")
			print("{:<25}".format("Primal variables:"), optimal_x)
			print("{:<25}".format("Dual slacks:"), optimal_s)
			print("{:<25}".format("Dual variables:"), optimal_y)
			print()
			print("{:<25}{:+.8e}".format("Optimal objective:", np.dot(c, optimal_x)) )
			
		
		elif parameters.Problem_Type == "SDO":
			print("{:<25}{:.2e}".format("Norm of C:",np.linalg.norm(c)) )
			print()

			print(20*"-", "Semi-Definite optimization problem", 19*"-")
			for (ind, matrix) in enumerate(coeff_matrices_A):
				print("Norm of matrix A{:<9}{:.2e}".format(str(ind) + ":",np.linalg.norm(matrix, 2)) )
				print("Condition number{:<9}{:.2e}".format(str(ind) + ":", np.linalg.cond(matrix)) )
				print("Coeff matrix A{:<11}[{:}".format(str(ind) + ":", matrix[0]))
				for row in matrix[1:-1]:
					print("{:<25} {:}".format(" ", row))
				print("{:<25} {:}]".format(" ", matrix[-1]))
				print()

			print("{:<25}{:}".format("RHS vector b:", b))
			print()

			print("matrix C{:<17}[{:}".format(":", c[0]))
			for row in c[1:-1]:
				print("{:<25} {:}".format(" ", row))
			print("{:<25} {:}]".format(" ", c[-1]))
			print()


			print(28*"-", "Optimal solution", 29*"-")

			print("Primal variables{:<9}[{:}".format(":", optimal_x[0]))
			for row in optimal_x[1:-1]:
				print("{:<25} {:}".format(" ", row))
			print("{:<25} {:}]".format(" ", optimal_x[-1]))
			print()

			print("Dual slacks{:<14}[{:}".format(":", optimal_s[0]))
			for row in optimal_s[1:-1]:
				print("{:<25} {:}".format(" ", row))
			print("{:<25} {:}]".format(" ", optimal_s[-1]))
			print()

			print("{:<25}".format("Dual variables:"), optimal_y)
			print()
			print("{:<25}{:+.8e}".format("Optimal objective:", np.trace(np.dot(c, optimal_x))) )
			

		print(75*"=")
		print()
	return results




