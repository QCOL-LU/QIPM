import numpy as np
from copy import deepcopy
import sys
import mosek.fusion as mosek
import time
import threading


from generate_lo_problem_with_int_and_opt import *
from generate_lo_problem_with_opt import *
from generate_lo_problem_with_int import *

from generate_sdo_problem_with_int_and_opt import *
from generate_sdo_problem_with_opt import *
from generate_sdo_problem_with_int import *

from generate_soco_problem_with_int_and_opt import *
from generate_soco_problem_with_opt import *
from generate_soco_problem_with_int import *

#===========================================================================
# Generate conic optimization problems
#===========================================================================
def generate_problem(m, n, parameters):

	if parameters.has_interior == False and parameters.has_optimal == False: 
		sys.exit("Both has_interior and has_optimal cannot be \"False\".")

	np.random.seed(parameters.seed)

	#-----------------------------------------------------------------------
	# Generate a Linear Optimization problems
	#-----------------------------------------------------------------------
	if parameters.Problem_Type == "LO":
		n 				= m if n == -1 else n

		if parameters.has_interior == True and parameters.has_optimal == True:
			if m < 1 or n < 1:
				sys.exit("m and n must be larger than 1.")

			results 	= generate_lo_problem_with_int_and_opt(m - 1, n - 1, parameters)

		elif parameters.has_interior == False and parameters.has_optimal == True:
			results 	= generate_lo_problem_with_opt(m, n, parameters)

		elif parameters.has_interior == True and parameters.has_optimal == False:
			results 	= generate_lo_problem_with_int(m, n, parameters)


	#-----------------------------------------------------------------------
	# Generate a Semidefinite Optimization problems
	#-----------------------------------------------------------------------
	elif parameters.Problem_Type == "SDO":
		if parameters.has_interior == True and parameters.has_optimal == True:
			results 	= generate_sdo_problem_with_int_and_opt(m, n, parameters)

		elif parameters.has_interior == False and parameters.has_optimal == True:
			results 	= generate_sdo_problem_with_opt(m, n, parameters)

		elif parameters.has_interior == True and parameters.has_optimal == False:
			results 	= generate_sdo_problem_with_int(m, n, parameters)


	#-----------------------------------------------------------------------
	# Generate a Second Order Cone Optimization problems
	#-----------------------------------------------------------------------
	elif parameters.Problem_Type == "SOCO":

		if parameters.has_interior == True and parameters.has_optimal == True:
			results 	= generate_soco_problem_with_int_and_opt(m, n, parameters)

		elif parameters.has_interior == False and parameters.has_optimal == True:
			results 	= generate_soco_problem_with_opt(m, n, parameters)

		elif parameters.has_interior == True and parameters.has_optimal == False:
			results 	= generate_soco_problem_with_int(m, n, parameters)	

	#-----------------------------------------------------------------------
	# Print the problem data
	#-----------------------------------------------------------------------
	if parameters.do_print == True:
		np.set_printoptions(precision=parameters.decimals, formatter={'float': ('{:+.'+str(parameters.decimals)+'f}').format})
		print()

		if parameters.Problem_Type == "LO":
			print(15*"-", "Linear optimization problem characteristics", 15*"-")
		elif parameters.Problem_Type == "SDO":
			print(12*"-", "Semidefinite optimization problem characteristics", 11*"-")


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

			print(20*"-", "Semidefinite optimization problem", 19*"-")
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




