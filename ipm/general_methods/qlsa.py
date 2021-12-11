from qiskit import Aer
from qiskit.circuit.library import QFT
from qiskit.aqua import QuantumInstance
from qiskit.aqua.algorithms import HHL as old_HHL
# from qiskit.utils.algorithm_globals import HHL as old_HHL
from qiskit.aqua.components.eigs import EigsQPE
from qiskit.aqua.components.reciprocals import LookupRotation
from qiskit.aqua.operators import MatrixOperator
from qiskit.aqua.components.initial_states import Custom

from qiskit.algorithms.linear_solvers.numpy_linear_solver import NumPyLinearSolver
from qiskit.algorithms.linear_solvers.hhl import HHL
from qiskit.quantum_info import Statevector
from qiskit import BasicAer

from ipm.general_methods.ParametersDefault import Parameters

import numpy as np
import sys
import time






#===========================================================================
# Eigenvalues using Quantum Phase Estimation
#===========================================================================
def create_eigs(matrix, parameters, negative_evals):
	ne_qfts 			= [None, None]
	num_ancillae 		= parameters.num_ancillae
	if (negative_evals == True):
		num_ancillae 	= 1 + num_ancillae
		ne_qfts 		= [QFT(num_ancillae - 1), QFT(num_ancillae - 1).inverse()]
    
	max_eig 			= max(np.linalg.eigvals(matrix))

	eigenvalues 		=  EigsQPE( MatrixOperator(matrix= matrix),
								QFT(parameters.num_ancillae).inverse(),
								num_time_slices		= parameters.num_time_slices,
								num_ancillae		= num_ancillae,
								expansion_mode		= 'suzuki',
								expansion_order		= parameters.expansion_order,
								evo_time 			= np.pi/max_eig.real,  
								negative_evals 		= negative_evals,
								ne_qfts 			= ne_qfts)

	return eigenvalues



#===========================================================================
# Quantum linear system solver
#===========================================================================
def QLSA(cofficent_matrix, rhs_vector, parameters):

	#-----------------------------------------------------------------------
	# Initialize the paramters of the algorithm
	#-----------------------------------------------------------------------
	orig_size 			= len(rhs_vector)
	shots 				= 1e3
	seed_simulator 		= 541376
	
	#-----------------------------------------------------------------------
	# Resize the cofficent matrix 
	#-----------------------------------------------------------------------
	matrix, vector, truncate_powerdim, truncate_hermitian 	= old_HHL.matrix_resize(cofficent_matrix, rhs_vector)


	#-----------------------------------------------------------------------
	# Clear the QLSA message from the screen
	#-----------------------------------------------------------------------
	if truncate_powerdim:
		sys.stdout.write("\033[F")
		sys.stdout.write("\033[K")
		print(" ", end="\r")
	if truncate_hermitian:
		sys.stdout.write("\033[F")
		sys.stdout.write("\033[K")
		print(" ", end="\r")


		# Ax = b => A.T Ax = A.T b
	#-----------------------------------------------------------------------
	# Normalize the linear system based on the norm of RHS  
	#-----------------------------------------------------------------------
	normlize_coef 	= np.linalg.norm(vector)
	matrix 			= matrix / normlize_coef
	vector 			= vector / normlize_coef

	#-----------------------------------------------------------------------
	# Check whether the matrix has a negative eigenvalue
	#-----------------------------------------------------------------------
	min_eigval 		= min(np.linalg.eigvals(cofficent_matrix))
	has_neg_eigval 	= True if min_eigval < 0 else False


	start_time 		= time.time()
	
	#-----------------------------------------------------------------------
	# The old implementation of HHL algorithm
	#-----------------------------------------------------------------------
	if parameters.HHL_Method == 1:
		#-------------------------------------------------------------------
		# Initialize eigenvalue finding module
		#-------------------------------------------------------------------
		eigenvalues 	= create_eigs(matrix, parameters, has_neg_eigval)
		num_q, num_a 	= eigenvalues.get_register_sizes()


		#-------------------------------------------------------------------
		# Initialize initial state module
		#-------------------------------------------------------------------
		init_state 		= Custom(num_q, state_vector=vector)


		#-------------------------------------------------------------------
		# Initialize reciprocal rotation module
		#-------------------------------------------------------------------
		reciprocal 		= LookupRotation(negative_evals=eigenvalues._negative_evals, evo_time=eigenvalues._evo_time)

		algorithm 		= old_HHL(matrix, vector, truncate_powerdim, truncate_hermitian, eigenvalues,
								init_state, reciprocal, num_q, num_a, orig_size)
		

		#-------------------------------------------------------------------
		# Solve the linear system 
		#-------------------------------------------------------------------
		if(parameters.Is_Simulator == True):
			results 	= algorithm.run(QuantumInstance(Aer.get_backend('statevector_simulator'), shots=shots, seed_simulator=seed_simulator))
		else:
			print("Please add quantum computer setting.")

		for _ in range(10):
			sys.stdout.write("\033[F")
			sys.stdout.write("\033[K")
			print(" ", end="\r")

		solution 		= results['solution'].real
		probability 	= results['probability_result']


	#-----------------------------------------------------------------------
	# The new implementation of HHL algorithm
	#-----------------------------------------------------------------------
	elif parameters.HHL_Method == 2:

		backend 		= BasicAer.get_backend('statevector_simulator')
		hhl 			= HHL(parameters.qlsa_precision,  quantum_instance=backend)
		
		hhl_solution 	= hhl.solve(matrix, vector)
		naive_sv 		= np.real(Statevector(hhl_solution.state).data)

		ind 			= int(len(naive_sv)/2)
		ind2 			= 0 if len(rhs_vector) >= int(len(vector)/2) else int(len(vector)/2)

		solution 		=  np.array(naive_sv[ind: ind+len(vector)])[ind2: ind2+len(rhs_vector)]
	

	end_time 			= time.time()

	
	#---------------------------------------------------------------
	# Make sure that ||Ax|| == ||b||
	#---------------------------------------------------------------
	div  				= np.linalg.norm(np.dot(cofficent_matrix, solution))
	if div != 0:
		normlize_coef 	= np.linalg.norm(rhs_vector) / div
		solution 		= normlize_coef * solution



	#---------------------------------------------------------------
	# Checking the quantum_sol sign || rhs - M q_sol|| >= || rhs + M q_sol||
	#---------------------------------------------------------------
	is_sign_changed 	= False

	if rhs_vector.dot(cofficent_matrix).dot(solution) < 0:
		is_sign_changed	= True
		solution 		*= - 1
	
	exact_solution 		= np.linalg.solve(cofficent_matrix, rhs_vector)
	norm_of_difference 	= np.linalg.norm(solution - exact_solution)
	norm_of_residual 	= np.linalg.norm(rhs_vector - cofficent_matrix.dot(solution))

	#-----------------------------------------------------------------------
	# Print the results
	#-----------------------------------------------------------------------
	if parameters.do_print == True:

		np.set_printoptions(precision=3, formatter={'float': '{:+.3f}'.format})
		print()
		print(17*"-", "Linear system Information", 17*"-")
		print("{:<25}{:.2e}".format("Norm of RHS vector:", np.linalg.norm(rhs_vector)))
		print("{:<25}{:.2e}".format("Norm of matrix:", np.linalg.norm(cofficent_matrix, 2)))
		print("{:<25}{:.2e}".format("Condition number:", np.linalg.cond(cofficent_matrix)))
		print()
		print("{:<25}{:+.2e}".format("Minimum eigenvalue:", min_eigval))
		print()
		print("{:<25}{:}".format("RHS vector transpose:", rhs_vector))
		print()
		print("{:<25}[{:}".format("Coffcient matrix:", cofficent_matrix[0]), sep='')
		for row in cofficent_matrix[1:-1]:
			print("{:<25} {:}".format(" ", row), sep='')
		print("{:<25} {:}]".format(" ", cofficent_matrix[-1]), sep='')
		print()

		print(13*"-", "Modified linear system Information", 12*"-")
		print("{:<25}{:.2e}".format("Norm of RHS vector:", np.linalg.norm(vector)))
		print("{:<25}{:.2e}".format("Norm of matrix:", np.linalg.norm(matrix, 2)))
		print("{:<25}{:.2e}".format("Condition number:", np.linalg.cond(matrix)))
		print()
		print("{:<25}{:.2e}".format("Norm of difference:", norm_of_difference))
		print("{:<25}{:.2e}".format("Norm of residual:", norm_of_residual))

		
		print()
		print(19*"-", "Parameters of the QLSA", 18*"-")
		print("{:<25}{:d}".format("HHL method:", parameters.HHL_Method))
		if parameters.HHL_Method == 1:
			print("{:<25}{:d}".format("Num of ancillae qubits:", parameters.num_ancillae))
			print("{:<25}{:d}".format("Expansion order:", parameters.expansion_order))
			print("{:<25}{:d}".format("Num of time slices:", parameters.num_time_slices))
		elif parameters.HHL_Method == 2:
			print("{:<25}{:d}".format("num_qubits:", Statevector(hhl_solution.state).num_qubits))

		print("{:<25}{:}".format("Is simulator:", parameters.Is_Simulator))


		print()
		print(20*"-", "Results of the QLSA", 20*"-")

		print("{:<25}{:}".format("Quantum solution:", solution))
		print("{:<25}{:}".format("Exact solution:", exact_solution))
		print()
		print("{:<25}{:.2e}".format("Norm of difference:", norm_of_difference))
		print("{:<25}{:.2e}".format("Norm of residual:", norm_of_residual))
		print("{:<25}{:}".format("Is sign changed:", is_sign_changed))
		
		if parameters.HHL_Method == 1:
			print("{:<25}{:.2e}".format("Probability:", probability))

		print("{:<25}{:.2e}".format("Time (s):", end_time - start_time))
		
		
		if parameters.HHL_Method == 1:
			print()
			print("{:<25}{:d}".format("Circuit depth:", results['circuit_info']['depth']))
			print("{:<25}{:d}".format("Circuit width:", results['circuit_info']['width']))
		print(61*"=")

	return (solution, norm_of_residual, is_sign_changed)
	


