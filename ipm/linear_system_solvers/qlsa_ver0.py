from qiskit import Aer
from qiskit.circuit.library import QFT
from qiskit.aqua import QuantumInstance
from qiskit.aqua.algorithms import HHL
from qiskit.aqua.components.eigs import EigsQPE
from qiskit.aqua.components.reciprocals import LookupRotation
from qiskit.aqua.operators import MatrixOperator
from qiskit.aqua.components.initial_states import Custom

import numpy as np
import sys
import time


#===========================================================================
# Eigenvalues using Quantum Phase Estimation
#===========================================================================
def create_eigs(matrix, num_ancillae, num_time_slices, negative_evals, expansion_order):
	ne_qfts 			= [None, None]

	if (negative_evals == True):
		num_ancillae 	+= 1
		ne_qfts 		= [QFT(num_ancillae - 1), QFT(num_ancillae - 1).inverse()]
    
	max_eig 			= max(np.linalg.eigvals(matrix))

	eigenvalues 		=  EigsQPE( MatrixOperator(matrix= matrix),
								QFT(num_ancillae).inverse(),
								num_time_slices		= num_time_slices,
								num_ancillae		= num_ancillae,
								expansion_mode		= 'suzuki',
								expansion_order		= expansion_order,
								evo_time 			= np.pi/max_eig.real,  
								negative_evals 		= negative_evals,
								ne_qfts 			= ne_qfts)

	return eigenvalues



#===========================================================================
# Quantum linear system solver
#===========================================================================
def QLSA(cofficent_matrix, rhs_vector, num_ancillae=7, num_time_slices=1, expansion_order=2, is_simulator=True, do_print=False):


	start_time 		= time.time()
	#-----------------------------------------------------------------------
	# Initialize the paramters of the algorithm
	#-----------------------------------------------------------------------
	orig_size 		= len(rhs_vector)
	shots 			= 1e3
	seed_simulator 	= 541376

	#-----------------------------------------------------------------------
	# Check whether the matrix has a negative eigenvalue
	#-----------------------------------------------------------------------
	min_eigval 		= min(np.linalg.eigvals(cofficent_matrix))
	has_neg_eigval 	= True if min_eigval < 0 else False


	#-----------------------------------------------------------------------
	# Resize the cofficent matrix 
	#-----------------------------------------------------------------------
	matrix, vector, truncate_powerdim, truncate_hermitian 	= HHL.matrix_resize(cofficent_matrix, rhs_vector)


	#-----------------------------------------------------------------------
	# Normalize the linear system based on the norm of RHS 
	#-----------------------------------------------------------------------
	normlize_coef 	= np.linalg.norm(vector)
	matrix 			= matrix / normlize_coef
	vector 			= vector / normlize_coef

	# normlize_coef 	= np.linalg.norm(matrix, 2)
	# matrix 			= matrix / normlize_coef
	# vector 			= vector / normlize_coef

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

	#-----------------------------------------------------------------------
	# Initialize eigenvalue finding module
	#-----------------------------------------------------------------------
	eigenvalues 	= create_eigs(matrix, num_ancillae, num_time_slices, has_neg_eigval, expansion_order)
	num_q, num_a 	= eigenvalues.get_register_sizes()


	#-----------------------------------------------------------------------
	# Initialize initial state module
	#-----------------------------------------------------------------------
	init_state 		= Custom(num_q, state_vector=vector)


	#-----------------------------------------------------------------------
	# Initialize reciprocal rotation module
	#-----------------------------------------------------------------------
	reciprocal 		= LookupRotation(negative_evals=eigenvalues._negative_evals, evo_time=eigenvalues._evo_time)

	algorithm 		= HHL(matrix, vector, truncate_powerdim, truncate_hermitian, eigenvalues,
							init_state, reciprocal, num_q, num_a, orig_size)
	

	#-----------------------------------------------------------------------
	# Solve the linear system 
	#-----------------------------------------------------------------------
	if(is_simulator == True):
		results 	= algorithm.run(QuantumInstance(Aer.get_backend('statevector_simulator'), shots=shots, seed_simulator=seed_simulator))
	else:
		print("Please add quantum computer setting.")


	solution 		= results['solution'].real
	m_solution 		= solution * normlize_coef
	probability 	= results['probability_result']

	

	#---------------------------------------------------------------
	# Checking the quantum_sol sign || rhs - M q_sol|| >= || rhs + M q_sol||
	#---------------------------------------------------------------
	is_sign_changed 	= False

	if rhs_vector.dot(cofficent_matrix).dot(m_solution) < 0:
		is_sign_changed	= True
		m_solution 		*= - 1
		solution 		*= - 1

	m_exact_solution 		= np.linalg.solve(cofficent_matrix, rhs_vector)
	m_norm_of_difference 	= np.linalg.norm(m_solution - m_exact_solution)
	m_norm_of_residual 		= np.linalg.norm(rhs_vector - cofficent_matrix.dot(m_solution))


	#-----------------------------------------------------------------------
	# Print the results
	#-----------------------------------------------------------------------
	if do_print == True:

		exact_solution 		= np.linalg.solve(cofficent_matrix, rhs_vector)
		norm_of_difference 	= np.linalg.norm(solution - exact_solution)
		norm_of_residual 	= np.linalg.norm(rhs_vector - cofficent_matrix.dot(solution))

		end_time 			= time.time()
		np.set_printoptions(precision=3, formatter={'float': '{:+.3f}'.format})
		print()
		print(17*"=", "Linear system Information", 17*"=")
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

		print(13*"=", "Modified linear system Information", 12*"=")
		print("{:<25}{:.2e}".format("Norm of RHS M-vector:", np.linalg.norm(vector)))
		print("{:<25}{:.2e}".format("Norm of M-matrix:", np.linalg.norm(matrix, 2)))
		print("{:<25}{:.2e}".format("M-condition number:", np.linalg.cond(matrix)))
		print()
		print("{:<25}{:.2e}".format("M-norm of difference:", norm_of_difference))
		print("{:<25}{:.2e}".format("M-norm of residual:", norm_of_residual))

		
		print()
		print(19*"=", "Parameters of the QLSA", 18*"=")
		print("{:<25}{:d}".format("Num of ancillae qubits:", num_ancillae))
		print("{:<25}{:d}".format("Expansion order:", expansion_order))
		print("{:<25}{:d}".format("Num of time slices:", num_time_slices))
		print("{:<25}{:}".format("Is simulator:", is_simulator))

		print()
		print(20*"=", "Results of the QLSA", 20*"=")

		print("{:<25}{:}".format("Quantum solution:", m_solution))
		print("{:<25}{:}".format("Exact solution:", m_exact_solution))
		print()
		print("{:<25}{:.2e}".format("Norm of difference:", m_norm_of_difference))
		print("{:<25}{:.2e}".format("Norm of residual:", m_norm_of_residual))
		print("{:<25}{:}".format("Is sign changed:", is_sign_changed))
		print("{:<25}{:.2e}".format("Probability:", probability))
		print("{:<25}{:.2e}".format("Time (s):", end_time - start_time))
		
		print()
		print("{:<25}{:d}".format("Circuit depth:", results['circuit_info']['depth']))
		print("{:<25}{:d}".format("Circuit width:", results['circuit_info']['width']))
		print(61*"=")

	return (m_solution, probability, m_norm_of_residual, is_sign_changed)
	


