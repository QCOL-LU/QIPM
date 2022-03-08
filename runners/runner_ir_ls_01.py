import sys
sys.path.insert(1, '..')
from ipm import *
import sys


parameters 					= Parameters()

parameters.seed  			= int(sys.argv[1])
parameters.norm_b 			= 4
parameters.norm_A 			= 4
parameters.condition_number = 20
	
parameters.have_interior 	= False
parameters.Problem_Type 	="LO"
parameters.symmetry 		= True


# A, b, c 	= generate_problem(m=3, norm_x= 1, seed=seed, condition_number=cond,do_print=False, symmetry=True, have_interior=False)

A, b, c	= generate_problem(m=4, n=4, parameters=parameters)


parameters.LS_Method 			= "IR-LS"
parameters.Is_Quantum 			= True
parameters.Is_Noisy 			= True
parameters.IR_LS_Precision 		= float(sys.argv[2])
parameters.do_print				= True
parameters.is_simulator			= True


parameters.num_ancillae			= 5
parameters.expansion_order		= 1
parameters.symmertric_method 	= 1
parameters.HHL_Method 			= 2
# parameters.qlsa_precision 		= 1e-3

linear_system_solver(A, b, parameters)
