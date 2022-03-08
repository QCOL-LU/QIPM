import sys
sys.path.insert(1, '..')
from ipm import *
import sys


parameters 						= Parameters()

parameters.seed  				= int(sys.argv[1])
dimension 						= int(sys.argv[2])
parameters.norm_b 				= 2
parameters.condition_number		= 4
	
parameters.do_print 			= True

parameters.have_interior 		= False
parameters.Problem_Type 		="LO"
parameters.symmetry 			= True

A, b, c	= generate_problem(m=dimension, n=dimension, parameters=parameters)



parameters.do_print				= True
parameters.is_simulator			= True
parameters.num_ancillae			= 5
parameters.expansion_order		= 1
parameters.symmertric_method 	= 1
parameters.HHL_Method 			= 2

QLSA(A, b, parameters)