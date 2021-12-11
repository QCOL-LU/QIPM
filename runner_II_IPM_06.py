from ipm import *
import numpy as np
import sys

parameters 					= Parameters()


parameters.seed  			= int(sys.argv[1])
parameters.norm_b 			= 2

parameters.have_interior 	= False
parameters.do_print 		= True


A, b, c						= generate_problem(m=2, n=4, parameters=parameters)

model 						= Model(A, b, c)

model.Params.HHL_Method 	= 2

model.Params.Method 		= "II-IPM"
model.Params.LO_Precision	= 1e-1
model.Params.IR_Precision 	= 1e-2

model.Params.LO_Verbosity 	= 2
model.Params.IR_Verbosity 	= 2

model.Params.Omega 			= float(sys.argv[2])
model.Params.Stop_Precision = 1e-8
model.Params.Stop_Cond_Num 	= 1e10
model.Params.qlsa_precision = 1e0


# model.Params.num_ancillae 	= 4
# model.Params.num_time_slices= 1
# model.Params.expansion_order= 2

# model.Params.Beta_2			= 1e-2
# model.Params.Beta_1			= 0.5



model.solve()