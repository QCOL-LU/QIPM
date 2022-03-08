import sys
sys.path.insert(1, '..')
from ipm import *
import numpy as np
import sys


parameters 	= Parameters()

parameters.seed  		= int(sys.argv[1])
parameters.norm_rhs 	= 2
parameters.cond 		= float(sys.argv[2])

A, b, c, int_x, int_y, int_s	= generate_problem(m=3, n=4, parameters=parameters)

model 		= Model(A, b, c)
model.x 	= int_x
model.y 	= int_y
model.s 	= int_s
model.Params.HHL_Method 		= 2

model.Params.Method 		= "IF-QIPM"
model.Params.LO_Precision	= 1e-1
model.Params.IR_Precision 	= 1e-2
model.Params.LO_Verbosity 	= 2
model.Params.IR_Verbosity 	= 2
model.Params.Omega 			= 1e1
model.Params.Stop_Precision = 1e-3

model.Params.num_ancillae 	= 4
model.Params.num_time_slices= 1
model.Params.expansion_order= 2

# model.Params.Beta_2			= 1e-2
# model.Params.Beta_1			= 0.5



model.solve()