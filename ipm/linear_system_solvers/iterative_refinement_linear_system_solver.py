import numpy as np
import sys
import time


from ipm.general_methods.ParametersDefault import Parameters
from ipm.linear_system_solvers.qlsa import *


#===========================================================================
# Quantum linear system solver
#===========================================================================
def iterative_refinement_linear_system_solver(cofficent_matrix, rhs_vector, parameters):

	#-----------------------------------------------------------------------
	# Initialize the paramters of the algorithm
	#-----------------------------------------------------------------------
    nabla             	= parameters.LS_ScalFact                    # Scaling factor (1e2)
    rho             	= parameters.LS_IncScalLim                  # Incremental scaling

    IRprecision         = parameters.IR_LS_Precision


    start_time          = time.time() 

    d                   = len(cofficent_matrix)
    iteration           = 0
    solution            = np.zeros(d)
    r                   = rhs_vector
    res                 = []

    failed              = False

    def noise():        
        noise_vector    = np.random.random(d)
        noise_vector    = noise_vector / np.linalg.norm(noise_vector) * parameters.LS_Precision
        return noise_vector

    if parameters.do_print:
        print("The IR-LS algorithm has started.\n")
        if parameters.Is_Quantum == False:
            print("{:}{:7}{:}{:9}{:}".format("iter", "","||r||", "","nabla"))
        else:
            print("{:}{:5}{:}{:7}{:}{:9}{:}{:9}{:}".format("iter", "","diff-QLSA", "","||r||", "","nabla", "", "time"))


    while (np.linalg.norm(r) > IRprecision):
        old_r           = r
        if parameters.Is_Quantum == False:
            c           = np.linalg.solve(cofficent_matrix, nabla*r) + (noise() if parameters.Is_Noisy else np.zeros(d))
        else:
            c           = QLSA(cofficent_matrix, nabla*r, parameters)[0]
            exact_solution  = np.linalg.solve(cofficent_matrix, nabla*r)
            QLSA_Precision  = np.linalg.norm(c - exact_solution) / np.linalg.norm(exact_solution)

        solution        = solution + (1/nabla)*c
        r               = rhs_vector - np.dot(cofficent_matrix, solution)
        nabla           = min(rho*nabla, 1/(np.linalg.norm(r)))
        res.append(np.linalg.norm(r))
        iteration       = iteration + 1
        runtime         = time.time() - start_time

        if parameters.Is_Quantum == False:
            print("{:>4d}{:5}{:.3e}{:5}{:.3e}{:5}{:>5.0f}s".format(iteration, " ",np.linalg.norm(r), " ",nabla, "", runtime))
        else:
            print("{:>4d}{:5}{:.3e}{:5}{:.3e}{:5}{:.3e}{:5}{:>5.0f}s".format(iteration," " ,QLSA_Precision , " ",np.linalg.norm(r), " ",nabla, "", runtime))

        if np.linalg.norm(old_r - r) < 1e-8 and parameters.do_print: 
            failed      = True
            print(61*"-")
            print("\nThe algorithm terminated before reaching the desired precision.")
            break 

    if failed == False and parameters.do_print: 
        print(61*"-")
        print("\nThe algorithm successfully terminated.")

    if parameters.Is_Quantum == True: parameters.LS_Precision = QLSA_Precision
    return solution, parameters.LS_Precision, iteration
	


