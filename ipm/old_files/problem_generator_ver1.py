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
	# LO problem with Interior 
	#-----------------------------------------------------------------------
	if parameters.have_interior == True and parameters.have_optimal == False:
		norm_x = 10.
		norm_s = 10.
		norm_y = 10.
		int_x = np.random.rand(n)
		int_s = np.random.rand(n)
		int_y = np.random.rand(m) - 0.5
		int_x = (norm_x/np.linalg.norm(int_x))*int_x
		int_s = (norm_s/np.linalg.norm(int_s))*int_s
		int_y = (norm_y/np.linalg.norm(int_y))*int_y
		A = np.random.rand(m, n) - .5 
		A = np.matmul(A, A.T) if parameters.make_psd == True else A

		u, s, v 	= np.linalg.svd(A, full_matrices=False)
		s = np.linspace(parameters.norm_A, parameters.norm_A/parameters.condition_number, min(m, n))
	
		A = np.dot(u * s, v) if parameters.symmetry == False else np.dot(u * s, u.T)

		b = np.matmul(A, int_x)
		c = np.matmul(A.T, int_y) + int_s
		if parameters.norm_b != -1:
			temp_norm_b = np.linalg.norm(b)
			coef = parameters.norm_b / temp_norm_b
			int_x = coef * int_x
			b = coef * b

		results 	= (A, b, c, int_x, int_y, int_s)
	elif parameters.have_interior == False and parameters.have_optimal == True:

        norm_x = 10.
		norm_s = 10.
		norm_y = 10.
        B=m
        mask             = [1 if ind < B else 0 for ind in range(n)]
        np.random.shuffle(mask)
        opt_x = np.random.rand(n)
        opt_x = np.multiply(opt_x,mask)
		opt_s = np.random.rand(n)
        opt_s = opt_s-np.multiply(opt_s,mask)
		opt_y = np.random.rand(m) - 0.5
		opt_x = (norm_x/np.linalg.norm(opt_x))*opt_x
		opt_s = (norm_s/np.linalg.norm(opt_s))*opt_s
		opt_y = (norm_y/np.linalg.norm(opt_y))*opt_y
		A = np.random.rand(m, n) - 0.5
		A = np.matmul(A, A.T) if parameters.make_psd == True else A

		u, s, v 	= np.linalg.svd(A, full_matrices=False)
		s = np.linspace(parameters.norm_A, parameters.norm_A/parameters.condition_number, min(m, n))
	
		A = np.dot(u * s, v) if parameters.symmetry == False else np.dot(u * s, u.T)

		b = np.matmul(A, opt_x)
		c = np.matmul(A.T, opt_y) + opt_s
		if parameters.norm_b != -1:
			temp_norm_b = np.linalg.norm(b)
			coef = parameters.norm_b / temp_norm_b
			opt_x = coef * opt_x
			b = coef * b

		results 	= (A, b, c, opt_x, opt_y, opt_s)
    elif parameters.have_interior == True and parameters.have_optimal == True:
        norm_x = 10.
        norm_s = 10.
        norm_y = 10.
        B=m
        mask= [1 if ind < B else 0 for ind in range(n)]
        np.random.shuffle(mask)
        opt_x = np.random.rand(n)
        opt_x = np.multiply(opt_x,mask)
        opt_s = np.random.rand(n)
        opt_s = opt_s-np.multiply(opt_s,mask)
        opt_y = np.random.rand(m) - 0.5
        A = np.random.rand(m, n) - 0.5
        int_x = np.random.rand(n+1)
        int_s = np.random.rand(n+1)
        int_y = np.random.rand(m+1) - 0.5
        int_y[m]= 0.5
        delta= np.dot((int_x[:n]-opt_x).T,(int_s[:n]-opt_s))
        de=max(0,-delta/int_x[n])
        int_s[n] = de if int_s[n]<= de else int_s[n]
        opt_x = np.append(opt_x,0)
        opt_s = np.append(opt_s,(delta/int_x[n])+int_s[n])
        opt_y = np.append(opt_y,0)
        opt_x = (norm_x/np.linalg.norm(opt_x))*opt_x
        opt_s = (norm_s/np.linalg.norm(opt_s))*opt_s
        opt_y = (norm_y/np.linalg.norm(opt_y))*opt_y
        int_x = (norm_x/np.linalg.norm(int_x))*int_x
        int_s = (norm_s/np.linalg.norm(int_s))*int_s
        int_y = (norm_y/np.linalg.norm(int_y))*int_y
        ahat=(1/int_x[n])*(np.dot(A,(opt_x[:n]-int_x[:n])))
        dhat=(1/int_y[m])*(np.dot(A.T,(opt_y[:m]-int_y[:m]))+opt_s[:n]-int_s[:n])
        dlast=(1/int_x[n])*(np.dot(dhat,(opt_x[:n]-int_x[:n])))
        AA=np.blocks([[A,np.matrix(ahat).T],[np.matrix(dhat),np.matrix(dlast)]])
        
        AA = np.matmul(AA, AA.T) if parameters.make_psd == True else AA

        u, s, v     = np.linalg.svd(A, full_matrices=False)
        s = np.linspace(parameters.norm_A, parameters.norm_A/parameters.condition_number, min(m, n))
    
        AA = np.dot(u * s, v) if parameters.symmetry == False else np.dot(u * s, u.T)

        b = np.matmul(AA, opt_x)
        c = np.matmul(AA.T, opt_y) + opt_s
        if parameters.norm_b != -1:
            temp_norm_b = np.linalg.norm(b)
            coef = parameters.norm_b / temp_norm_b
            opt_x = coef * opt_x
            int_x = coef * int_x
            b = coef * b
            
        results     = (A, b, c, opt_x, opt_y, opt_s, int_x, int_y, int_s)

	else: print("Wrong input!")
        
	if parameters.do_print == True:
		np.set_printoptions(precision=parameters.decimals, formatter={'float': ('{:+.'+str(parameters.decimals)+'f}').format})
		print()

		print(17*"=", "Linear optimization problem Information", 17*"=")
		print("{:<25}{:d}".format("Number of variables:",n) )
		print("{:<25}{:d}".format("Number of constraints:",m) )
		print()
		print("{:<25}{:.2e}".format("Norm of primal sol:",np.linalg.norm(opt_x)) )
		print("{:<25}{:.2e}".format("Norm of dual sol:",np.linalg.norm(opt_y)) )
		print("{:<25}{:.2e}".format("Norm of dual slacks:",np.linalg.norm(opt_s)) )
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
		print("{:<25}".format("Primal variables:"), opt_x)
		print("{:<25}".format("Dual slacks:"), opt_s)
		print("{:<25}".format("Dual variables:"), opt_y)
		print()
		print("{:<25}{:+.8e}".format("Optimal objective:", np.dot(c, opt_x)) )
		print(75*"=")
		print()


	return results




