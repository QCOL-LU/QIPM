from ipm.linear_system_solvers.linear_system_solvers import *
from ipm.print_methods.print_IPM import *

import numpy as np
from time import time
from copy import deepcopy
import sys




#===========================================================================
# Inexact Infeasible Interior Point Method
#===========================================================================
def inexact_feasible_IPM(LO):
	#-----------------------------------------------------------------------
	# Initialize the paramters of the algorithm
	#-----------------------------------------------------------------------
	beta_1 			= 0.5	#LO.Params.Beta_1					# Default value: 0.1
	eta 			= 0.5
	precision 		= LO.Params.LO_Precision

	iteration 		= 0
	
	start_time 		= time()
	

	if LO.Params.LO_Verbosity > -1:
		if LO.Params.Is_Quantum == False:
			print("\nThe Inexact Feasible Interior-Point Method starts running ...\n")
		else:
			print("\nThe Quantum Inexact Feasible Interior-Point Method starts running ...\n")
	

	if np.linalg.norm(LO.x) == np.linalg.norm(LO.s) == 0: 			# self-dual embedding standard form LO
		#-----------------------------------------------------------------------
		# Building embeding self-dual formulation
		#-----------------------------------------------------------------------
		LO.x 		= np.ones(LO.n)  				# primal solution
		LO.y 		= np.ones(LO.m)					# dual solution
		LO.s 		= np.ones(LO.n) 	 			# dual slack
		xx 		= np.ones(LO.n)  				# primal solution
		yy 		= np.ones(LO.m)					# dual solution
		ss 		= np.ones(LO.n) 	 			# dual slack
		tau			= 1
		u 			= 1
		theta		= 1

		bb			= LO.b - np.dot(LO.A, np.ones(LO.n))
		cc 			= LO.c - np.ones(LO.n) - np.dot(LO.A.T, np.ones(LO.m))
		oo			= 1 + np.dot(LO.c.T, np.ones(LO.n)) - np.dot(LO.b.T, np.ones(LO.m))

		A1			= np.concatenate((np.concatenate( (np.concatenate(([1], np.zeros(LO.n)), axis=0),np.concatenate((LO.c.T,[0]), axis=0)), axis=0),np.concatenate((-1*LO.b.T,[-oo]), axis=0)), axis=0)
		A2			= np.concatenate((np.concatenate( (np.concatenate((np.zeros((LO.n,1)), np.identity(LO.n)), axis=1),np.concatenate((np.zeros((LO.n, LO.n)),np.matrix(-1*LO.c).T), axis=1)), axis=1),np.concatenate((LO.A.T,np.matrix(cc).T), axis=1)), axis=1)
		A3			= np.concatenate((np.concatenate( (np.concatenate((np.zeros((LO.m,1)), np.zeros((LO.m, LO.n))), axis=1),np.concatenate((LO.A,np.matrix(-1*LO.b).T), axis=1)), axis=1),np.concatenate((np.zeros((LO.m, LO.m)),np.matrix(bb).T), axis=1)), axis=1)
		A4			= np.concatenate((np.concatenate( (np.concatenate((np.zeros((1,1)), np.zeros((1,LO.n))), axis=1),np.concatenate((np.matrix(cc),np.matrix([-1*oo])), axis=1)), axis=1),np.concatenate((np.matrix(-1*bb),np.zeros((1,1))), axis=1)), axis=1)
		
		AA 			= np.concatenate((np.concatenate((np.matrix(A1),A2), axis=0),np.concatenate((A3,A4), axis=0)), axis=0)
		Acopy		= deepcopy(AA)
		
		#-----------------------------------------------------------------------
		# finding basis for null space
		#-----------------------------------------------------------------------
		for i in range(LO.n+1, LO.n+LO.m+2):
			AA[i,]	= AA[i,]/AA[i,i]
			for j in range(LO.n+LO.m+2):
				if j !=i:
					AA[j,]	= AA[j,] - AA[j,i] * AA[i,]
		
		V 		= np.concatenate((AA[:,LO.n+LO.m+2:],-1*np.identity(LO.n+1)),axis=0)
		lam 	= np.zeros(LO.n+1)

		#-----------------------------------------------------------------------
		# Loop with 
		#-----------------------------------------------------------------------
		while (True):

			#-------------------------------------------------------------------
			# Calculate complementarity 
			#-------------------------------------------------------------------
			mu 		= (np.dot(xx.T,ss) +tau*u )/ (LO.n+1)


			#-------------------------------------------------------------------
			# Calculate the intermediate terms
			#-------------------------------------------------------------------
			X		= np.diag(xx)
			S 		= np.diag(ss)
		
			D1 		= np.concatenate((np.concatenate(([tau],np.zeros(2*LO.n)),axis=0),np.concatenate(([u],np.zeros(LO.m+1)),axis=0)),axis=0)
			D2		= np.concatenate((np.concatenate((np.zeros((LO.n,1)),X),axis=1),np.concatenate((S,np.zeros((LO.n,LO.m+2))),axis=1)),axis=1)
			
			D       = np.concatenate((np.matrix(D1),D2),axis=0)		
									

			M 		= np.dot(D,V) 	

			
			r 		= (beta_1*mu)*np.ones(LO.n+1)-np.concatenate(([tau*u],np.dot(X,ss)),axis=0)	
							
			#-------------------------------------------------------------------
			# Liner system solver solution
			#-------------------------------------------------------------------
			lam, norm_of_residual, is_sign_changed 	= linear_system_solver(M, r, LO.Params) 

			if LO.Params.Is_Quantum == False:
					
				# add random error to solution	
				eps 				= (eta*mu)/np.linalg.norm(M)
				err                 = eps/len(lam)
				#for i in range(len(lam)):
				#	lam[i]	= lam[i] + (-err+2*np.random.rand()*err)
				

			#-------------------------------------------------------------------
			# Calculate delta_x and delta_s
			#-------------------------------------------------------------------
			
			delta_X		= np.dot(V, lam).T
			
			delta_u 	= delta_X[0:1]
			delta_s 	= delta_X[1:LO.n+1]
			delta_x 	= delta_X[LO.n+1:2*LO.n+1]
			delta_tau 	= delta_X[2*LO.n+1:2*LO.n+2]
			delta_y 	= delta_X[2*LO.n+2:2*LO.n+2+LO.m]
			delta_theta = delta_X[2*LO.n+2+LO.m:2*LO.n+3+LO.m]
			
			u 			= u + delta_u[0,0]
			ss 			= ss + np.squeeze(np.asarray(delta_s)) 
			xx			= xx + np.squeeze(np.asarray(delta_x))
			tau 		= tau + delta_tau[0,0]
			yy 			= yy + np.squeeze(np.asarray(delta_y)) 
			theta 		= theta + delta_theta[0,0]
			

			if tau>10^(-10):
				LO.x=xx/tau
				LO.s=ss/tau
				LO.y=yy/tau

			#-------------------------------------------------------------------
			# Calculate the run-time
			#-------------------------------------------------------------------
			end_time 	= time()
			run_time 	= end_time - start_time

			#-------------------------------------------------------------------
			# Print the results
			#-------------------------------------------------------------------
			print_iteration_IPM(LO, M, r, iteration, norm_of_residual, is_sign_changed, 1, 1, run_time)

			iteration 	+= 1

			#-------------------------------------------------------------------
			# Stopping condition
			#-------------------------------------------------------------------
			if (LO.complementarity()  <= precision  ): break
		
		if tau<=10^(-10):
			sys.exit("The problem is infeasible.")
		else:
			LO.x=xx/tau
			LO.s=ss/tau
			LO.y=yy/tau

		print_final_IPM(LO, iteration, run_time)
		

		return LO
	else:				# standard LO with interior feasible solution

		AA = deepcopy(LO.A)
		#-----------------------------------------------------------------------
		# finding basis for null space
		#-----------------------------------------------------------------------
		for i in range(LO.m):
			AA[i,]=AA[i,]/AA[i,i]
			for j in range(LO.m):
				if j!=i:
					AA[j,]=AA[j,]-AA[j,i]*AA[i,]
		V 	= np.concatenate((AA[:,LO.m:],-1*np.identity(LO.n-LO.m)),axis=0)

		
		#-----------------------------------------------------------------------
		# Loop with 
		#-----------------------------------------------------------------------
		while (True):

			#-------------------------------------------------------------------
			# Calculate complementarity 
			#-------------------------------------------------------------------
			mu 		= np.dot(LO.x, LO.s) / (LO.n)


			#-------------------------------------------------------------------
			# Calculate the intermediate terms
			#-------------------------------------------------------------------
			X		= np.diag(LO.x)
			S 		= np.diag(LO.s)


			M 		= np.concatenate((np.matmul(-X, LO.A.T), np.matmul(S,V)),axis=1) 	
			r 		= (beta_1*mu)*np.ones(LO.n) - np.dot(X,LO.s)
								
			#-------------------------------------------------------------------
			# Liner system solver solution
			#-------------------------------------------------------------------
			lam, norm_of_residual, is_sign_changed 	= linear_system_solver(M, r, LO.Params) 
			
			if LO.Params.Is_Quantum == False:
				# # add random error to solution	
				eps 				= (eta*mu)/np.linalg.norm(M)
				err                 = eps/len(z)
				z  					= z + err * (np.ones(len(z)) - 2 * np.random.rand(len(z)) )

						

			#-------------------------------------------------------------------
			# Calculate delta_x and delta_s
			#-------------------------------------------------------------------
			
			delta_y		= z[:LO.m]
			lam			= z[LO.m:]
			delta_x 	= np.dot(V,lam)
			delta_s 	= - np.dot(LO.A.T, delta_y)

			
			LO.s 		= LO.s + delta_s
			LO.x 		= LO.x + delta_x
			LO.y 		= LO.y + delta_y

			

			#-------------------------------------------------------------------
			# Calculate the run-time
			#-------------------------------------------------------------------
			end_time 	= time()
			run_time 	= end_time - start_time

			#-------------------------------------------------------------------
			# Print the results
			#-------------------------------------------------------------------
			print_iteration_IPM(LO, M, r, iteration , norm_of_residual, is_sign_changed, 1, 1, run_time)

			iteration 	+= 1

			#-------------------------------------------------------------------
			# Stopping condition
			#-------------------------------------------------------------------
			if (LO.complementarity()  <= precision  ): break
		

		print_final_IPM(LO, iteration, run_time)
		
		return LO




