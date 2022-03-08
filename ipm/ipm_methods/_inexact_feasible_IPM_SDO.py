from ipm.linear_system_solvers.linear_system_solvers import *
from ipm.print_methods.print_IPM import *

import numpy as np
from time import time
from copy import deepcopy
import sys
def hvec(A):
	n=len(A)
	s=[]
	for i in range(n):
		for j in range(i,n):
			s.append(A[i,j])
	return np.array(s)
def svec(A):
	n=len(A)
	s=[]
	for i in range(n):
		s.append(A[i,i])
		for j in range(i+1,n):
			s.append(np.sqrt(2)*A[i,j])
	return np.array(s)
def smat(s):
  n=int(np.floor(np.sqrt(2*len(s))))
  A=np.zeros((n,n))
  t=0
  for i in range(n):
    A[i,i]=s[t]
    t=t+1
    for j in range(i+1,n):
      A[i,j]=s[t]/np.sqrt(2)
      A[j,i]=s[t]/np.sqrt(2)
      t=t+1
  return A
def symet(M,type):
	if type=="AHO":
		return 0.5*(M+M.T)
	elif type=="NT":
		return M
	else:
		return 0.5*(M+M.T)


#===========================================================================
# Inexact Infeasible Interior Point Method
#===========================================================================
def inexact_feasible_IPM_SDO(SDO):
	#-----------------------------------------------------------------------
	# Initialize the paramters of the algorithm
	#-----------------------------------------------------------------------
	beta_1 			= 0.5	#SDO.Params.Beta_1					# Default value: 0.1
	eta 			= 0.5
	precision 		= SDO.Params.LO_Precision

	iteration 		= 0
	
	start_time 		= time()
	nn=int(SDO.n*(SDO.n+1)/2)

	if SDO.Params.LO_Verbosity > -1:
		if SDO.Params.Is_Quantum == False:
			print("\nThe Inexact Feasible Interior-Point Method starts running ...\n")
		else:
			print("\nThe Quantum Inexact Feasible Interior-Point Method starts running ...\n")
	

	if np.linalg.norm(SDO.x) == np.linalg.norm(SDO.s) == 0: 			# self-dual embedding standard form SDO
		#-----------------------------------------------------------------------
		# Building embeding self-dual formulation
		#-----------------------------------------------------------------------
		SDO.x 		= np.ones(SDO.n)  				# primal solution
		SDO.y 		= np.ones(SDO.m)					# dual solution
		SDO.s 		= np.ones(SDO.n) 	 			# dual slack
		xx 		= np.ones(SDO.n)  				# primal solution
		yy 		= np.ones(SDO.m)					# dual solution
		ss 		= np.ones(SDO.n) 	 			# dual slack
		tau			= 1
		u 			= 1
		theta		= 1

		bb			= SDO.b - np.dot(SDO.A, np.ones(SDO.n))
		cc 			= SDO.c - np.ones(SDO.n) - np.dot(SDO.A.T, np.ones(SDO.m))
		oo			= 1 + np.dot(SDO.c.T, np.ones(SDO.n)) - np.dot(SDO.b.T, np.ones(SDO.m))

		A1			= np.concatenate((np.concatenate( (np.concatenate(([1], np.zeros(SDO.n)), axis=0),np.concatenate((SDO.c.T,[0]), axis=0)), axis=0),np.concatenate((-1*SDO.b.T,[-oo]), axis=0)), axis=0)
		A2			= np.concatenate((np.concatenate( (np.concatenate((np.zeros((SDO.n,1)), np.identity(SDO.n)), axis=1),np.concatenate((np.zeros((SDO.n, SDO.n)),np.matrix(-1*SDO.c).T), axis=1)), axis=1),np.concatenate((SDO.A.T,np.matrix(cc).T), axis=1)), axis=1)
		A3			= np.concatenate((np.concatenate( (np.concatenate((np.zeros((SDO.m,1)), np.zeros((SDO.m, SDO.n))), axis=1),np.concatenate((SDO.A,np.matrix(-1*SDO.b).T), axis=1)), axis=1),np.concatenate((np.zeros((SDO.m, SDO.m)),np.matrix(bb).T), axis=1)), axis=1)
		A4			= np.concatenate((np.concatenate( (np.concatenate((np.zeros((1,1)), np.zeros((1,SDO.n))), axis=1),np.concatenate((np.matrix(cc),np.matrix([-1*oo])), axis=1)), axis=1),np.concatenate((np.matrix(-1*bb),np.zeros((1,1))), axis=1)), axis=1)
		
		AA 			= np.concatenate((np.concatenate((np.matrix(A1),A2), axis=0),np.concatenate((A3,A4), axis=0)), axis=0)
		Acopy		= deepcopy(AA)
		
		#-----------------------------------------------------------------------
		# finding basis for null space
		#-----------------------------------------------------------------------
		for i in range(SDO.n+1, SDO.n+SDO.m+2):
			AA[i,]	= AA[i,]/AA[i,i]
			for j in range(SDO.n+SDO.m+2):
				if j !=i:
					AA[j,]	= AA[j,] - AA[j,i] * AA[i,]
		
		V 		= np.concatenate((AA[:,SDO.n+SDO.m+2:],-1*np.identity(SDO.n+1)),axis=0)
		lam 	= np.zeros(SDO.n+1)

		#-----------------------------------------------------------------------
		# SDOop with 
		#-----------------------------------------------------------------------
		while (True):

			#-------------------------------------------------------------------
			# Calculate complementarity 
			#-------------------------------------------------------------------
			mu 		= (np.dot(xx.T,ss) +tau*u )/ (SDO.n+1)


			#-------------------------------------------------------------------
			# Calculate the intermediate terms
			#-------------------------------------------------------------------
			X		= np.diag(xx)
			S 		= np.diag(ss)
		
			D1 		= np.concatenate((np.concatenate(([tau],np.zeros(2*SDO.n)),axis=0),np.concatenate(([u],np.zeros(SDO.m+1)),axis=0)),axis=0)
			D2		= np.concatenate((np.concatenate((np.zeros((SDO.n,1)),X),axis=1),np.concatenate((S,np.zeros((SDO.n,SDO.m+2))),axis=1)),axis=1)
			
			D       = np.concatenate((np.matrix(D1),D2),axis=0)		
									

			M 		= np.dot(D,V) 	

			
			r 		= (beta_1*mu)*np.ones(SDO.n+1)-np.concatenate(([tau*u],np.dot(X,ss)),axis=0)	
							
			#-------------------------------------------------------------------
			# Liner system solver solution
			#-------------------------------------------------------------------
			lam, norm_of_residual, is_sign_changed 	= linear_system_solver(M, r, SDO.Params) 
			if SDO.Params.Is_Quantum == False:
				
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
			delta_s 	= delta_X[1:SDO.n+1]
			delta_x 	= delta_X[SDO.n+1:2*SDO.n+1]
			delta_tau 	= delta_X[2*SDO.n+1:2*SDO.n+2]
			delta_y 	= delta_X[2*SDO.n+2:2*SDO.n+2+SDO.m]
			delta_theta = delta_X[2*SDO.n+2+SDO.m:2*SDO.n+3+SDO.m]
			
			u 			= u + delta_u[0,0]
			ss 			= ss + np.squeeze(np.asarray(delta_s)) 
			xx			= xx + np.squeeze(np.asarray(delta_x))
			tau 		= tau + delta_tau[0,0]
			yy 			= yy + np.squeeze(np.asarray(delta_y)) 
			theta 		= theta + delta_theta[0,0]
			

			if tau>10^(-10):
				SDO.x=xx/tau
				SDO.s=ss/tau
				SDO.y=yy/tau

			#-------------------------------------------------------------------
			# Calculate the run-time
			#-------------------------------------------------------------------
			end_time 	= time()
			run_time 	= end_time - start_time

			#-------------------------------------------------------------------
			# Print the results
			#-------------------------------------------------------------------
			print_iteration_IPM(SDO, M, r, iteration, norm_of_residual, is_sign_changed, 1, 1, run_time)

			iteration 	+= 1

			#-------------------------------------------------------------------
			# Stopping condition
			#-------------------------------------------------------------------
			if (SDO.complementarity()  <= precision  ): break
		
		if tau<=10^(-10):
			sys.exit("The problem is infeasible.")
		else:
			SDO.x=xx/tau
			SDO.s=ss/tau
			SDO.y=yy/tau

		print_final_IPM(SDO, iteration, run_time)
		

		return SDO
	else:				# standard SDO with interior feasible solution



		
		AA = np.zeros((SDO.m,nn))
		for i in range(SDO.m):
			t=svec(SDO.A[i])
			for j in range(len(t)):
				AA[i,j]=t[j]

		#-----------------------------------------------------------------------
		# finding basis for null space
		#-----------------------------------------------------------------------
		for i in range(SDO.m):
			AA[i,]=AA[i,]/AA[i,i]
			for j in range(SDO.m):
				if j!=i:
					AA[j,]=AA[j,]-AA[j,i]*AA[i,]
		v 	= np.concatenate((AA[:,SDO.m:],-1*np.identity(nn-SDO.m)),axis=0)
		V=np.zeros((nn-SDO.m,SDO.n,SDO.n))
		for i in range(nn-SDO.m):
			V[i]=smat(v[:,i])



		
		#-----------------------------------------------------------------------
		# SDOop with 
		#-----------------------------------------------------------------------
		while (True):

			#-------------------------------------------------------------------
			# Calculate complementarity 
			#-------------------------------------------------------------------
			mu 		= SDO.complementarity() / (SDO.n)


			#-------------------------------------------------------------------
			# Calculate the intermediate terms
			#-------------------------------------------------------------------
			X		= svec(SDO.x)
			S 		= svec(SDO.s)

			HVS=np.zeros((nn-SDO.m,SDO.n,SDO.n))
			for i in range(nn-SDO.m):
				HVS[i]=symet(np.matmul(V[i],SDO.s),"AHO")
			HXA=np.zeros((SDO.m,SDO.n,SDO.n))
			for i in range(SDO.m):
				HXA[i]=symet(np.matmul(SDO.x,SDO.A[i]),"AHO")
			M=np.zeros((nn,nn))
			for i in range(nn-SDO.m):
					M[:,i]=hvec(HVS[i])
			for i in range(SDO.m):
					M[:,i+nn-SDO.m]=-1*hvec(HXA[i])
	
			r 		= hvec((beta_1*mu)*np.eye(SDO.n) - symet(np.matmul(SDO.x,SDO.s),"AHO"))
			
			MTM=np.matmul(M.T,M)
			MTr=np.dot(M.T,r)					
			#-------------------------------------------------------------------
			# Liner system solver solution
			#-------------------------------------------------------------------
			z, norm_of_residual, is_sign_changed  	= linear_system_solver(MTM, MTr, SDO.Params) 
			if SDO.Params.Is_Quantum == False:		
				# # add random error to solution	
				eps 				= (eta*mu)/np.linalg.norm(M)
				err                 = eps/len(z)
				z  					= z + err * (np.ones(len(z)) - 2 * np.random.rand(len(z)) )

						

			#-------------------------------------------------------------------
			# Calculate delta_x and delta_s
			#-------------------------------------------------------------------
			
			delta_y		= z[nn-SDO.m:]
			lam			= z[:nn-SDO.m]
			delta_x=np.zeros((SDO.n,SDO.n))
			for i in range(nn-SDO.m):
				delta_x=delta_x+lam[i]*V[i]
			delta_s=np.zeros((SDO.n,SDO.n))
			for i in range(SDO.m):
				delta_s=delta_s-delta_y[i]*SDO.A[i]

			
			SDO.s 		= SDO.s + delta_s
			SDO.x 		= SDO.x + delta_x
			SDO.y 		= SDO.y + delta_y

			

			#-------------------------------------------------------------------
			# Calculate the run-time
			#-------------------------------------------------------------------
			end_time 	= time()
			run_time 	= end_time - start_time

			#-------------------------------------------------------------------
			# Print the results
			#-------------------------------------------------------------------
			print_iteration_IPM(SDO, M, r, iteration , norm_of_residual, is_sign_changed, 1, 1, run_time)

			iteration 	+= 1

			#-------------------------------------------------------------------
			# Stopping condition
			#-------------------------------------------------------------------
			if (SDO.complementarity()  <= precision  ): break
		

		print_final_IPM(SDO, iteration, run_time)
		
		return SDO




