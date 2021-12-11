from ipm import *
import numpy as np
from numpy import linalg as LA
import sys


parameters 	= Parameters()

#parameters.seed  		= int(sys.argv[1])
#parameters.norm_rhs 	= 2
#parameters.cond 		= float(sys.argv[2])

#A, b, c, int_x, int_y, int_s	= generate_problem(m=3, n=4, parameters=parameters)
m=1
n=2
A=np.zeros((m,n,n))
for i in range(m):
  R=1*np.random.rand(n,n)
  A[i]=np.matmul(R.T,R)
A[0]=np.eye(n)
R=1*np.random.rand(n,n)
x=np.matmul(R.T,R)
v,Q=LA.eigh(x)
for i in range(len(v)):
  if v[i]==0:
    v[i]=1
x=np.dot(np.dot(Q,np.diag(v)),Q.T)
R=1*np.random.rand(n,n)
s=np.matmul(R.T,R)
v,Q=LA.eigh(s)
for i in range(len(v)):
  if v[i]==0:
    v[i]=1
s=np.dot(np.dot(Q,np.diag(v)),Q.T)
y=np.random.rand(m)
b=np.zeros(m)
c=s
for i in range(m):
  b[i]=np.trace(np.matmul(A[i],x))
  c=c+y[i]*A[i]

model 		= SDO(A, b, c)
model.x 	= x
model.y 	= y
model.s 	= s
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