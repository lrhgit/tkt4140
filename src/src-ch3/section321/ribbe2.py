# src-ch3/section321/ribbe2.py

import numpy as np
import scipy as sc
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import time
from numpy import cosh

#import matplotlib.pyplot as plt
from matplotlib.pyplot import *
# Change some default values to make plots more readable on the screen
LNWDT=3; FNT=20
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT

# Set simulation parameters
beta = 3.0
h = 0.001               # element size
L =1.0               # length of domain
n = int(round(L/h))  # # of unknowns, assuming known bndry values at outlet
x=np.arange(n+1)*h      # x includes min and max at boundaries were bc are imposed.


#Define useful functions

def tri_diag_setup(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def theta_analytical(beta,x):
    return np.cosh(beta*x)/np.cosh(beta)

#Create matrix for linalg solver
a=np.ones(n-1)                  # sub-diagonal
b=-np.ones(n)*(2+(beta*h)**2)   # diagonal
c=np.ones(n-1)                  # sub-diagonal
#c=a.copy()                      # super-diagl, copy as elts are modified later
#c=a
# particular diagonal values due to derivative bc
version=2
if (version==1):
    c[0]=2.0
else:      
    b[0]=2.0
    c[0]=-(2-(beta*h)**2)
    print 'version 2' 
                     
A=tri_diag_setup(a,b,c)

#Create matrix for sparse solver
diagonals=np.zeros((3,n))
diagonals[0,:]= 1.0               # all elts in first row is set to 1
diagonals[0,0]= 1.0               
diagonals[1,:]= -(2+(beta*h)**2)  
diagonals[2,:]=1.0

if (version==1):
    diagonals[2,1]= 2.0               # index 1 as the superdiagonal of spdiags is not used,  
else:
    diagonals[1,0]=2.0                # Sets the first element in the main diagonal  
    diagonals[2,1]= -(2+(beta*h)**2)  # index 1 as the superdiagonal of spdiags is not used,  
 super-diagonal

A_sparse = sc.sparse.spdiags(diagonals, [-1,0,1], n, n,format='csc') #sparse matrix instance

#Crete rhs array
d=np.zeros(n)
d[-1]=-1

#Solve linear problems
tic=time.clock()
theta = sc.sparse.linalg.spsolve(A_sparse,d) #theta=sc.linalg.solve_triangular(A,d)
toc=time.clock()
print 'sparse solver time:',toc-tic

tic=time.clock()
theta2=sc.linalg.solve(A,d,)
toc=time.clock()
print 'linalg solver time:',toc-tic

# Plot solutions
plot(x[0:-1],theta,x[0:-1],theta2,'-.',x,theta_analytical(beta,x),':')
xlabel('x')
ylabel(r'Dimensionless temperature $\mathregular{\theta}$')
legend(['sparse','linalg','analytical'])
show()
close()
print 'done'