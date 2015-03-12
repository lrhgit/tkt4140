#
# Solve y'' = 3/2 y**2
# y(0) = 4 and y(1) = 1
#
# Implicit tri-diagonal solution
# Linearized with Newton-linearization
#

import numpy as np
import scipy as sc
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import time
from math import sinh

#import matplotlib.pyplot as plt
from matplotlib.pyplot import *
# Change some default values to make plots more readable on the screen
LNWDT=3; FNT=20
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT

# Set simulation parameters
h=0.1               # element size
L=1.0                  # length of domain
n=int(round(L/h)) -1  # number of unknowns, assuming known boundary values
x=np.arange(n+2)*h      # x includes min and max at boundaries were bc are imposed.

def ya(x):
     return 4.0/(1 + x)**2;

#ym=np.zeros(n)
ym=-20.0*(x[1:-1]-x[1:-1]**2)
#ym=-20.0*(x[1:-1])
#ym=-20.0*np.ones_like(x[1:-1])

diagonals=np.zeros((3,n))
d=np.zeros(n)

#Create matrix for sparse solver
diagonals[0,:]= 1.0                       #all elts in first row is set to 1
diagonals[2,:]= 1.0 

Nmax = 10
TOL = 1.0E-10

print 'Iter. \t Error'
for i in range(Nmax-1):
    # compute main diagonal of eqn system
    diagonals[1,:]= -(2.0+3.0*h**2*ym)  
    A = sc.sparse.spdiags(diagonals, [-1,0,1], n, n,format='csc') #sparse matrix instance
    
    # Crete rhs and impose bounary conditions
    d[:]= -3.0/2.0*h**2*ym**2
    d[0]+=-4.0
    d[-1]+=-1.0
    
    # Solve the linearized system
    y=sc.sparse.linalg.spsolve(A,d) 

    dy_max = np.max(np.abs((y-ym)/y))
    
    if (dy_max<TOL): 
        print 'Solution converged after',i, 'iterations.'
        break
    else:
        print '{0:3d} {1:10.7f}'.format(i, dy_max) 
    ym=y

# Make solution vector with boundary conditions
yn = np.zeros_like(x)
yn[0]=ya(x[0]); yn[-1]=ya(x[-1]); yn[1:-1]=y 

plot(x,yn,'-o',x,ya(x),':')    
legend(['numerical','analytical'])
show()
close()
