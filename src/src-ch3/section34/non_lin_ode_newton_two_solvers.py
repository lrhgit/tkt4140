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
h = 0.1               # element size
L =1.0                  # length of domain
n = int(round(L/h)) -1  # number of unknowns, assuming known boundary values
x=np.arange(n+2)*h      # x includes min and max at boundaries were bc are imposed.


def tri_diag_setup(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def ya(x):
     return 4.0/(1 + x)**2;

#ym=3.0/(2-x[1:-1]);
ym=np.zeros(n)
ym_linalg=np.zeros(n)
diagonals=np.zeros((3,n))
d=np.zeros(n)
diagonals[0,:]= 1.0                       #all elts in first row is set to 1
diagonals[2,:]= 1.0 
#Create matrix for sparse solver

Nmax = 50
TOL = 1.0E-5

# for linalg
a=np.ones(n-1)
b=np.zeros(n)
c=np.ones(n-1)

for i in range(Nmax-1):
    # compute main diagonal of eqn system
    diagonals[1,:]= -(2.0+3.0*h**2*ym)  
    A = sc.sparse.spdiags(diagonals, [-1,0,1], n, n,format='csc') #sparse matrix instance
    
    b=-(2.0+3.0*h**2*ym_linalg)  
    A_linalg=tri_diag_setup(a,b,c)
                                                                        
    # Crete rhs and impose bounary conditions
    d[:]= -3.0/2.0*h**2*ym**2
    d[0]=d[0]-4.0
    d[-1]=d[-1]-1.0
    
    # Solve the linearized system
    y=sc.sparse.linalg.spsolve(A,d) 
    
    y_linalg=sc.linalg.solve(A_linalg,d,)
    
    dy_max = np.max(np.abs((y-ym)))
    if (dy_max<TOL): 
        print 'solution converged after ', i, ' iterations.'
        break

    ym=y
    ym_linalg = y_linalg

Adense=sc.sparse.spdiags(diagonals, [-1,0,1], n, n,format='csc').todense()

plot(x[1:-1],y,x[1:-1],y_linalg,x,ya(x),':')    
legend(['sparse','linalg','analytical'])
show()
close()
