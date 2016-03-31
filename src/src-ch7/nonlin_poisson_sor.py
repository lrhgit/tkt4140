# coding: utf-8
# src-ch7/laplace_Diriclhet2.py; Visualization.py @ git@lrhgit/tkt4140/src/src-ch7/Visualization.py;
import numpy as np
import scipy 
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pylab as plt
import time
import gauss_seidel_sor as gs


#import matplotlib.pyplot as plt

# Change some default values to make plots more readable on the screen
LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT

xmax=1.0
ymax=1.0
 
# Set simulation parameters
#need hx=(1/nx)=hy=(1.5/ny)
Nx = 40
h=xmax/Nx
Ny = int(ymax/h)

nx = Nx-1
ny = Ny-1
n = (nx)*(ny) #number of unknowns

# surfaceplot:
x = np.linspace(0, xmax, Nx + 1)
y = np.linspace(0, ymax, Ny + 1)

X, Y = np.meshgrid(x, y)

U = np.zeros_like(X)
U2 = U.copy()
U3 = U.copy()

reltol=1.0e-3
 
rho=np.cos(np.pi/nx)
 
omega = 2/(1+np.sqrt(1-rho**2))
print 'omega={0:.2f}'.format(omega)
 
iteration = 0
rel_res=1.0

# First method 
# Python Gauss-Seidel method
tic=time.clock()
while (rel_res > reltol):
    du_max=0.0
    for j in range(1,ny+1):
        for i in range(1, nx + 1):
            R = (U[j,i-1]+U[j-1,i]+U[j,i+1]+U[j+1,i]-4.0*U[j,i]) + h**2*(U[j,i]**2+1.0)
            df=4-2*h**2*U[j,i]
            dU = omega*R/df
            U[j,i]+=dU
            du_max=np.max([np.abs(dU),du_max])

         
    rel_res=du_max/np.max(np.abs(U))

    iteration+=1
     
toc=time.clock()
print "Python Gauss-Seidel CPU-time:\t{0:0.2f}. \t Nr. iterations {1}".format(toc-tic,iteration)
 
 
iteration = 0
rel_res=1.0

# Second method
# Jacobi iterative solution
tic=time.clock()
while (rel_res > reltol):
    R2 = (U2[1:-1,0:-2]+U2[0:-2,1:-1]+U2[1:-1,2:]+U2[2:,1:-1]-4.0*U2[1:-1,1:-1]) + h**2*(U2[1:-1,1:-1]**2+1.0)
    df=4-2*h**2*U2[1:-1,1:-1]
    dU2 = R2/df
    U2[1:-1,1:-1]+=dU2
    rel_res=np.max(dU2)/np.max(U2)
    iteration+=1    
 
toc=time.clock()
print "Jacobi CPU-time:\t\t{0:0.2f}. \t Nr. iterations {1}".format(toc-tic,iteration)

# Third method
# Cython Gauss-Seidel method 
rel_res=1.0
tic=time.clock() 
U3, relreturn, itsused=gs.gauss(U3,reltol,h, omega)  
toc=time.clock()
print "Cython Gauss-Seidel CPU-time:\t{0:0.2f}. \t Nr. iterations {1}".format(toc-tic,itsused)

# Plot the solutions 
from Visualization import plot_Surface_yx_3subplots
plot_Surface_yx_3subplots(X,Y,U,U2,U3,['Python Gauss-Seidel','Jacobi and numpy','Cython Gauss-Seidel']) 
 
plt.show()
