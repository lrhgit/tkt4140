# coding: utf-8
# src-ch7/laplace_Diriclhet2.py
import numpy as np
import scipy 
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pylab as plt
import time
from math import sinh
from astropy.units import dT

#import matplotlib.pyplot as plt

# Change some default values to make plots more readable on the screen
LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT

# Set temperature at the top
Ttop=10
Tbottom=0.0
Tleft=0.0
Tright=0.0

xmax=1.0
ymax=1.5
 
# Set simulation parameters
#need hx=(1/nx)=hy=(1.5/ny)
Nx = 20
h=xmax/Nx
Ny = int(ymax/h)

nx = Nx-1
ny = Ny-1
n = (nx)*(ny) #number of unknowns

# surfaceplot:
x = np.linspace(0, xmax, Nx + 1)
y = np.linspace(0, ymax, Ny + 1)

X, Y = np.meshgrid(x, y)

T = np.zeros_like(X)

# set the imposed boudary values
T[-1,:] = Ttop
T[0,:] = Tbottom
T[:,0] = Tleft
T[:,-1] = Tright

T2 = T.copy()

reltol=1.0e-3

omega = 1.5
iteration = 0
rel_res=1.0

# Gauss-Seidel iterative solution
tic=time.clock()
while (rel_res > reltol):
    dTmax=0.0
    for j in range(1,ny+1):
        for i in range(1, nx + 1):
            R = (T[j,i-1]+T[j-1,i]+T[j,i+1]+T[j+1,i]-4.0*T[j,i])
            dT = 0.25*omega*R
            T[j,i]+=dT
            dTmax=np.max([np.abs(dT),dTmax])
        
    rel_res=dTmax/np.max(np.abs(T))
    iteration+=1
    
toc=time.clock()
print "Gauss-Seidel solver time:\t{0:0.2f}. \t Nr. iterations {1}".format(toc-tic,iteration)


iteration = 0
rel_res=1.0
# Jacobi iterative solution
tic=time.clock()
while (rel_res > reltol):
    R2 = (T2[1:-1,0:-2]+T2[0:-2,1:-1]+T2[1:-1,2:]+T2[2:,1:-1]-4.0*T2[1:-1,1:-1])
    dT2 = 0.25*R2
    T2[1:-1,1:-1]+=dT2
    rel_res=np.max(dT2)/np.max(T2)
    iteration+=1    

toc=time.clock()
print "Jacobi solver time:\t\t{0:0.2f}. \t Nr. iterations {1}".format(toc-tic,iteration)

# Comment out below to plot Gauss-Seidel solution     
# T = T2

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, T, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_zlim(0, Ttop+10)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('T [$^o$C]')


nx=4
xticks=np.linspace(0.0,xmax,nx+1)
ax.set_xticks(xticks)

ny=8
yticks=np.linspace(0.0,ymax,ny+1)
ax.set_yticks(yticks)

nTicks=5
dT=int(Ttop/nTicks)
Tticklist=range(0,Ttop+1,dT)
ax.set_zticks(Tticklist)

#fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
