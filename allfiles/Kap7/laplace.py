'''
Created on Mar 5, 2015

@author: leifh
'''
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import scipy as sc
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import time


# Geometry
width = 1.0
height = width

Nx = 4 # number of points in x-direction
Ny = Nx # number of points in y-direction
h = width/Nx
N = Nx*Ny

h=width/Nx

diagonals=np.zeros((5,N))
diagonals[0,:]=1.0               # all elts in first row is set to 1

diagonals[1,:]=1.0               # all elts in first row is set to 1
diagonals[1,Nx-1::Nx]=0.0            # every 4th elt set to 2,staring at idx 

diagonals[2,:]=-4.0               # all elts in first row is set to 1

diagonals[3,:]=1.0               # all elts in first row is set to 1
diagonals[3,1::Nx]=2.0           # every 4th elt set to 2,staring at idx 
diagonals[3,Nx::Nx]=0.0            # every 4th elt set to 2,staring at idx 3

diagonals[4,:]=1.0               # all elts in first row is set to 1n
diagonals[4,0:Nx-1]=2.0               # all elts in first row is set to 1

A=sc.sparse.spdiags(diagonals, [-Nx,-1,0,1,Nx], N, N,format='csc') #sparse matrix instance
#print A.todense()

d=np.zeros(N)
d[-Nx:]=-1.0 # set the last Nx elts to -1.0
Tvector = sc.sparse.linalg.spsolve(A,d) #theta=sc.linalg.solve_triangular(A,d)

Tmatrix=np.reshape(Tvector, (Nx, Ny))

T = np.zeros((Nx+1,Ny+1))
T[-1,:]=1.0

T[0:-1,0:-1] = Tmatrix

x = np.linspace(0, width, Nx+1)
y = np.linspace(0, height, Ny+1)
X,Y = np.meshgrid(x, y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, T, rstride=1, cstride=1, cmap=cm.coolwarm)
plt.xlabel('x-values')
plt.ylabel('y-values')
plt.show()
plt.close()






