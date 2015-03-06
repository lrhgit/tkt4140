'''
Created on Mar 5, 2015

@author: leifh
'''
import numpy as np
import scipy as sc
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import time


# Geometry
width = 1.0
height = width

Nx = 30 # number of points in x-direction
Ny = Nx # number of points in y-direction
h = width/Nx
N = Nx*(Nx+1)/2#Symmetric problem

h=width/Nx

diagonals=np.zeros((5,N))
diagonals[0,:]=1.0               # all elts in first row is set to 1
diagonals[1,:]=1.0               # all elts in first row is set to 1
diagonals[2,:]=-4.0               # all elts in first row is set to 1

diagonals[3,:]=1.0               # all elts in first row is set to 1
diagonals[3,0::4]=2.0            # every 4th elt set to 2,staring at idx 
diagonals[3,3::4]=0.0            # every 4th elt set to 2,staring at idx 3

diagonals[4,:]=1.0               # all elts in first row is set to 1n
diagonals[4,0:3]=2.0               # all elts in first row is set to 1

A=sc.sparse.spdiags(diagonals, [-4,-1,0,1,4], N, N,format='csc') #sparse matrix instance

d=np.zeros(N)
d[-4:]=-1.0 # set the last 4 elts to -1.0
T = sc.sparse.linalg.spsolve(A,d) #theta=sc.linalg.solve_triangular(A,d)

print 'done'
print T



