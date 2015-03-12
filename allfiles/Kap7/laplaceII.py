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
from numpy import cosh, cos

# Geometry
width = 1.0
height = width

Nx = 20 # number of points in x-direction
Ny = Nx # number of points in y-direction
h = width/Nx
N = Nx*Ny

h=width/Nx

def laplace_direct_xdir(N):
    diagonals=np.zeros((5,N))
    
    diagonals[0,:]=1.0   # -Nx diagonal, the last Nx-elt are not used         
 
    diagonals[1,:]=1.0          
    diagonals[1,Nx-1::Nx]=0.0  # Left dT/dy           
    
    diagonals[2,:]=-4.0         
    
    diagonals[3,:]=1.0
    diagonals[3,1::Nx]=2.0           # every 4th elt set to 2,staring at idx           
    diagonals[3,Nx::Nx]=0.0  # Left dT/dy           
    
    diagonals[4,:]=1.0       # Nx-diagonal,the frst Nx-elts are not used
    diagonals[4,Nx:2*Nx]=2.0 

    
    A=sc.sparse.spdiags(diagonals, [-Nx,-1,0,1,Nx], N, N,format='csc') #sparse matrix instance
#     print A.todense()
    
    d=np.zeros(N)
    d[-Nx:]=-1.0 # set the last Nx elts to -1.0
    
#     print d
    Tvector = sc.sparse.linalg.spsolve(A,d) #theta=sc.linalg.solve_triangular(A,d)
    
    Tmatrix=np.reshape(Tvector, (Ny, Nx), order='C')
    
    T = np.zeros((Ny+1,Nx+1))
    T[-1,:]=1.0 # T=1 for y=1
    
    T[0:-1,0:-1] = Tmatrix
    
    return T



def T_analytical(x,y):
    N_inf=220
    
    #T = np.zeros((len(x),len(y)))
    T = np.zeros_like(x)
    
    for n in range(1,N_inf+1):
        lambda_n = (2*n-1)*np.pi/2.0
        An = 2*(-1)**(n-1)/(lambda_n*cosh(lambda_n))
        T+=An*cosh(lambda_n*y)*cos(lambda_n*x)
        
    return T

    
def T2(x,y):
    N_inf=30
    T = np.zeros_like(x)
    n=1.0
    for n in range(1,N_inf+1):
        lambda_n = (2*n-1)*np.pi/2.0
        An = 2*(-1)**(n-1)/(lambda_n*cosh(lambda_n))
        T+=An*np.cosh(lambda_n*x)*np.cos(lambda_n*y)
    
    return T


def plot2D(x, y, p):
    fig = plt.figure(figsize=(11,7), dpi=100)
    ax = fig.gca(projection='3d')
    X,Y = np.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,p[:], rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    ax.view_init(30,225)
    plt.xlabel('x-values')
    plt.ylabel('y-values')

T = np.zeros((Ny+1,Nx+1))
T=laplace_direct_xdir(N)

x = np.linspace(0, width, Nx+1)
y = np.linspace(0, height, Ny+1)
X,Y = np.meshgrid(x, y)

Ta=T_analytical(X,Y)

plot2D(x,y,Ta)
plt.title('Analytical solution')

plot2D(x,y,T)
plt.title('numerical')


plt.show()
plt.close()

