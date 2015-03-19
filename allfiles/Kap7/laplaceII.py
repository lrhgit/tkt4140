'''
Created on Mar 5, 2015

Solution of Heat equation d2T/dx2 + d2T/dy2 = 0

Robin BCs
          dT/dy = 0 at y=0 
          dT/dx = 0 at x=0 
          T=1 at y=1
          T=0 at x=0
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

Nx = 4 # number of points in x-direction
Ny = Nx # number of points in y-direction
h = width/Nx
N = Nx*Ny

h=width/Nx

def laplace_direct_xdir(N):
    diagonals=np.zeros((5,N))
    
    diagonals[0,:]=1.0         # -Nx diagonal, the last Nx-elt are not used         
 
    diagonals[1,:]=1.0         # 1 diagonal, the last elt is not used                
    diagonals[1,Nx-1::Nx]=0.0  # dT/dx=0 at x=0,every Nx-th elt to 0  
    
    diagonals[2,:]=-4.0        # main diagonal    
    
    diagonals[3,:]=1.0         # 1 diagonal, the last elt is not used                
    diagonals[3,1::Nx]=2.0     # dT/dx=0 at x=0, every Nx-th elt to 2 
    diagonals[3,Nx::Nx]=0.0    # T=0  at x=1, every Nx-th elt 0 
    
    diagonals[4,:]=1.0       
    diagonals[4,Nx:2*Nx]=2.0   # dT/dy=0 at y= 0

    
    A=sc.sparse.spdiags(diagonals, [-Nx,-1,0,1,Nx], N, N,format='csc') #sparse matrix instance
    #print A.todense()
    
    d=np.zeros(N)
    d[-Nx:]=-1.0 # set the last Nx elts to -1.0
    
#     print d
    Tvector = sc.sparse.linalg.spsolve(A,d) #theta=sc.linalg.solve_triangular(A,d)
    
    Tmatrix=np.reshape(Tvector, (Ny, Nx), order='C')
    
    T = np.zeros((Ny+1,Nx+1))
    T[-1,:]=1.0 # T=1 for y=1
    
    T[0:-1,0:-1] = Tmatrix
    
    return T

def laplace2d(T, dx, dy, l1_eps):
    l1norm = 1.0
    Tn = np.empty_like(T)
    
    while l1norm > l1_eps:
        Tn = T.copy()

        T[1:-1,1:-1] = (dx**2*(Tn[2:,1:-1]+Tn[0:-2,1:-1])+dy**2*(Tn[1:-1,2:]+Tn[1:-1,0:-2]))/(2*(dx**2+dy**2)) 
        T[1:-1,0] = (dx**2*(Tn[2:,0]+Tn[0:-2,0])+dy**2*(2*Tn[1:-1,1]))/(2*(dx**2+dy**2)) 
        T[0,1:-1] = (dx**2*(2*Tn[1,1:-1])+dy**2*(Tn[0,2:]+Tn[0,0:-2]))/(2*(dx**2+dy**2)) #dT/dy = 0 @ y=0
        T[0,0] = (dx**2*(2*Tn[1,0])+dy**2*(2*Tn[0,1]))/(2*(dx**2+dy**2)) #dT/dy = 0 @ y=0
        T[:,-1] = 0       ##T = 0 @ x = 1.0
        T[-1,:] = 1.0    ##T = 1 @ y = 1.0
#        l1norm = (sum(abs(T[:])-abs(Tn[:])))/sum(abs(Tn[:]))
        l1norm = sum(sum(abs(T[:,:]-Tn[:,:])))/sum(sum(abs(Tn[:,:])))
    
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

    
def plot3D(x, y, p):
    fig = plt.figure(figsize=(11,7), dpi=100)
    ax = fig.gca(projection='3d')
    X,Y = np.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,p[:], rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    ax.view_init(30,225)
    plt.xlabel('x-values')
    plt.ylabel('y-values')
    
def subplot3D(x,y,p,Npx=1,Npy=1,Cp=1, title=''):
    if (Cp==1):
        fig = plt.figure()
    else:
        fig=plt.gcf()
    
    ax = fig.add_subplot(Npx, Npy, Cp, projection='3d')
        
    X,Y = np.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,p[:], rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False)
    plt.title(title)
    plt.xlabel('x-values')
    plt.ylabel('y-values')
    ax.view_init(30,225)
    

T = np.zeros((Ny+1,Nx+1))
T=laplace_direct_xdir(N)

Ti = np.zeros((Ny+1,Nx+1))
Ti=laplace2d(Ti, h, h, 0.00001)


x = np.linspace(0, width, Nx+1)
y = np.linspace(0, height, Ny+1)
X,Y = np.meshgrid(x, y)

Ta=T_analytical(X,Y)

# plot3D(x,y,Ta)
# plt.title('Analytical solution')
# 
# plot3D(x,y,T)
# plt.title('numerical')
# 
# plot3D(x,y,Ti)
# plt.title('iterative')


subplot3D(x,y,Ta,Npx=2,Npy=2,Cp=1,title='analytic')
subplot3D(x,y,T,Npx=2,Npy=2,Cp=2,title='numerical')
subplot3D(x,y,Ti,Npx=2,Npy=2,Cp=3,title= 'iterative solver')


plt.show()
plt.close()

