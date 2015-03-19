'''
Created on Mar 5, 2015
 
Solution of Heat equation d2T/dx2 + d2T/dy2 = 0
Robin BCs
          dT/dy = 0 at y=0 and y=1
          T=0 at x=0
          T=y at x=2.0

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
from numpy import cosh, cos, sum, abs, log2
from pprint import pprint
from matplotlib.pyplot import rcParams,rc
LNWDT=3; FNT=6
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 12}; rc('font', **font)



def T_analytical(x,y):
    Nmax = 15
    T = x/4.0
    
    for n in range(1,Nmax,2):
        npi=n*np.pi
        T+=-4*np.sinh(npi*x)*np.cos(npi*y)/(npi**2*np.sinh(2*npi))
        
    return T

    
def laplace2d(T, y, dx, dy, l1_eps):
   l1norm = 1
   Tn = np.empty_like(T)

   while l1norm > l1_eps:
       Tn = T.copy()
       T[1:-1,1:-1] = (dx**2*(Tn[2:,1:-1]+Tn[0:-2,1:-1])+dy**2*(Tn[1:-1,2:]+Tn[1:-1,0:-2]))/(2*(dx**2+dy**2)) 
       T[0,1:-1] = (dx**2*(2*Tn[1,1:-1])+dy**2*(Tn[0,2:]+Tn[0,0:-2]))/(2*(dx**2+dy**2)) #dT/dy = 0 @ y=0
       T[-1,1:-1] = (dx**2*(2*Tn[-2,1:-1])+dy**2*(Tn[-1,2:]+Tn[-1,0:-2]))/(2*(dx**2+dy**2)) #dT/dy = 0 @ y=1
   
       T[:,0] = 0        ##T = 0 @ x = 0
       T[:,-1] = y        ##T = y @ x = 2
       l1norm = (sum(abs(T[:])-abs(Tn[:])))/sum(abs(Tn[:]))
    
   return T
        


def laplace_directsolver_x_order(T,N,y):
    diagonals=np.zeros((5,N))
    diagonals[0,:]=1.0              # -Nx diagonal, the last Nx-elt are not used 
    diagonals[0,-2*Nx:]=2.0         # set the last Nx elt to 2, dT/dy=0 at y=1
      
    
    diagonals[1,:]=1.0               
    diagonals[1,Nx-1::Nx]=0.0       # every Nxth elt set to 0,staring at eqn Nx 
    
    diagonals[2,:]=-4.0              
    
    diagonals[3,:]=1.0               
    diagonals[3,Nx::Nx]=0.0          # every Nxth elt set to 0,staring at idx 
    
    diagonals[4,:]=1.0               # Nx-diagonal,the frst Nx-elts are not used 
    diagonals[4,Nx:2*Nx]=2.0         # set ther first Nx elt to 2, dT/dy=0 at y=0
    
    A=sc.sparse.spdiags(diagonals, [-Nx,-1,0,1,Nx], N, N,format='csc') #sparse matrix instance
#    print A.todense() 
    
    d=np.zeros(N)
    d[Nx-1::Nx]=-y # set the last Nx elts to -1.0
    
#     print d
    Tvector = sc.sparse.linalg.spsolve(A,d) #theta=sc.linalg.solve_triangular(A,d)
    
    Tmatrix=np.reshape(Tvector, (Ny,Nx), order='C')
    
    T = np.zeros((Ny,Nx+2))
    T[:,0]=0.0 # Left BC
    T[:,-1]=y  # Right BC
    
    T[:,1:-1] = Tmatrix
    
    return T

def laplace_directsolver_y_order(T,Nx,Ny,y):
    N=Nx*Ny
    diagonals=np.zeros((5,N))
    
    diagonals[0,:]=1.0   # -Ny diagonal => the last Ny elts are neglected            

    diagonals[1,:]=1.0          # -1 diagonal => the last elt is neglected
    diagonals[1,Ny-2::Ny]=2.0   #  dT/dy=0 at y=1, 1st elt to 2nd eqn                            
    diagonals[1,Ny-1::Ny]=0.0   #  dT/dy=0 at y=0, 
    
    
    diagonals[2,:]=-4.0 # Main diagonal                 
    
    diagonals[3,:]=1.0      # 1 diagonal => the first elt is neglected                      
    diagonals[3,1::Ny]=2.0  # dT/dy=0 at y=0                
    diagonals[3,Ny::Ny]=0.0 # dT/dy=0 at y=1.0                
    
    diagonals[4,:]=1.0     # Ny diagonal => the first Ny elts are neglected            

    #sparse matrix instance
    A=sc.sparse.spdiags(diagonals, [-Ny,-1,0,1,Ny], N, N,format='csc') 
#    print A.todense()

    
    d=np.zeros(N)
    d[-Ny:]=-y # set the last Ny elts to -y
     
    Tvector = sc.sparse.linalg.spsolve(A,d) #theta=sc.linalg.solve_triangular(A,d)
    
    # Reshape the Tvector to a matrix, Ny-vls first to comply with meshgrid
    Tmatrix=np.reshape(Tvector, (Ny,Nx), order='F') #first index changes fastest
    
    # Set the BC for the resulting array
    T[:,0]=0.0 # BC at x=0
    T[:,-1]=y  # BC at x=xmax
    T[:,1:-1] = Tmatrix # Fill in the computed solution in the matrix
    
    return T

        
def subplot3D(x,y,p,Npx=1,Npy=1,Cp=1, title=''):
    if (Cp==1):
        fig = plt.figure(figsize=(11,7), dpi=100)
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
    
#     ax.set_xlim(0,2)
#     ax.set_ylim(0,1)


def convergence_test(h=0.25,Lx=2.0,Ly=1.0,Nhs=6):
    
    order_approx=[]
    
    for i in range(Nhs):
        Nx=int((np.rint(Lx/h)-1)) # number of points in x-direction
        Ny=int((np.rint(Ly/h)+1)) # number of points in y-direction
        N=Nx*Ny

        x = np.linspace(0, Lx, Nx+2)
        y = np.linspace(0, Ly, Ny)

        # analytical
        X,Y = np.meshgrid(x,y)
        Ta = np.zeros((Ny,Nx+2))
        Ta=T_analytical(X,Y)
        T4 = np.zeros((Ny,Nx+2))
        T4=laplace_directsolver_y_order(T4,Nx,Ny,y)
        log2L1 = log2(abs(sum(abs(T4[:]-Ta[:])))/sum(abs(Ta[:])))

        if i > 0:
            order_approx.append(previous_L1norm-log2L1)
        
        previous_L1norm=log2L1
        h=h/2.0 
        
    plt.plot(range(1,Nhs), order_approx)
        
    return 

# Geometry
Lx=2.0
Ly=1.0
h=0.1


Nx=int((np.rint(Lx/h)-1)) # number of points in x-direction
Ny=int((np.rint(Ly/h)+1)) # number of points in y-direction
#print 'Nx=',Nx,' Ny=', Ny
N=Nx*Ny

x = np.linspace(0, Lx, Nx+2)
y = np.linspace(0, Ly, Ny)

# analytical
X,Y = np.meshgrid(x,y)
Ta=T_analytical(X,Y)

# iterative solver
T3 = np.zeros((Ny,Nx+2))
##boundary conditions
T3[:,0] = 0        ##p = 0 @ x = 0
T3[:,-1] = y        ##p = y @ x = 2
T3[0,:] = T3[1,:]        ##dp/dy = 0 @ y = 0
T3[-1,:] = T3[-2,:]    ##dp/dy = 0 @ y = 1
 
 
T3=laplace2d(T3, y, h, h, 0.0001)

T2 = np.zeros((Ny,Nx+2))
T2=laplace_directsolver_x_order(T2,N,y)
 
T4 = np.zeros((Ny,Nx+2))
T4=laplace_directsolver_y_order(T4,Nx,Ny,y)
# L1norm = (sum(abs(T4[:])-abs(Ta[:])))/sum(abs(Ta[:]))
# print 'L1n =', L1norm
#  
subplot3D(x,y,Ta,Npx=2,Npy=2,Cp=1,title='analytic')
subplot3D(x,y,T2,Npx=2,Npy=2,Cp=3,title='directsolver x-dir')
subplot3D(x,y,T4,Npx=2,Npy=2,Cp=2,title='directsolver y-dir')
subplot3D(x,y,T3,Npx=2,Npy=2,Cp=4,title= 'iterative solver')


plt.figure()
convergence_test()
plt.xlabel('h-level')
plt.ylabel('order-approx')

plt.show()
plt.close()






