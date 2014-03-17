# The equation solved is the parabolic equaiton
#
#       du     d  du
#       -- = k -- --
#       dt     dx dx
#
# along with boundary conditions

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import scipy.sparse
import scipy.sparse.linalg


class Grid1d:
    """A simple grid class for grid information and the solution."""
    def __init__(self, nx=10, xmin=0.0, xmax=1.0):
        self.xmin, self.xmax = xmin, xmax
        self.dx = float(xmax-xmin)/(nx)
        self.nx = nx                           # Number of dx  
        self.u = np.zeros((nx+1, 1), 'd')      # Number of x-values is nx+1
        self.x = np.linspace(xmin,xmax,nx+1)


class HeatSolver1d:
    """A simple 1dheat equation solver that can use different schemes to solve the problem."""
    def __init__(self, grid, scheme='explicit',k=1.0,r=0.5):
        self.grid = grid
        self.setSolverScheme(scheme)
        self.k = k
        self.r = r

    def setSolverScheme(self, scheme='explicit'):
        """Sets the scheme to be which should be one of ['slow', 'explicit', 'implicit']."""
        if scheme == 'slow':
            self.solver = self.pythonExplicit
        elif scheme == 'explicit':
            self.solver = self.numpyExplicit
        elif scheme == 'implicit':
            self.solver = self.numpyImplicit
        else:
            self.solver = self.numpyExplicit

    def numpyExplicit(self, tmin, tmax):
        """Solve equation for all t in time step using a NumPy expression."""
        g = self.grid
        k = self.k       #Diffusivity
        r = self.r       #Numerical Fourier number
        u, x, dx  = g.u, g.x, g.dx
        xmin, xmax = g.xmin, g.xmax
        
        dt=r*dx**2/k     #Compute timestep based on Fourier number, dx and diffusivity

        m=round((tmax-tmin)/dt) # Number of temporal intervals
        time=np.linspace(tmin,tmax,m)
        
        nOutputInt=5 #output every nOutputInt iteration
        i = 0        #iteration counter

        #Plot initial solution
        fig = plt.figure()
        ax=fig.add_subplot(111)
        Curve, = ax.plot( x, u[:], '-')
        ax.set_xlim([xmin,xmax])
#        ax.set_ylim([umin,umax])
        plt.xlabel('x')
        plt.ylabel('Velocity')

        plt.ion()
        plt.show()

        for t in time:
            i+=1

            u[1:-1] =  r*(u[0:-2]+ u[2:]) + (1.0-2.0*r)*u[1:-1]
    
            if (np.mod(i,nOutputInt)==0): #output every nOutputInt iteration
                Curve.set_ydata(u)
                plt.pause(.005)
                plt.title( 'step = %3d; t = %f' % (i,t ) )
        
        g.u = u        
        plt.pause(2)
        plt.ion()
        plt.close()

#        return u

    def numpyImplicit(self, tmin, tmax, theta=1.0):
        g = self.grid
        k = self.k       #Diffusivity
        r = self.r       #Numerical Fourier number
        u, x, dx  = g.u, g.x, g.dx
        xmin, xmax = g.xmin, g.xmax
        
        dt=r*dx**2/k     #Compute timestep based on Fourier number, dx and diffusivity
        print 'timestep = ',dt

        m=round((tmax-tmin)/dt) # Number of temporal intervals
        print 'm = ',m
        time=np.linspace(tmin,tmax,m)


        #Create matrix for sparse solver. Solve for interior values only (nx-1)
        diagonals=np.zeros((3,g.nx-1))   
        diagonals[0,:] = -r*theta                       #all elts in first row is set to 1
        diagonals[1,:] = 1+2.0*r*theta  
        diagonals[2,:] = -r*theta 
        As = sc.sparse.spdiags(diagonals, [-1,0,1], g.nx-1, g.nx-1,format='csc') #sparse matrix instance

        #Crete rhs array
        d=np.zeros((g.nx-1,1),'d')
        
        
        nOutputInt=10 #output every nOutputInt iteration
        i = 0        #iteration counter

        #Plot initial solution
        fig = plt.figure()
        ax=fig.add_subplot(111)
        Curve, = ax.plot( x, u[:], '-')
        ax.set_xlim([xmin,xmax])
#        ax.set_ylim([umin,umax])
        plt.xlabel('x')
        plt.ylabel('Velocity')

        plt.ion()
        plt.show()

        #Advance in time an solve tridiagonal system for each t in time
        for t in time:
            i+=1
            d[:] = (1-theta)*(u[0:-2]-2*u[1:-1]+u[2:])  
            #u[1:-1] = sc.sparse.linalg.spsolve(As,d) #theta=sc.linalg.solve_triangular(A,d)
            w= sc.sparse.linalg.spsolve(As,d) #theta=sc.linalg.solve_triangular(A,d)
            g.u[1:-1] = w[:,None]
            u = g.u
            if (np.mod(i,nOutputInt)==0): #output every nOutputInt iteration
                Curve.set_ydata(u)
                plt.pause(.005)
                plt.title( 'step = %3d; t = %f' % (i,t ) )
        
        
        plt.pause(2)
        plt.ion()
        plt.close()

           
    def solve(self, tmin, tmax):
        return self.solver(tmin,tmax)

    def initialize(self,U0=1.0):
        self.grid.u[0] = U0
        

## Main program
print 'hello world'
## Make a grid
nx=30
L=1.0
mg=Grid1d(nx,0,L)
mg2=Grid1d(nx/2,0,L)
mg3=Grid1d(nx,0,L)

## Make a solver
ms=HeatSolver1d(mg,scheme='explicit',k=1.0,r=0.5)
ms2=HeatSolver1d(mg2,scheme='explicit',k=1.0,r=0.3)
ms3=HeatSolver1d(mg3,scheme='implicit',k=1.0,r=0.6)

## Find the solution(s)
U0=1.0
ms.initialize(U0=U0)
ms2.initialize(U0=U0)
ms3.initialize(U0=U0)

(tmin, tmax)=(0,0.025)
#ms.solve(tmin,tmax)
#ms2.solve(tmin,tmax)
ms3.solve(tmin,tmax)

#plt.plot(ms.grid.x,ms.grid.u,ms2.grid.x,ms2.grid.u,ms3.grid.x,ms3.grid.u)
#plt.plot(ms.grid.x,ms.grid.u,ms2.grid.x,ms2.grid.u)
#plt.legend(['explicit r=0.5','explicit r=0.3','implicit r=0.6'])
#plt.show()
#plt.pause(5)
#plt.close()
#print 'u',ms.grid.u
#print ms3.grid.u
print 'done'



