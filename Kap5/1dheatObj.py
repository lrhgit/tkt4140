# The equation solved is the parabolic equaiton
#
#       du     d  du
#       -- = k -- --
#       dt     dx dx
#
# along with boundary conditions

import matplotlib.pyplot as plt
import numpy as np


class Grid1d:
    """A simple grid class for grid information and the solution."""
    def __init__(self, nx=10, xmin=0.0, xmax=1.0):
        self.xmin, self.xmax = xmin, xmax
        self.dx = float(xmax-xmin)/(nx)
        self.nx = nx
        self.u = np.zeros((nx+1, 1), 'd')
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
        print 'timestep = ',dt

        m=round((tmax-tmin)/dt) # Number of temporal intervals
        print 'm = ',m
        time=np.linspace(tmin,tmax,m)

        print 'time = ', time[0], 'and ', time[-1]
        

        nOutputInt=2 #output every nOutputInt iteration
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
                
        plt.pause(2)
        plt.ion()
        plt.close()

#        return u

    def solve(self, tmin, tmax):
        return self.solver(tmin,tmax)

    def initialize(self,U0=1.0):
        self.grid.u[0] = U0
        

## Main program

## Make a grid
nx=30
L=1.0
mg=Grid1d(nx,0,L)

mg2=Grid1d(nx/2,0,L)

## Make a solver
ms=HeatSolver1d(mg,scheme='explicit',k=1.0,r=0.55)
ms2=HeatSolver1d(mg2,scheme='explicit',k=1.0,r=0.3)


## Find the solution
U0=1.0
ms.initialize(U0=U0)
ms2.initialize(U0=U0)

(tmin, tmax)=(0,0.025)
ms.solve(tmin,tmax)
ms2.solve(tmin,tmax)

plt.plot(ms.grid.x,ms.grid.u,ms2.grid.x,ms2.grid.u)

plt.show()
plt.pause(5)
plt.close()

print 'done'



