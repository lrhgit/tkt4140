# The equation solved is the parabolic equaiton
#
#       du     d  du
#       -- = k -- --
#       dt     dx dx
#
# along with boundary conditions

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy as sc
import scipy.sparse
import scipy.sparse.linalg
import time

class Grid1d:
    """A simple grid class for grid information and the solution."""
    def __init__(self, nx, xmin, xmax):
        self.xmin, self.xmax = xmin, xmax
        self.dx = float(xmax-xmin)/(nx)
        self.nx = nx                           # number of dx  
        self.u = np.zeros((nx+1, 1), 'd')      # number of x-values is nx+1
        self.u_old = np.zeros((nx+1, 1), 'd')  # for use in pythonExplicit
        self.x = np.linspace(xmin,xmax,nx+1)


class HeatSolver1d:
    """A simple 1dheat equation solver that can use different schemes to solve the problem."""
    def __init__(self, grid, k, r, scheme, theta=0.5):
        self.grid = grid
        self.setSolverScheme(scheme)
        self.k = k
        self.r = r
        self.theta = theta # used for implicit solver only

    def setSolverScheme(self, scheme):
        """Sets the scheme to be which should be one of ['explicit', 'implicit']."""
        #if scheme == 'explicit':
        self.solver = self.pythonExplicit
        self.name = 'explicit'
        #elif scheme == 'implicit':
            #self.solver = self.numpyImplicit
            #self.name   = 'implicit'
        #else:
            #self.solver = self.numpyImplicit
            #self.name   = 'implicit'
    

    def pythonExplicit(self, tmin, tmax):
        """Solve equation for all t in time step using a Python expression."""
        g = self.grid
        k = self.k       # diffusivity
        r = self.r       # numerical Fourier number
        u, u_old, x, dx, n = g.u, g.u_old, g.x, g.dx, g.nx
        xmin, xmax = g.xmin, g.xmax
        
        dt=r*dx**2/k     # compute timestep based on Fourier number, dx and diffusivity

        m=round((tmax-tmin)/dt) # number of temporal intervals
        time=np.linspace(tmin,tmax,m)
        
        # advance in time
        #for t in time:    # uncomment this
            #####################################
            ### write integration scheme here ###
            #####################################     
   
    def solve(self, tmin, tmax):
        return self.solver(tmin,tmax)


    def initialize(self,U0,Uleft,Uright):
        self.grid.u[1::-1] = U0
        self.grid.u[0] = Uleft
        self.grid.u[-1] = Uright
        
        
## Main program starts here

nx = 20 # number of nodes
L  = 1.0    # length of beam

## make various solvers.
solvers=[]
solvers.append(HeatSolver1d(Grid1d(nx,0,L), scheme = 'explicit', k=1.0, r=0.2))
solvers.append(HeatSolver1d(Grid1d(nx,0,L), scheme = 'explicit', k=1.0, r=0.5))


U0=0.0    # initial values
Uleft=100.0 # boundary value on the left side
Uright=0.0  # boundary value on the right side
(tmin, tmax)=(0,0.025) # time period

## compute a solution for all solvers
for solver in solvers:
    solver.initialize(U0,Uleft,Uright)
    tic = time.time()
    solver.solve(tmin,tmax)
    toc = time.time()
    cputime=toc-tic
    print cputime, 'seconds process time for', solver.name, 'solver with r=', solver.r
lstyle=['r-',':','.','-.','--']
mylegends=[]
plt.figure()
i=0

for solver in solvers:
    plt.plot(solver.grid.x,solver.grid.u,lstyle[i])
    mylegends.append(str('%s r = %3.1f' % (solver.name, solver.r)))
    i+=1

# change some default values to make plots more readable on the screen
LNWDT=2; FNT=15
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT

plt.legend(mylegends)
#plt.savefig('1dheat.pdf')
plt.show()