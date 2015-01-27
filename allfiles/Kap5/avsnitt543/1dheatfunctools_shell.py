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

# change some default values to make plots more readable on the screen
LNWDT = 5; FNT = 15
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT

def explicit_python_solver(u_left=1.0, u_right=0.0, nx=20, r=0.5, xmin=0.0, xmax=1.0, tmin=0.0, tmax=1.0, k=1.0):
    dx = float(xmax-xmin)/nx
    u = np.zeros((nx+1, 1), 'd')
    u_old = np.zeros((nx+1, 1), 'd')
    x = np.linspace(xmin,xmax,nx+1)

    # set boundary conditions
    u[0] = u_left
    u_old[-1] = u_right
    
    dt = r*dx**2/k     # compute timestep based on Fourier number, dx and diffusivity
    
    m = round((tmax-tmin)/dt) # number of temporal intervals
    time = np.linspace(tmin,tmax,m)
    
    # advance in time
    #for t in time:    # uncomment this
        #####################################
        ### write integration scheme here ###
        #####################################
            
    return x, u

import functools        

## Main program starts here

nx = 20 # number of nodes
L  = 1.0    # length of beam
tmax = 0.025    # time length
theta = 0.75    # parameter for implicitness: theta=0.5 Crank-Nicholson, theta=1.0 fully implicit

call = functools.partial
solverlist = [explicit_python_solver]
solvers = [
    [call(solverlist[0], u_left=100.0, u_right=0.0, nx=nx, r=0.5, xmin=0.0, xmax=L, tmin=0.0, tmax=tmax, k=1.0)],
    ['explicit python solver']
    ]

lstyle = ['r-', ':', '.', '-.', '--']
legends = solvers[1]
i = 0
for solve in solvers[0]:
    tic = time.time()
    x, u = solve()
    toc = time.time()
    cputime = toc - tic
    print cputime, 'seconds process time for', legends[i]
    plt.plot(x,u,lstyle[i])
    i += 1
plt.legend(legends)
plt.title('Temperature field')
plt.xlabel('Position on beam')
plt.ylabel('Temperature')
#plt.savefig('1dheat0_025.pdf')
plt.show()