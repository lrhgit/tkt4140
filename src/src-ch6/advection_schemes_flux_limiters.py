# ../Kap6/advection_schemes.py
import Advection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy import interpolate
from numpy import where
from math import sin

LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT


# function defining the initial condition

def init_step(x):
    """Assigning a value of 1.0 for values less than 0.1"""
    f = np.zeros_like(x)
    f[np.where(x <= 0.1)] = 1.0
    return f

def init_sine2(x):
    """A smooth sin^2 function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = 0.05
    x_right = 0.3
    f = where((x>x_left) & (x<x_right), np.sin(np.pi*(x-x_left)/(x_right-x_left))**2,f) 
    return f

def init_box(x):
    """A smooth sin^2 function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = dx
    x_right = 0.2
    f = where((x>x_left) & (x<x_right), 1,f) 
    return f

def init_box_sine(x):
    """A smooth sin^2 function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = dx
    x_right = 0.2
    f = where((x>x_left) & (x<x_right), 1,f)
    x_left2 = 0.5
    x_right2 = 0.75
    f = where((x>x_left2) & (x<x_right2), np.sin(np.pi*(x-x_left2)/(x_right2-x_left2))**2,f) 
    return f


# Discretize

a = 1.0 # wave speed
tmin, tmax = 0.0, 2.0 # start and stop time of simulation
xmin, xmax = 0.0, 2.0 # start and end of spatial domain
Nx = 200 # number of spatial points
c = 0.8 # courant number, need c<=1 for stability





x = np.linspace(xmin, xmax, Nx+1) # discretization of space
dx = float((xmax-xmin)/Nx) # spatial step size
dt = c/a*dx # stable time step calculated from stability requirement
Nt = int((tmax-tmin)/dt) # number of time steps
time = np.arange(tmin, tmax, dt) # discretization of time


f = init_box
limiters = ['lax_wendroff', 'superbee'] # 'upwind', 'minmod', 'van_leer', 
solvers = []
for limiter in limiters:
    solvers.append(Advection.FluxLimiters(a, dx, dt, c, limiter))
    #solvers.append(Advection.FluxLimiters2(a, dx, dt, c, limiter))


u_solutions=np.zeros((len(solvers),len(time),len(x)))
uanalytical = np.zeros((len(time), len(x))) # holds the analytical solution


    
for k, solver in enumerate(solvers): # Solve for all solvers in list
    u = f(x)
    un = np.zeros((len(time), len(x))) # holds the numerical solution
    un[0:] = f(x)
    uanalytical[0,:] = f(x)
    for i, t in enumerate(time[1:]):
        i = i + 1
        if k==0:
            uanalytical[i,:] = f(x-a*t) # compute analytical solution for this time step
            
        u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
        
        u[1:-1] = solver.solve(u[:]) # calculate numerical solution of interior
        u[-1] = u_bc(x[-1] - a*dt) # interpolate along a characteristic to find the boundary value
        
        un[i,:] = u[:] # storing the solution for plotting
    
    u_solutions[k,:,:] = un



### Animation 
 
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(xmin,xmax), ylim=(np.min(u_solutions), np.amax(u_solutions)*1.1))

lines=[]     # list for plot lines for solvers and analytical solutions
legends=[]   # list for legends for solvers and analytical solutions

for solver in solvers:
    line, = ax.plot([], [])
    lines.append(line)
    legends.append(solver.name())

line, = ax.plot([], []) #add extra plot line for analytical solution
lines.append(line)
legends.append('Analytical')

plt.xlabel('x-coordinate [-]')
plt.ylabel('Amplitude [-]')
plt.legend(legends, loc='best', frameon=False)
 
# initialization function: plot the background of each frame
def init():
    for line in lines:
        line.set_data([], [])
    return lines,

# animation function.  This is called sequentially
def animate(i):
    for k, line in enumerate(lines):
        if (k==0):
            line.set_data(x, un[i,:])
        else:
            line.set_data(x, uanalytical[i,:])
    return lines,

def animate_alt(i):
    for k, line in enumerate(lines):
        if (k==len(lines)-1):
            line.set_data(x, uanalytical[i*jump,:])
        else:
            line.set_data(x, u_solutions[k,i*jump,:])
    return lines,

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
jump = 2
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate_alt, init_func=init, frames=Nt/jump, interval=100, blit=False)
#anim.save('../figs/step_flux_limiter.mov', writer=writer)
 
plt.show()
