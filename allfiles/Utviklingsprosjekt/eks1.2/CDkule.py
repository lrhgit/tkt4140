from DragCoefficientGeneric import cd_sphere    
import odespy
from matplotlib.pyplot import *
import numpy as np
#change some default values to make plots more readable on the screen
LNWDT=5; FNT=25
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT

g = 9.81      # Gravity m/s^2
d = 41.0e-3     # Diameter of the sphere
rho_f = 1.22  # Density of fluid [kg/m^3]
rho_s = 1275  # Density of sphere [kg/m^3]
nu = 1.5e-5   # Kinematical viscosity [m^2/s]

def f(z, t):
    """2x2 syst for sphere with constant drag."""
    zout = np.zeros_like(z)
    CD = 0.5
    alpha = 3.0*rho_f/(4.0*rho_s*d)*CD
    zout[:] = [z[1], g - alpha*z[1]**2]
    return zout 


def f2(z, t):
    """2x2 syst for sphere with Re-dependent drag."""
    zout = np.zeros_like(z)
    v = abs(z[1]) 
    Re = v*d/nu
    CD = cd_sphere(Re)
    alpha = 3.0*rho_f/(4.0*rho_s*d)*CD
    zout[:] = [z[1], g - alpha*z[1]**2]
    return zout

# Main program starts here
from numpy import linspace
T = 30  # end of simulation
N = 50  # no of time steps
time = linspace(0, T, N+1)

solvers=[]
solvers.append(odespy.RK3(f2)) 
solvers.append(odespy.RK4(f2)) 

def euler(func,z0, time):
    dt = time[1]-time[0]
    z = np.zeros((np.size(time),2))
    z[0,:] = z0

    for i, t in enumerate(time[1:]):
        z[i+1,:]=z[i,:] + func(z[i,:],t)*dt

    return z
        

legends=[]
linet=['r-',':','.','-.','--']

z0 = [2.0, 0.0]
for i, solver in enumerate(solvers):
    solver.set_initial_condition(z0)
    z, t = solver.solve(time)
    plot(t,z[:,1],linet[i])
    legends.append(str(solver))

z0e=np.zeros(2)
z0e[0] = 2.0
time2 = linspace(0, T, N+1)
ze = euler(f2,z0e,time)
plot(time,ze[:,1])
legends.append('Euler')
legend(legends)

show()

