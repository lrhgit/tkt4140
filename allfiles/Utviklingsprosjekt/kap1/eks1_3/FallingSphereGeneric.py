import odespy
from matplotlib.pyplot import legend, plot, show
import matplotlib
import numpy as np
#change some default values to make plots more readable on the screen
LNWDT=5; FNT=11
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


def euler(func,z0, time):
    """The Euler scheme for solution of systems of of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""
    
    dt = time[1]-time[0]
    z = np.zeros((np.size(time),2))
    z[0,:] = z0

    for i, t in enumerate(time[1:]):
        z[i+1,:]=z[i,:] + func(z[i,:],t)*dt

    return z

def heun(func,z0, time):
    """The Heun scheme for solution of systems of of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""
    
    dt = time[1]-time[0]
    z = np.zeros((np.size(time),2))
    z[0,:] = z0
    zp = np.zeros_like(z0)
    
    for i, t in enumerate(time[1:]):
        zp = z[i,:] + func(z[i,:],t)*dt   # Predictor step
        z[i+1,:] = z[i,:] + (func(z[i,:],t) + func(zp,t+dt))*dt/2.0 # Corrector step

    return z
        
# Main program starts here
from numpy import linspace
T = 30  # end of simulation
N = 25  # no of time steps
time = linspace(0, T, N+1)

solvers=[]
solvers.append(odespy.RK3(f2)) 
solvers.append(odespy.RK4(f2)) 

legends=[]

z0=np.zeros(2)
z0[0] = 2.0

for i, solver in enumerate(solvers):
    solver.set_initial_condition(z0)
    z, t = solver.solve(time)
    plot(t,z[:,1])
    legends.append(str(solver))


scheme_list  = [euler, heun]

for scheme in scheme_list:
    z = scheme(f2,z0,time)
    plot(time,z[:,1])
    legends.append(scheme.func_name)
    
# 
# A more direct but less generic approach 
#   
# ze = euler(f2,z0,time)
# zh = heun(f2,z0,time)
#     
# plot(time,ze[:,1])
# legends.append('Euler')
# 
# plot(time,zh[:,1])
# legends.append('Heun')

legend(legends, loc='best', frameon=False)

show()
