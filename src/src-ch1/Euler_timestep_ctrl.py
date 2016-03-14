# coding: utf-8
# src-ch1/Euler_timestep_ctrl.py;DragCoefficientGeneric.py @ git@lrhgit/tkt4140/src/src-ch1/DragCoefficientGeneric.py;
from DragCoefficientGeneric import cd_sphere    
from matplotlib.pyplot import *
import numpy as np

# change some default values to make plots more readable 
LNWDT=2; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT



g = 9.81      # Gravity m/s^2
d = 41.0e-3   # Diameter of the sphere
rho_f = 1.22  # Density of fluid [kg/m^3]
rho_s = 1275  # Density of sphere [kg/m^3]
nu = 1.5e-5   # Kinematical viscosity [m^2/s]
CD = 0.4      # Constant drag coefficient

def f(z, t):
    """2x2 system for sphere with constant drag."""
    zout = np.zeros_like(z)
    alpha = 3.0*rho_f/(4.0*rho_s*d)*CD
    zout[:] = [z[1], g - alpha*z[1]**2]
    return zout 

# define euler scheme
def euler(func,z0, time):
    """The Euler scheme for solution of systems of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""
    
    z = np.zeros((np.size(time),np.size(z0)))
    z[0,:] = z0

    for i in range(len(time)-1):
        dt = time[i+1]-time[i]
        z[i+1,:]=z[i,:] + np.asarray(func(z[i,:],time[i]))*dt
    return z

def v_taylor(t):
#    z = np.zeros_like(t)
    v = np.zeros_like(t)
    alpha = 3.0*rho_f/(4.0*rho_s*d)*CD
    v=g*t*(1-alpha*g*t**2)
    return v
     
# main program starts here

T = 10  # end of simulation
N = 10  # no of time steps


z0=np.zeros(2)
z0[0] = 2.0

# Prms for the analytical solution 
k1 = np.sqrt(g*4*rho_s*d/(3*rho_f*CD))
k2 = np.sqrt(3*rho_f*g*CD/(4*rho_s*d))

Ndts = 4  # Number of times to divide the dt by 2 
legends=[]
error_diff = []

for i in range(Ndts+1):
    time = np.linspace(0, T, N+1)
    ze = euler(f, z0, time)     # compute response with constant CD using Euler's method
    v_a = k1*np.tanh(k2*time)   # compute response with constant CD using analytical solution
    abs_error=np.abs(ze[:,1] - v_a)
    log_error = np.log2(abs_error[1:])
    max_log_error = np.max(log_error)
    
    plot(time[1:], log_error)
    legends.append('Euler scheme: N ' + str(N) + ' timesteps' )
    N*=2
    if i > 0:
        error_diff.append(previous_max_log_err-max_log_error)

    previous_max_log_err = max_log_error
    
print 'Approximate order of scheme n =', np.mean(error_diff)
print 'Approximate error reuduction by dt=dt/2:', 1/2**(np.mean(error_diff))


# plot analytical solution
# plot(time,v_a)
# legends.append('analytical')

# fix plot 
legend(legends, loc='best', frameon=False)
xlabel('Time [s]')
#ylabel('Velocity [m/s]')
ylabel('log2-error')
#savefig('example_euler_timestep_study.png', transparent=True)
show()
