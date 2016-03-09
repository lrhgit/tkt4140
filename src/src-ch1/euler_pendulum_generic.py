# src-ch1/euler_pendulum_generic.py

import numpy as np
import matplotlib.pylab as plt
from math import pi

# define Euler solver
def euler(func, y_0, time):
    """ Generic implementation of the euler scheme for solution of systems of ODEs:
            y0' = y1
            y1' = y2
                .
                .
            yN' = f(yN-1,..., y1, y0, t)
            
            method:
            y0^(n+1) = y0^(n) + h*y1
            y1^(n+1) = y1^(n) + h*y2
                .
                .
            yN^(n + 1) = yN^(n) + h*f(yN-1, .. y1, y0, t)

        Args:
            func(function): func(y, t) that returns y' at time t; [y1, y2,...,f(yn-1, .. y1, y0, t)] 
            y_0(array): initial conditions
            time(array): array containing the time to be solved for
        
        Returns:
            y(array): array/matrix containing solution for y0 -> yN for all timesteps"""

    y = np.zeros((np.size(time), np.size(y_0)))
    y[0,:] = y_0

    for i in range(len(time)-1):
        dt = time[i+1] - time[i]
        y[i+1,:]=y[i,:] + np.asarray(func(y[i,:], time[i]))*dt

    return y

def pendulum_func(y, t):
    """ function that returns the RHS of the mathematcal pendulum ODE:    
        Reduction of higher order ODE:
        theta = y0
        theta' = y1
        theta'' = - theta = -y0
    
        y0' = y1
        y1' = -y0
        
        Args:
            y(array): array [y0, y1] at time t
            t(float): current time
        
        Returns:
            dy(array): [y0', y1'] = [y1, -y0]
        """
        
    dy = np.zeros_like(y)
    dy[:] = [y[1], -y[0]]
    return dy

N = 100
time = np.linspace(0, 2*pi, N + 1)
thetha_0 = [0.1, 0]

theta = euler(pendulum_func, thetha_0, time)
thetha = theta[:, 0]

thetha_analytic = thetha_0[0]*np.cos(time)

# change default values of plot to make it more readable
LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT
plt.rcParams['font.size'] = FNT

plt.figure()
plt.plot(time, thetha_analytic, 'b')
plt.plot(time, thetha, 'r--')

plt.legend([r'$\theta_0 \cdot cos(t)$', 'euler'], loc='best', frameon=False)
plt.xlabel('t')
plt.ylabel(r'$\theta$')
#plt.savefig('../fig-ch1/euler_pendulum.png', transparent=True)
plt.show()


