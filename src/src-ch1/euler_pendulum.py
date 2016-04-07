# src-ch1/euler_pendulum.py

import numpy as np
import matplotlib.pylab as plt
from math import pi

""" example using eulers method for solving the ODE:
    theta''(t) + thetha(t) = 0
    thetha(0) = theta_0
    thetha'(0) = dtheta_0
    
    Reduction of higher order ODE:
    theta = y0
    theta' = y1
    theta'' = - theta = -y0
    
    y0' = y1
    y1' = -y0
    
    eulers method:
    y0^(n + 1) = y0^(n) + h*y1, h = dt
    y1^(n + 1) = y1^(n) + h*(-y0), h = dt
    """

N = 100
t = np.linspace(0, 2*pi, N + 1)
h = t[1] - t[0] # steplength
thetha_0 = 0.1
y0_0 = thetha_0 # initial condition
y1_0 = 0
Y = np.zeros((2, N + 1)) # 2D array for storing y values

Y[0, 0] = y0_0 # apply initial conditions
Y[1, 0] = y1_0

for n in range(N):
    y0_n = Y[0, n]
    y1_n = Y[1, n]
    
    Y[0, n + 1] = y0_n + h*y1_n
    Y[1, n + 1] = y1_n - h*y0_n

thetha = Y[0, :]
thetha_analytic = thetha_0*np.cos(t)

# change default values of plot to make it more readable
LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT
plt.rcParams['font.size'] = FNT

plt.figure()
plt.plot(t, thetha_analytic, 'b')
plt.plot(t, thetha, 'r--')

plt.legend([r'$\theta_0 \cdot cos(t)$', 'euler'], loc='best', frameon=False)
plt.xlabel('t')
plt.ylabel(r'$\theta$')
#plt.savefig('../fig-ch1/euler_pendulum.png', transparent=True)
plt.show()


