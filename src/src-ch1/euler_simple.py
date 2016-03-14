# src-ch1/euler_simple.py

import numpy as np
import matplotlib.pylab as plt

""" example using eulers method for solving the ODE 
    y'(x) = f(x, y) = y
    y(0) = 1
    
    Eulers method:
    y^(n + 1) = y^(n) + h*f(x, y^(n)), h = dx
    """

N = 100
x = np.linspace(0, 1, N + 1)
h = x[1] - x[0] # steplength
y_0 = 1 # initial condition 
Y = np.zeros_like(x) # vector for storing y values
Y[0] = y_0 # first element of y = y(0)

for n in range(N):
    f = Y[n]
    Y[n + 1] = Y[n] + h*f

Y_analytic = np.exp(x)

# change default values of plot to make it more readable
LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT
plt.rcParams['font.size'] = FNT

plt.figure()
plt.plot(x, Y_analytic, 'b', linewidth=2.0)
plt.plot(x, Y, 'r--', linewidth=2.0)
plt.legend(['$e^x$', 'euler'], loc='best', frameon=False)
plt.xlabel('x')
plt.ylabel('y')
#plt.savefig('../fig-ch1/euler_simple.png', transparent=True)
plt.show()


