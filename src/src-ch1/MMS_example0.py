# src-ch1/MMS_example1.py; ODEschemes.py @ git@lrhgit/tkt4140/src/src-ch1/ODEschemes.py;

import numpy as np
import matplotlib.pyplot as plt
from ODEschemes import euler, heun, rk4
from math import pi
from sympy import symbols, diff, lambdify, sin

def func(y,t):
    """ Function that returns the du/dt of the differential equation
            du/dt = g
            
        Args:
            y(float): solutian array u(t)
            t(float): current time

        Returns:
            yout(float): g(t)
    """
    return gfunc(t)

#### Main program starts here
t = symbols('t')
u = sin(t)
g = diff(u, t)

ufunc = lambdify(t, u, np) # create python function of the manufactured solution
gfunc = lambdify(t, g, np) # create python function of the source term g

N = 100
t0, tend = 0, 2*pi # domain
time = np.linspace(t0, tend, N)

u_0 = ufunc(t0) # initial value

uNum = euler(func, u_0, time)
uM = ufunc(time)

#plotting:
plt.figure()
plt.plot(time, uM, 'k')
plt.plot(time, uNum, 'r--')
plt.legend(['Manufactured', 'Numerical'], frameon=False)
plt.xlabel('t')
plt.ylabel('u')
#plt.savefig('../fig-ch1/MMSExample0.png', transparent=True) # transparent=True
plt.show()




