# src-ch1/MMS_example1.py; ODEschemes.py @ git@lrhgit/tkt4140/src/src-ch1/ODEschemes.py;

import numpy as np
import matplotlib.pyplot as plt
from ODEschemes import euler, heun, rk4
from math import sqrt, pi
from sympy import symbols, diff, lambdify, sin

def func(y,x):
    """ Function that returns the du/dt of the differential equation
            dy/dx = f(y,x)
            
        Args:
            y(float): solutian array y(x)
            x(float): current time

        Returns:
            out(float): f(y,x)
    """
    return ffunc(y,x)

#### Use sympy to calculate the source term g based on the operator d/dx, and um ####
x = symbols('x')
yt = symbols('yt')
y_m = sin(x)
beta = -1.
f = beta*yt
f_m = diff(y_m, x) - beta*y_m
yfunc = lambdify(x, y_m, np) # create python function of the manufactured solution.
ffunc = lambdify((yt,x), f+f_m, np) # create python function of the source term g


x0, xend = 0, 2*pi # domain
y0 = yfunc(x0) # initial value

schemes = [euler, heun, rk4] # list of solvers imported from ODEschemes.py. 
schemes_order = {} # dict to be filled in with order approximations for all schemes

Ndts = 5 # Number of times to refine timestep in convergence test

for scheme in schemes: # iterate through all schemes
    N = 10    # initial number of time steps
    order_approx = [] # start of with empty list of orders for all schemes
    
    for i in range(Ndts+1):
        time = np.linspace(x0, xend, N+1)
        y = scheme(func, y0, time) # . e.g: euler(func, z0, time) 
        max_error = max(abs(y[:,0] - yfunc(time))) # calculate infinity norm 
        
        if i > 0: # Compute the observed order of accuracy based on two error norms
            order = np.log(previous_max_error/max_error)/np.log(2)
            order_approx.append(order)
            
        previous_max_error = max_error
        N *=2 # double number of timestep (halve dt)
    
    schemes_order[scheme.func_name] = order_approx #listof orders assigned to scheme

#### Plotting ####
N = N/2**Ndts
N_list = [N*2**i for i in range(1, Ndts+1)]
N_list = np.asarray(N_list)

plt.figure()
for key in schemes_order: # Iterate all keys in dictionary schemes_order 
    plt.plot(N_list, (np.asarray(schemes_order[key])))

#Plot theoretical n for 1st, 2nd and 4th order schemes
plt.axhline(1.0, xmin=0, xmax=N, linestyle=':', color='k')
plt.axhline(2.0, xmin=0, xmax=N, linestyle=':', color='k')
plt.axhline(4.0, xmin=0, xmax=N, linestyle=':', color='k')
plt.xticks(N_list, rotation=-70)
legends = schemes_order.keys()
legends.append('theoretical') 
plt.legend(legends, loc='best', frameon=False)
plt.title('Observed order of accuracy and MMS')
plt.xlabel('Number of time_steps')
plt.ylabel('Scheme order approximation')
plt.axis([0, max(N_list), 0, 5])
#plt.savefig('../figs/MMSExample1.png') # transparent=True
plt.show()

