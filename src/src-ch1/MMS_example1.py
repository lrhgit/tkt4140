# chapter1/src-ch1/MMS_example1.py; ODEschemes.py @ git@lrhgit/tkt4140/allfiles/digital_compendium/chapter2/src-ch1/ODEschemes.py;

import numpy as np
from ODEschemes import euler, heun, rk4
from math import sqrt, pi
from sympy import symbols, diff, lambdify, sin
from matplotlib.pyplot import plot, show, legend, hold,rcParams,rc, figure, axhline, close,\
    title, xticks, xlabel, ylabel, savefig, axis, grid
#### Use sympy to calculate the source term g, based on the differential operator d/dx, and the manufactured solution u = sin(t)####
t = symbols('t')
u = sin(t)
g = diff(u, t)
ufunc = lambdify(t, u, np) # create python function of the manufactured/chosen solution. "np" makes it able to handle arrays as input
gfunc = lambdify(t, g, np) # create python function of the source term g

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

t0, tend = 0, 2*pi # domain
z0 = ufunc(t0) # initial value

schemes = [euler, heun, rk4] # list of solvers imported from ODEschemes.py. each of which is a function
schemes_order = {} # empty dictionary. to be filled in with order approximations for all schemes

Ndts = 5 # Number of times to refine timestep in convergence test

for scheme in schemes: # iterate through all schemes
    N = 10    # initial number of time steps
    order_approx = [] # start of with empty list of orders for all schemes
    
    for i in range(Ndts+1):
        time = np.linspace(t0, tend, N+1)
        z = scheme(func, z0, time) # Solve the ODE by calling the scheme with arguments. e.g: euler(func, z0, time) 
        max_error = max(abs(z[:,0] - ufunc(time))) # calculate infinity norm of the error
        
        if i > 0: # Compute the observed order of accuracy based on error from two grids; error1(h), error2(h/2)
            order = np.log(previous_max_error/max_error)/np.log(2)
            order_approx.append(order)
            
        previous_max_error = max_error
        N *=2 # double number of timestep (halve dt)
    
    schemes_order[scheme.func_name] = order_approx # Add a key:value pair to the dictionary. e.g: "euler":[orderapprox1, orderapprox2, ..., orderapproxNtds]

#### Plotting ####
N = N/2**Ndts
N_list = [N*2**i for i in range(1, Ndts+1)]
N_list = np.asarray(N_list)

figure()
for key in schemes_order: # Iterate all keys in dictionary schemes_order. e.g: key = "euler" 
    plot(N_list, (np.asarray(schemes_order[key])))

# Plot theoretical n for 1st, 2nd and 4th order schemes
axhline(1.0, xmin=0, xmax=N, linestyle=':', color='k')
axhline(2.0, xmin=0, xmax=N, linestyle=':', color='k')
axhline(4.0, xmin=0, xmax=N, linestyle=':', color='k')
xticks(N_list, rotation=-70)
legends = schemes_order.keys()
legends.append('theoretical') 
legend(legends, loc='best', frameon=False)
title('Observed order of accuracy and MMS')
xlabel('Number of time_steps')
ylabel('Scheme order approximation')
axis([0, max(N_list), 0, 5])
#savefig('../figs/MMSExample1.png') # transparent=True
show()

