# src-ch2/beam_deflect_shoot_constant.py;ODEschemes.py @ git@lrhgit/tkt4140/src/src-ch2/ODEschemes.py;

from ODEschemes import euler, heun, rk4
from numpy import cos, sin
import numpy as np
from matplotlib.pyplot import *

# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT

def Sys1(y, x):
    """ system 1 of Governing differential equation on beam bending with constant cross section
    Args:
        y(array): an array containg y and its derivatives up to second order. (RHS)
        x(array): an array 
    Returns:
        yout(array): dydx RHS of reduced system of 1st order equations
    """
    yout = np.zeros_like(y)
    yout[:] = [y[1],-theta2*(y[0]+0.5*(1-x**2))]
    
    return yout

def Sys2(y, x):
    """ system 2 of Governing differential equation on beam bending with constant cross section
    Args:
        y(array): an array containg y and its derivatives up to second order. (RHS)
        x(array): an array 
    Returns:
        yout(array): dydx RHS of reduced system of 1st order equations
    """
    yout = np.zeros_like(y)
    yout[:] = [y[1], -theta2*y[0]]
    return yout

# === main program starts here ===

N = 20 # number of elements
L = 1.0 # half the length of the beam
x = np.linspace(-L, L, N + 1) # allocate space
theta = 1 # PL**2/EI
theta2 = theta**2

solverList = [euler, heun, rk4] 
solver = solverList[2] 

s = [0, 1] # guessed values
# === shoot ===
y0Sys1 = [0, s[0]] # initial values of u and u'
y0Sys2 = [0, s[1]] 

u0 = solver(Sys1, y0Sys1,x)
u1 = solver(Sys2, y0Sys2,x)

u0 = u0[:,0] # extract deflection from solution data
u1 = u1[:,0] 

u = u0 -(u0[-1]/u1[-1])*u1 # interpolate to find correct solution
ua = (1/theta2)*(cos(theta*x)/cos(theta) - 1)-(1 - x**2)/2 # analytical solution

legendList=[] # empty list to append legends as plots are generated

plot(x,u,'y',) 
plot(x,ua,'r:')
legendList.append('shooting technique')
legendList.append('analytical solution')
## Add the labels
legend(legendList,loc='best',frameon=False) 
ylabel('u')
xlabel('x')
grid(b=True, which='both', color='0.65',linestyle='-')

show()
