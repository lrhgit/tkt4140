# src/src-ch2/beam_deflect_shoot_varying.py;ODEschemes.py @ git@lrhgit/tkt4140/src/src-ch2/ODEschemes.py;

from ODEschemes import euler, heun, rk4
from numpy import cos, sin
import numpy as np
from matplotlib.pyplot import *

# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)

def f(y, x):
    """Governing differential equation on beam bending as a system of 1. order equations.
    Args:
        y(array): an array containg y and its derivatives up to second order. (RHS)
        x(array): space array
    Returns:
        dydx(array): dydx. RHS of reduced system of 1st order equations
    """
    yout = np.zeros_like(y)
    yout[:] = [y[1],-(1+(1+x**n)*y[0])]
    
    return yout


# === main program starts here ===

N = 10 # number of elements
L = 1.0 # half the length of the beam
x = np.linspace(0,L,N+1) # vector of 
theta = 1 # PL**2/EI
theta2 = theta**2
h0=1 # height of beam at x = 0
n=2 # polynomial order to go into h(x)
h=h0/(1+(x/L)**n)**(1/3.) # height of beam assuming constant width

solvers = [euler, heun, rk4]
solver = solvers[2] 

s = [0, 1] # guessed values
# === shoot ===
y01 = [s[0], 0] # initial values
y02 = [s[1], 0] 

m0 = solver(f, y01, x)
m1 = solver(f, y02, x)

phi0 = m0[-1,0]
phi1 = m1[-1,0]

sfinal = phi0/(phi0 - phi1) # find correct initial value of moment m(x=0) using secant method

y03 = [sfinal, 0] 

mfinal = solver(f,y03,x)
m = mfinal[:, 0] # extract moment from solution data

u = m -0.5*(1 - x**2) # final solution of dimensionless deflection
ua = (1/theta2)*(cos(theta*x)/cos(theta)-1)-(1-x**2)/2 #analytical solution with constant stifness

legendList=[] # empty list to append legends as plots are generated

plot(x, m, 'b',) 
legendList.append('m')
plot(x, u, 'r',) 
legendList.append('u')
## Add the labels
legend(legendList,loc='best',frameon=False) 
ylabel('m, u')
xlabel('x')
grid(b=True, which='both', axis='both',linestyle='-')
show()
