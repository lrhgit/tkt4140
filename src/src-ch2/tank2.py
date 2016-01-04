# chapter2/src-ch2/tank2.py;ODEschemes.py @ git@lrhgit/tkt4140/allfiles/digital_compendium/chapter2/src-ch2/ODEschemes.py;

from ODEschemes import euler, heun, rk4
from numpy import cos, sin
import numpy as np
from matplotlib.pyplot import *

# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)


def tank2(y, x):
    """Differential equation for the displacement w in a cylindrical tank with linearly varying wall-thickness
    Args:
        y(array): an array containg w and its derivatives up to third order.
        x(array): independent variable
    Returns:
        dwdx(array): RHS of the system of first order differential equations 
    """
    z = 1-alpha*x
    dwdx = np.zeros_like(y)
    dwdx[0] = y[1]
    dwdx[1] = y[2]
    dwdx[2] = y[3]
    temp = (6*alpha/z)*y[3]- (6*alpha**2/z**2)*y[2]
    dwdx[3] = temp - 4*beta4*y[0]/z**2 - 4*beta4*(1-x)/z**3
    
    return dwdx

R = 8.5 # radius [m]
H = 7.95 # height [m]
t0 = 0.35 # thickness [m]
t1 = 0.1 # thickness [m]
ny = 0.2 # poissons number 
beta = H*(3*(1-ny**2)/(R*t0)**2)**0.25
beta4 = beta**4
alpha = (t0-t1)/t0
N = 100
print "beta: ", beta, "alpha", alpha

X = 1.0
x = np.linspace(0,X,N + 1) 

solverList = [euler, heun, rk4] #list of solvers
solver = solverList[2] # select specific solver

# shoot:
s = np.array([0, 0, 1])
r = np.array([0, 1, 0])
phi = np.zeros(3)
psi = np.zeros(3)
for k in range(3):
    y0 = np.array([0,0,s[k],r[k]])
    y = solver(tank2,y0,x)
    phi[k] = y[-1, 2]
    psi[k]=y[-1, 3]
    
phi1, phi2, phi3, psi1, psi2, psi3 = phi[0], phi[1], phi[2], psi[0], psi[1], psi[2]
# calculate correct r and s
nev = (psi3-psi1)*(phi2-phi1) - (phi3-phi1)*(psi2-psi1)
rstar = (phi3*psi1-psi3*phi1)/nev
sstar = (psi2*phi1 - phi2*psi1)/nev

print 'rstar', rstar, 'sstar', sstar

# solve with correct initial values
y0 = np.array([0, 0, sstar, rstar])
y = solver(tank2,y0,x)

legends=[] # empty list to append legends as plots are generated
plot(x,-y[:,3]/beta**2)
plot(x,-y[:,2]/beta)
legends.append('v(x)/B2')
legends.append('m(x)/B')

# Add the labels
legend(legends,loc='best',frameon=False) # Add the legends
ylabel('v, m')
xlabel('x')
grid(b=True, which='both', color='0.65',linestyle='-')

show()
