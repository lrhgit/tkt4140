# src-ch2/beam_deflect_shoot.py;ODEschemes.py @ git@lrhgit/tkt4140/src/src-ch2/ODEschemes.py;

from ODEschemes import euler, heun, rk4
from numpy import cos, sin
import numpy as np
from matplotlib.pyplot import *

# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT



N=20
L = 1.0
y = np.linspace(0,L,N+1)

def dsfunction(phi0,phi1,s0,s1):
    if (abs(phi1-phi0)>0.0):   
        return    -phi1 *(s1 - s0)/float(phi1 - phi0)
    else:
        return 0.0


def f(z, t):
    """RHS for deflection of beam"""
    zout = np.zeros_like(z)
    zout[:] = [z[1],-alpha2*cos(z[0]),sin(z[0])]
    return zout 

alpha2 = 5.0 
beta=0.0 # Boundary value at y = L

solvers = [euler, heun, rk4] #list of solvers
solver=solvers[2] # select specific solver

# Guessed values
s=[2.5, 5.0]

z0=np.zeros(3)

z0[1] = s[0]
z = solver(f,z0,y)
phi0 = z[-1,1] - beta

nmax=10
eps = 1.0e-10
for n in range(nmax):
    z0[1] = s[1]
    z = solver(f,z0,y)
    phi1 = z[-1,1] - beta
    ds = dsfunction(phi0,phi1,s[0],s[1])
    s[0]  = s[1]
    s[1] +=  ds
    phi0 = phi1
    print 'n = {}  s1 = {} and ds = {}'.format(n,s[1],ds)
    
    if (abs(ds)<=eps):
        print 'Solution converged for eps = {} and s1 ={} and ds = {}. \n'.format(eps,s[1],ds)
        break

legends=[] # empty list to append legends as plots are generated

plot(y,z[:,0])
legends.append('theta')

plot(y,z[:,1])
legends.append('dtheta/dl')

plot(y,z[:,2])
legends.append('deflection y')



# Add the labels
legend(legends,loc='best',frameon=False) # Add the legends
ylabel('theta theta')
xlabel('y/L')
#grid(b=True, which='both', axis='both',linestyle='-')
grid(b=True, which='both', color='0.65',linestyle='-')

show()
