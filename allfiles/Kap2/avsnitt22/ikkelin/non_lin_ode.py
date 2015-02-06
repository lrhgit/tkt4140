from ODEschemes import euler, heun, rk4
from numpy import cos, sin
import numpy as np
from matplotlib.pyplot import *

# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)


N=20
L = 1.0
x = np.linspace(0,L,N+1)

def dsfunction(phi0,phi1,s0,s1):
    if (abs(phi1-phi0)>0.0):   
        return    -phi1 *(s1 - s0)/float(phi1 - phi0)
    else:
        return 0.0


def f(z, t):
    zout = np.zeros_like(z)
    zout[:] = [z[1],3.0*z[0]**2/2.0]
    return zout 


def y_analytical(x):
    return 4.0/(1.0+x)**2
  

beta=1.0 # Boundary value at x = L

solvers = [euler, heun, rk4] #list of solvers
solver=solvers[2] # select specific solver

# Guessed values
s=[-5.0, 5.0]

z0=np.zeros(2)
z0[0] = 4.0
z0[1] = s[0]

z = solver(f,z0,x)
phi0 = z[-1,0] - beta

nmax=10
eps = 1.0e-3
for n in range(nmax):
    z0[1] = s[1]
    z = solver(f,z0,x)
    phi1 = z[-1,0] - beta
    ds = dsfunction(phi0,phi1,s[0],s[1])
    s[0]  = s[1]
    s[1] +=  ds
    phi0 = phi1
    print 'n = {}  s1 = {} and ds = {}'.format(n,s[1],ds)
    
    if (abs(ds)<=eps):
        print 'Solution converged for eps = {} and s1 ={} and ds = {}. \n'.format(eps,s[1],ds)
        break

legends=[] # empty list to append legends as plots are generated

plot(x,z[:,0])
legends.append('y')

plot(x,y_analytical(x),':')
legends.append('analytical')

# Add the labels
legend(legends,loc='best',frameon=False) # Add the legends
ylabel('y')
xlabel('x/L')
show()
