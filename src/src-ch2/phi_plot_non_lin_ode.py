# src/src-ch2/phi_plot_non_lin_ode.py;ODEschemes.py @ git@lrhgit/tkt4140/src/src-ch2/ODEschemes.py;

from ODEschemes import euler, heun, rk4
from numpy import cos, sin
import numpy as np
from matplotlib.pyplot import *

# change some default values to make plots more readable 
LNWDT=2; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)


N=20
L = 1.0
x = np.linspace(0,L,N+1)

def f(z, t):
    zout = np.zeros_like(z)
    zout[:] = [z[1],3.0*z[0]**2/2.0]
    return zout 

beta=1.0 # Boundary value at x = L

solvers = [euler, heun, rk4] #list of solvers
solver=solvers[2] # select specific solver

smin=-45.0
smax=1.0
s_guesses = np.linspace(smin,smax,20)

# Guessed values
#s=[-5.0, 5.0]

z0=np.zeros(2)
z0[0] = 4.0
    
z = solver(f,z0,x)
phi0 = z[-1,0] - beta

nmax=10
eps = 1.0e-3
phi = []
for s in s_guesses:
    z0[1] = s
    z = solver(f,z0,x)
    phi.append(z[-1,0] - beta)
    

legends=[] # empty list to append legends as plots are generated

plot(s_guesses,phi)
ylabel('phi')
xlabel('s')
grid(b=True, which='both', color='0.65',linestyle='-')

show()
close()
