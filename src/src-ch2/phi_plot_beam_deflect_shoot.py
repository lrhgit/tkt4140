# src/src-ch2/phi_plot_beam_shoot.py;ODEschemes.py @ git@lrhgit/tkt4140/src/src-ch2/ODEschemes.py;

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
y = np.linspace(0,L,N+1)

def f(z, t):
    """RHS for deflection of beam"""
    zout = np.zeros_like(z)
    zout[:] = [z[1],-alpha2*cos(z[0]),sin(z[0])]
    return zout 

solvers = [euler, heun, rk4] #list of solvers
solver=solvers[2] # select specific solver

alpha2 = 5.0 
beta=0.0 # Boundary value at y = L

N_guess = 30
s_guesses=np.linspace(1,5,N_guess)

z0=np.zeros(3)

phi = []
for s_guess in s_guesses:
    z0[1] = s_guess
    z = solver(f,z0,y)
    phi.append(z[-1,1] - beta)

legends=[] # empty list to append legends as plots are generated
plot(s_guesses,phi)

# Add the labels
legend(legends,loc='best',frameon=False) # Add the legends
title('alpha2 = ' + str(alpha2))
ylabel('phi')
xlabel('s')
grid(b=True, which='both')
show()
