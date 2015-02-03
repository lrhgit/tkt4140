from ODEschemes import euler, heun, rk4
import numpy as np
from matplotlib.pyplot import *

# change some default values to make plots more readable 
LNWDT=5; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)


N=200
L = 1.0
y = np.linspace(0,L,N+1)

def f(z, t):
    """RHS for Couette-Posieulle flow"""
    zout = np.zeros_like(z)
    zout[:] = [z[1], -dpdx]
    return zout 

def u_a(y,dpdx):
    return y*(1.0 + dpdx*(1.0-y)/2.0);

beta=1.0 # Boundary value at y = L


# Guessed values
s=[1.0, 1.5]

z0=np.zeros(2)

dpdx_list=[-5.0, -2.5, -1.0, 0.0, 1.0,2.5, 5.0]
legends=[]

for dpdx in dpdx_list:
    phi = []
    for svalue in s:
        z0[1] = svalue
        z = rk4(f, z0, y)
        phi.append(z[-1,0] - beta)
    
    # Compute correct initial guess 
    s_star = (s[0]*phi[1]-s[1]*phi[0])/(phi[1]-phi[0])
    z0[1] = s_star
    
    # Solve the initial value problem which is a solution to the boundary value problem
    z = rk4(f, z0, y)

    plot(z[:,0],y,'-.')
    legends.append('rk4: dp='+str(dpdx))
    
    # Plot the analytical solution
    plot(u_a(y, dpdx),y,':')
    legends.append('exa: dp='+str(dpdx))

# Add the labels
legend(legends,loc='best',frameon=False) # Add the legends
xlabel('u/U0')
ylabel('y/L')
show()

