# src-ch3/trapes.py;

import numpy as np
from matplotlib.pyplot import *
from numpy import sqrt as sqrt
from scipy.special import kv as besselk
from scipy.special import iv as besseli


# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)

z0 = 4*sqrt(2)
z1 = 8
g = sqrt(2)/8

k1z0, i0z1 = besselk(1, z0), besseli(0, z1)
k0z1, i1z0 = besselk(0, z1), besseli(1, z0)
k0z0, i0z0 = besselk(0, z0), besseli(0, z0)

J = k1z0*i0z1 + k0z1*i1z0 + g*(k0z0*i0z1 - k0z1*i0z0)
A = (k1z0 + g*k0z0)/J
B = (i1z0 - g*i0z0)/J
x = np.linspace(0, 1,  11)
z = sqrt(32.*(1 + x))
theta = A* besseli(0, z) + B* besselk(0, z)
dtheta = 16*(A*besseli(1, z) - B*besselk(1, z))/z

legends=[] # empty list to append legends as plots are generated

plot(x,theta,'b',) 

legends.append('$\theta$')
plot(x,dtheta,'r')
legends.append("$\theta$'")

## Add the labels
legend(legends,loc='best',frameon=False) # Add the legends
ylabel("$\theta$ / $\theta$'")
xlabel('x')

show()