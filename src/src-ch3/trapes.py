# coding: utf-8
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

def trapes_analytical(h):
    
    N=int(round(1.0/h))      # number of unknowns, assuming the RHS boundary value is known
    x=np.arange(N+1)*h   
#    x = np.linspace(0, 1,  N)
    z0 = 4*sqrt(2)
    z1 = 8
    g = sqrt(2)/8
    
    k1z0, i0z1 = besselk(1, z0), besseli(0, z1)
    k0z1, i1z0 = besselk(0, z1), besseli(1, z0)
    k0z0, i0z0 = besselk(0, z0), besseli(0, z0)
    
    J = k1z0*i0z1 + k0z1*i1z0 + g*(k0z0*i0z1 - k0z1*i0z0)
    A = (k1z0 + g*k0z0)/J
    B = (i1z0 - g*i0z0)/J
    
    z = sqrt(32.*(1 + x))
    theta = A* besseli(0, z) + B* besselk(0, z)
    dtheta = 16*(A*besseli(1, z) - B*besselk(1, z))/z
    return x, theta, dtheta 


if __name__ == '__main__':
    legends=[] # empty list to append legends as plots are generated
    
    h=0.2
    x, theta, dtheta = trapes_analytical(h)

    plot(x,theta,'b',) 
    legends.append(r'$\theta$')
    plot(x,dtheta,'r')
    legends.append(r'$\theta^{\prime}$')
    
    ## Add the labels
    legend(legends,loc='best',frameon=False) # Add the legends
    ylabel(r'$\theta$ and $\theta^{\prime}$')
    xlabel('x')
    
    
    
    show()
    
    close()
