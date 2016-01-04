from ODEschemes import euler, heun, rk4
import numpy as np

N=20
L = 1.0
y = np.linspace(0,L,N+1)

def f(z, t, dpdx=0.0):
    """RHS for Couette-Posieulle flow"""
    zout = np.zeros_like(z)
    zout[:] = [z[1], -dpdx]
    return zout 

s0=1.0
z0=np.zeros(2)
z0[1] = s0


z = rk4(f, z0, y)
phi0 =z[]