#/python_intro/solve_simple_equation.py;solver.py @ git@lrhgit/tkt4140/allfiles/digital_compendium/python_intro/doc/src/src-bumpy/solver.py
import numpy as np
from solver import solver
from numpy import sin, pi  # for nice math

def F(t):
    # Sinusoidal bumpy road
    return A*sin(pi*t)

def s(u):
    return k*u

A = 0.25
k = 2
t = np.linspace(0, 20, 2001)
u, t = solver(I=0.1, V=0, m=2, b=0.05, s=s, F=F, t=t)

# Show u(t) as a curve plot
import matplotlib.pyplot as plt
plt.plot(t, u)
plt.show()