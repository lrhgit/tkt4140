#-----------------------------------------------------------------------------
# Inspired by 
# Jonathan Senning Gordon College but...seriously modified by Leif Rune Hellevik march 2014
#
# Use FTCS (forward time, centered space) scheme to solve the heat equation
# in a thin rod.
#
# The equation solved is
#
#       du     d  du
#       -- = k -- --
#       dt     dx dx
#
# along with boundary conditions
#
#       u(xmin,t) = a(t)
#       u(xmax,t) = b(t)
#
# and initial conditions
#
#       u(x,tmin) = f(x)
#
#-----------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import time as timelib


def f(x):
    return np.sin(np.pi*x) + 3.0*np.sin(2*np.pi*x)

def f2(x):
    return np.sin(np.pi*x)*(1+6.0*np.cos(np.pi*x))

def uexact(x,t):
    return np.sin(np.pi*x)*np.exp(-np.pi**2*t) +3.0*np.sin(2*np.pi*x)*np.exp(-4*np.pi**2*t) 

def uexact2(x,t,U0,L):
    n=2
    omegas=np.pi*(2.0*np.linspace(1,n,n)+1)/2/L
    
    for omegan in omegas:
        bn=(alpha*omegan)**2
        u += np.cos(omegan*x)*np.exp(-bn*t) 
        
    return u*U0

k  = 0.25                       #Thermal diffusivity
(xmin, xmax) = (0,1)
n  = 20; #32;                    # Number of spatial intervals
dx = (xmax-xmin)/float(n)
x = np.linspace(xmin,xmax,n+1)

r=0.5                          #Numerical Fourier number
dt=r*dx**2/k**2                   #Compute timestep based on Fourier number, spatial discretization and thermal diffusivity
print 'timestep = ',dt
(tmin, tmax)=(0,2.0)

#m = 128; #4096;                 # Number of temporal intervals
m=round((tmax-tmin)/dt)
time=np.linspace(tmin,tmax,m)

print 'number of timesteps: ',m
u=np.zeros((n+1,1),float)
u[0]=1
u=uexact(x,0.0)

(umin,umax)=(min(u),max(u))
#Plot initial solution
fig = plt.figure()
ax=fig.add_subplot(111)
Curve, = ax.plot( x, u[:], '-')
ax.set_xlim([xmin,xmax])
ax.set_ylim([umin,umax])
plt.xlabel('x')
plt.ylabel('Temperature')

plt.ion()
plt.show()
nOutputInt=10
i = 0

for t in time:
    i+=1
    u[1:-1] =  r*(u[0:-2]+ u[2:]) + (1.0-2.0*r)*u[1:-1]
    
    if (np.mod(i,nOutputInt)==0):
        Curve.set_ydata(u)
        plt.pause(.05)
        plt.title( 'step = %3d; t = %f' % (i,t ) )
        
Curve.set_ydata(u)

plt.pause(2)
plt.ion()
plt.close()

print 'done'