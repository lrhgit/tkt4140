#-----------------------------------------------------------------------------
# The equation to be solved 
#
#       du     d  du
#       -- = k -- --
#       dt     dx dx
#
# along with boundary conditions
#
#       u(xmin,t) = U0
#       u(xmax,t) = 0
#
# and initial conditions
#
#       u(x,tmin) = 0
#-----------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import time as timelib

k  = 0.5                       #Thermal diffusivity
(xmin, xmax) = (0,1)
n  = 20; #32;                    # Number of spatial intervals
dx = (xmax-xmin)/float(n)
x = np.linspace(xmin,xmax,n+1)

r=0.6                       #Numerical Fourier number
dt=r*dx**2/k**2             #Compute timestep based on Fourier number, spatial discretization and thermal diffusivity
print 'timestep = ',dt
(tmin, tmax)=(0,0.5)


m=round((tmax-tmin)/dt)                  # Number of temporal intervals
time=np.linspace(tmin,tmax,m)
U0 = 4
print 'number of timesteps: ',m

u=np.zeros((n+1,1),float)
u[0]=U0
(umin,umax)=(min(u),max(u))


#Plot initial solution
fig = plt.figure()
ax=fig.add_subplot(111)
Curve, = ax.plot( x, u[:], '-')
ax.set_xlim([xmin,xmax])
ax.set_ylim([umin,umax])
plt.xlabel('x')
plt.ylabel('Velocity')

plt.ion()
plt.show()
nOutputInt=1
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
