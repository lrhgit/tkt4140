import matplotlib
import matplotlib.pyplot as plt
import numpy as np
LNWDT=2; FNT=15
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT

def f(x):
    f=np.ones(x.shape)
    f[np.where(x>0.5)] = 0.0
    
    return f

(xmin, xmax) = (0,1)
n=20  #Number of space intervals

x = np.linspace(xmin,xmax,n+1,'float')
dx = float(xmax-xmin)/n

a=1.0
C=1.01
dt = C*dx/a

(tmin,tmax) = (0, 0.4)

m=round((tmax-tmin)/dt) #Number of temporal intervals
time=np.linspace(tmin,tmax,m)


(umin,umax)=(min(f(x)),max(f(x)))

fig = plt.figure()
ax=fig.add_subplot(111)
u0=f(x)
Curve, = ax.plot(x,u0)
ax.set_xlim([xmin,xmax])
ax.set_ylim([umin,1.4*umax])
plt.xlabel('x')
plt.ylabel('Velocity')

plt.ion()
plt.show()

nOutputInt=2
i=0 

u = np.zeros(u0.shape)
u[:]=u0[:]

ua=np.zeros(u.shape)

for t in time:
    i+=1
    
    ua[:] = f(x-a*t)

    u[1:n-2] = C*(1+C)/2*u[0:n-3] + (1-C**2)*u[1:n-2] - C*(1-C)/2*u[2:n-1]
#    u[:] = ua[:]

    if (np.mod(i,nOutputInt)==0):
        Curve.set_ydata(u)
        plt.pause(0.0001)
        plt.title('step = %3d; t = %f' % (i,t))

#Curve.set_ydata(u)
plt.pause(10)
plt.ion()
plt.close()

print 'done'


