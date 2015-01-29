#from numpy import *
import numpy as np
import matplotlib.pyplot as plt

MG=1
x=1
y=0
u=0
v=.7

def a(x,y):
    r=np.sqrt(x**2+y**2)
    return -MG*x/r**3

tmin=0
tmax= 100
dt=0.01
plt.scatter(0,0,s=500)
plt.scatter(x,y,s=250)
plt.show()

for t in arange(tmin+dt,tmax,dt):
    xn = x
    yn = y

    up = u+dt*a(xn,yn)
    vp = v+dt*a(yn,xn)

    x=x+dt*(u+up)/2
    y=y+dt*(v+vp)/2
    
    u=u+dt*(a(xn,yn)+a(xp,yp))/2
    v=v+dt*(a(yn,xn)+a(yp,xp))/2

    
print 'done'
