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
h=plt.scatter(x,y,s=250)
#plt.show()

i=0

for t in np.arange(tmin+dt,tmax,dt):
    i = i+ 1
    xn = x
    yn = y

    up = u+dt*a(xn,yn)
    vp = v+dt*a(yn,xn)

    xp=x+dt*up
    yp=y+dt*vp

    x=x+dt*(u+up)/2
    y=y+dt*(v+vp)/2
    
    u=u+dt*(a(xn,yn)+a(xp,yp))/2
    v=v+dt*(a(yn,xn)+a(yp,xp))/2

    if np.mod(i,100)==0:
        plt.scatter(0,0,s=500)
        plt.scatter(x,y,s=250)
        plt.draw()
        plt.pause(.1)
        plt.clf()

        
print 'done'
