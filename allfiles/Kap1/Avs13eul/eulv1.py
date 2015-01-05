#from numpy import *
import numpy as np
import matplotlib.pyplot as plt

g    = 9.81
alfa = 7.0e-2
dt   = 0.051
Tend = 1
N    = np.rint(Tend/dt)
N = N.astype(np.int64)

Tend = N*dt #True end

def va(t):
    k1=np.sqrt(g/alfa)
    k2=np.sqrt(alfa*g)
    return k1*np.tanh(k2*t)

def vt(t):
   return g*t*(1-alfa*np.power(t,2)+2*alfa**2*g**2*np.power(t,4)/15);

t=np.linspace(0,Tend,N+1)
v=np.zeros((N+1))

for i in xrange(0,N):
    v[i+1] = v[i]+dt*(g-alfa*v[i]*v[i])

plt.plot(t,v,label='Euler')
plt.plot(t,va(t),label='Analytical')
plt.plot(t,vt(t),label='Taylor')
plt.legend(loc='upper left')
plt.show()
