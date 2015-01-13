import numpy as np
def s(t, v0, a):
    s =[]
    for te in t:
        s.append(v0*te + 0.5*a*te**2)

    return s

def s2(t, v0, a):
    return v0*t + 0.5*a*t**2


# L=[1.0, 2.0]

t = np.linspace(0.0,2.0)
L = np.array(t).tolist()

v0=1.0
a=0.5
import time
t0=time.clock()
print s(L,v0,a)
print time.clock()-t0

t0=time.clock()
print s2(t,v0,a)
print time.clock()-t0
