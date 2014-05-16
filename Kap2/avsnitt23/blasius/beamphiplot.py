from matplotlib.pyplot import *
import scipy as sc
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import time


# Change some default values to make plots more readable on the screen
LNWDT=3; FNT=20
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT
from numpy import linspace, exp, abs
import numpy.linalg as la

# Import odespy a module with many ode-solvers
import odespy

legends=[]
linet=['r-',':','.','-.','--']
Tinf = 200.0
h = 0.05
sigma = 2.7E-9
#sigma = 0.0

def fbeam(T, x):
    """ODE-system for the Blasius-equation"""
        
    return [T[1],-h*(Tinf-T[0])-sigma*(Tinf**4-T[0]**4)]



def dsfunction(phi0,phi1,s0,s1):
    if (abs(phi1-phi0)>0.0):   
        return    -phi1 *(s1 - s0)/float(phi1 - phi0)
    else:
        print 'we have problems'
        return 0.0


solvers=[]
solvers.append(odespy.RK4(fbeam))
solvers.append(odespy.RK2(fbeam))
solvers.append(odespy.RK3(fbeam))

solver=solvers[0]                         


xmin = 0.0
xmax =6.0

N = 30  # no x-values, including the boundaries
xspan, dx = linspace(xmin, xmax, N+1,retstep=True)

# From the blaplot.py we have two initial guesses
s0 = -35
s1 = -30
Ns = 50
srange = linspace(s0,s1,Ns)
phi=np.zeros(Ns)

#Boundary value for x=xmax
T0 =300
T10=500
i = 0
for s in srange:
	solver.set_initial_condition([T0,s])
	T, x = solver.solve(xspan)
	phi[i] = T[-1,0] - T10
	i+=1



#Plot the solutions
plot(srange,phi)

xlabel('s')
ylabel('p')
show()
close() #Close the window opened by show() 
print 'the end'







