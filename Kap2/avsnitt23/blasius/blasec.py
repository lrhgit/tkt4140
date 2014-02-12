
from matplotlib.pyplot import *
#change some default values to make plots more readable on the screen
LNWDT=3; FNT=20
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT
import odespy

legends=[]
linet=['r-',':','.','-.','--']

def fblasius(y, x):
    """ODE-system for the Blasius-equation"""
    return [y[1],y[2], -y[0]*y[2]]


def dsfunction(phi0,phi1,s0,s1):
    if (abs(phi1-phi0)>0.0):   
        return    -phi1 *(s1 - s0)/float(phi1 - phi0)
    else:
        return 0.0


solvers=[]
solvers.append(odespy.RK4(fblasius))
solvers.append(odespy.RK2(fblasius))
solvers.append(odespy.RK3(fblasius))

from numpy import linspace, exp, abs
xmin = 0
xmax = 5.75

N = 50  # no x-values
xspan = linspace(xmin, xmax, N+1)

# From the blaplot.py we have two initial guesses

s0 = 0.1
s1 = 0.8

solver=solvers[1]                         

beta=1.0
i=0


## Compute phi0
solver.set_initial_condition([0.0, 0.0, s0])
u, x = solver.solve(xspan)
phi0 = u[-1,1] - beta

nmax=10
eps = 1.0e-3
for n in range(nmax):
    solver.set_initial_condition([0.0, 0.0, s1])
    u, x = solver.solve(xspan)
    phi1 = u[-1,1] - beta
    ds = dsfunction(phi0,phi1,s0,s1)
    s0  = s1
    s1  = s1 + ds
    phi0 = phi1
    print 'n = {}  s1 = {} and ds = {}'.format(n,s1,ds)
    
    if (abs(ds)<=eps):
        print 'Solution converged for eps = {} and {}. \n'.format(eps,dsfunction(phi0,phi1,s0,s1))
        break


plot(u[:,1],x,u[:,2],x)
xlabel('eta')
ylabel('u og u\'')
legends.append(str(solver))
legend(legends)
show()
close() #Close the window opened by show() 
print 'the end'








