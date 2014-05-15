from matplotlib.pyplot import *
# Change some default values to make plots more readable on the screen
LNWDT=3; FNT=20
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT
from numpy import linspace, exp, abs

# Import odespy a module with many ode-solvers
import odespy

legends=[]
linet=['r-',':','.','-.','--']

def fbeam(T, x):
    """ODE-system for the Blasius-equation"""
    h = 0.05
    Tinf = 200.0
    sigma = 2.7E-9
        
    return [T[1],-h*(Tinf-T[0])-sigma*(Tinf**4-T[0]**4)]



def dsfunction(phi0,phi1,s0,s1):
    if (abs(phi1-phi0)>0.0):   
        return    -phi1 *(s1 - s0)/float(phi1 - phi0)
    else:
        return 0.0


solvers=[]
solvers.append(odespy.RK4(fbeam))
solvers.append(odespy.RK2(fbeam))
solvers.append(odespy.RK3(fbeam))

solver=solvers[0]                         


xmin = 0.0
xmax = 3.0

N = 100  # no x-values
xspan = linspace(xmin, xmax, N+1)

# From the blaplot.py we have two initial guesses
s0 = 10
s1 = -10

pout=50.0 #Boundary value for x=xmax
T0 =300
T10=400
## Compute phi0

solver.set_initial_condition([T0,s0])
T, x = solver.solve(xspan)
phi0 = T[-1,0] - T10



nmax=10
eps = 1.0e-3
for n in range(nmax):
    solver.set_initial_condition([T0,s1])
    T, x = solver.solve(xspan)
    phi1 = T[-1,0] - T10
    ds = dsfunction(phi0,phi1,s0,s1)
    s0  = s1
    s1  = s1 + ds
    phi0 = phi1
    print 'n = {}  s1 = {} and ds = {} and T(end) - T10 = {}'.format(n,s1,ds, T[-1,0]-T10)
    
    if (abs(ds)<=eps):
        print 'Solution converged for eps = {} and s1 ={} , ds = {}  and T(end) - T10 = {}. \n'.format(eps,s1,ds, T[-1,0]-T10)  
        break


plot(x,T[:,0])
xlabel('x')
ylabel('T')
legends.append(str(solver))
legend(legends)
show()
close() #Close the window opened by show() 
print 'the end'







