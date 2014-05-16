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
xmax =5.0

N = 30  # no x-values, including the boundaries
xspan, dx = linspace(xmin, xmax, N+1,retstep=True)

# From the blaplot.py we have two initial guesses
s0 = -35
s1 = -30

#Boundary value for x=xmax
T0 =300
T10=500
## Compute phi0

solver.set_initial_condition([T0,s0])
T, x = solver.solve(xspan)
phi0 = T[-1,0] - T10



nmax=20
eps = 1.0e-6
for n in range(nmax):
    solver.set_initial_condition([T0,s1])
    T, x = solver.solve(xspan)
    phi1 = T[-1,0] - T10
    ds = dsfunction(phi0,phi1,s0,s1)
    s0  = s1
    s1  = s1 + ds
    phi0 = phi1
    print 'n = {}  s1 = {} and ds = {} and T(end) - T10 = {}'.format(n,s1,ds, T[-1,0]-T10)
    
    if (abs(ds)<=eps and abs(T[-1,0]-T10)/T10 <= eps):
        print 'Solution converged for eps = {} and s1 ={} , ds = {}  and T(end) - T10 = {}. \n'.format(eps,s1,ds, T[-1,0]-T10)  
        break

#Crete rhs array
d=np.zeros(N-1)

#Crete subiteration array for Temperature
Tm=np.ones(N-1)*T0


#Create matrix for sparse solver
diagonals=np.zeros((3,N-1))
diagonals[0,:]= 1                       #all elts in first row is set to 1
diagonals[2,:]= 1 



#Solve linear problems
for n in range(nmax):
    diagonals[1,:]= -(2+dx**2*h + 3.0*dx**2*sigma*Tm[:]**3)
    As = sc.sparse.spdiags(diagonals, [-1,0,1], N-1, N-1,format='csc') #sparse matrix instance
    d[:]=-dx**2*(h*Tinf + sigma*Tinf**4 + 2*sigma*Tm**4)
    d[0]-=T0
    d[-1]-=T10
    Ti = sc.sparse.linalg.spsolve(As,d) #theta=sc.linalg.solve_triangular(A,d)
    if (la.norm(Ti-Tm)<eps):
        print 'Solution direct nonlinear solver after n={} iterations'.format(n)
        break
    Tm = Ti
      
#Plot the solutions
plot(x,T[:,0],xspan[1:-1],Ti[:],':')

xlabel('x')
ylabel('T')
legends.append(str(solver))
legends.append('Direct')
legend(legends)
show()
close() #Close the window opened by show() 
print 'the end'







