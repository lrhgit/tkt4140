from ODEschemes import euler, heun, rk4
from matplotlib.pyplot import *
# Change some default values to make plots more readable on the screen
LNWDT=3; FNT=20
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT

def fblasius(y, x):
    """ODE-system for the Blasius-equation"""
    return [y[1],y[2], -y[0]*y[2]]


def dsfunction(phi0,phi1,s0,s1):
    if (abs(phi1-phi0)>0.0):   
        return    -phi1 *(s1 - s0)/float(phi1 - phi0)
    else:
        return 0.0

solvers = [euler, heun, rk4] #list of solvers
solver=solvers[0] # select specific solver


from numpy import linspace, exp, abs
xmin = 0
xmax = 5.750

N = 400  # no x-values
x = linspace(xmin, xmax, N+1)

# Guessed values
s=[0.1,0.8]


z0=np.zeros(3)
z0[2] = s[0]


beta=1.0 #Boundary value for eta=infty

## Compute phi0

u = solver(fblasius, z0, x)
phi0 = u[-1,1] - beta

nmax=10
eps = 1.0e-3
for n in range(nmax):
    z0[2] = s[1]
    u = solver(fblasius, z0, x)
    phi1 = u[-1,1] - beta
    ds = dsfunction(phi0,phi1,s[0],s[1])
    s[0]  = s[1]
    s[1]  += ds
    phi0 = phi1
    print 'n = {}  s1 = {} and ds = {}'.format(n,s[1],ds)
    
    if (abs(ds)<=eps):
        print 'Solution converged for eps = {} and s1 ={} and ds = {}. \n'.format(eps,s[1],ds)
        break


plot(u[:,1],x,u[:,2],x)
xlabel('u og u\'')
ylabel('eta')


legends=[]
legends.append('velocity')
legends.append('wall shear stress')
legend(legends,loc='best',frameon=False)
title('Solution of the Blaisus eqn with '+str(solver.func_name)+'-shoot')
show()
close() #Close the window opened by show() 








