from matplotlib.pyplot import *
#change some default values to make plots more readable on the screen
LNWDT=5; FNT=25
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT
figure()
legends=[]
linet=['r-',':','.','-.','--']

def fblasius(y, x):
    """ODE-system for the Blasius-equation"""
    return [y[1],y[2], -y[0]*y[2]]

import odespy

solvers=[]
solvers.append(odespy.RK4(fblasius))
solvers.append(odespy.RK2(fblasius))
solvers.append(odespy.RK3(fblasius))
solvers.append(Kutta4(fblasius))

from numpy import linspace, exp
xmin = 0
xmax = 5.75

N = 150  # no x-values
xspan = linspace(xmin, xmax, N+1)

# From the blaplot.py we have two initial guesses

s0 = 0.3
s1 = 0.5

solver=solvers[3]                         
phi=np.zeros(srange.size)
beta=1
i=0

## Compute phi0
solver.set_initial_condition([0.0, 0.0, s0])
u, x = solver.solve(xspan)
phi0 = u[-1,1] -beta

nmax=10

for 

for s in srange:
    solver.set_initial_condition([0.0, 0.0, s])
    u, x = solver.solve(xspan)
    phi[i] = u[-1,1] -beta
    i+=1

# i=0
# for solver in solvers:
#     solver.set_initial_condition([2.0, 0.0])
#     u, t = solver.solve(time)
#     plot(t,u[:,0],linet[i])
#     legends.append(str(solver))
#     i+=1

plot(srange,phi)
xlabel('s')
ylabel('phi')


show()








