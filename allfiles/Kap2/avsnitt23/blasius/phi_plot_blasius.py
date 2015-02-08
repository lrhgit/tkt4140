from numpy import linspace, exp

def fblasius(y, x):
    """ODE-system for the Blasius-equation"""
    return [y[1],y[2], -y[0]*y[2]]

import odespy
from matplotlib.pyplot import *
#change some default values to make plots more readable on the screen
LNWDT=2; FNT=20
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT

solvers=[]
solvers.append(odespy.RK4(fblasius))
solvers.append(odespy.RK2(fblasius))
solvers.append(odespy.RK3(fblasius))

xmin = 0
xmax = 5.75

N = 150  # no x-values
xspan = linspace(xmin, xmax, N+1)

smin=0.1
smax=0.8

Ns=30
srange = linspace(smin,smax,Ns)

solver=solvers[2]                         
phi=np.zeros(srange.size)
beta=1
i=0

for s in srange:
    solver.set_initial_condition([0.0, 0.0, s])
    u, x = solver.solve(xspan)
    phi[i] = u[-1,1] -beta
    i+=1


figure()
legends=[]
linet=['r-',':','.','-.','--']


plot(srange,phi)
xlabel('s')
ylabel('phi')
grid(b=True, which='both', color='0.65',linestyle='-')


show()
close()







