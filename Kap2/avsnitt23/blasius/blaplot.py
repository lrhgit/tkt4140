def fblasius(y, x):
    """ODE-system for the Blasius-equation"""
    return [y[1],y[2], -y[0]*y[2]]

import odespy

class Kutta4(odespy.Solver):
    """
   Kutta 4 method by LR Hellevik::

       u[n+1] = u[n] + (K1 + 3*K2 + 3*K3 + K4)/8.0 

    where::
           dt3 = dt/3.0
           K1 = dt*f(u[n], t[n])
           K2 = dt*f(u[n] + K1/3.0, t[n] + dt3)
           K3 = dt*f(u[n] - K1/3 +  K2, t[n] + 2*dt3)
           K4 = dt*f(u[n] + K1 - K2 + K3, t[n] + dt)
      """
    quick_description = "Explicit 4th-order Kutta method. LR Hellevik implementation"

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        dt3 = dt/3.0
        K1 = dt*f(u[n], t[n])
        K2 = dt*f(u[n] + K1/3.0, t[n] + dt3)
        K3 = dt*f(u[n] - K1/3 +  K2, t[n] + 2*dt3)
        K4 = dt*f(u[n] + K1 - K2 + K3, t[n] + dt)
        u_new = u[n] + (K1 + 3*K2 + 3*K3 + K4)/8.0
        return u_new


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

smin=0.1
smax=0.8

Ns=30
srange = linspace(smin,smax,Ns)

from matplotlib.pyplot import *
#change some default values to make plots more readable on the screen
LNWDT=5; FNT=25
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT
figure()
legends=[]
linet=['r-',':','.','-.','--']

solver=solvers[0]
phi=np.zeros(srange.size)
beta=1
i=0

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

# plot(x,u[:,0],x,u[:,1],x,u[:,2])
# xlabel('x')
# ylabel('u')

# legend(legends)

show()








