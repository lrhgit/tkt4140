<<<<<<< HEAD:allfiles/Kap2/avsnitt23/blasius/blasius.py
mu = 3.0
def f(u, t):
    """2x2 system for a van der Pool oscillator."""
    return [u[1], mu*(1. - u[0]*u[0])*u[1] - u[0]]

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
solvers.append(odespy.RK4(f))
solvers.append(odespy.RK2(f))
solvers.append(odespy.RK3(f))
solvers.append(Kutta4(f))

from numpy import linspace, exp
T = 30  # end of simulation
N = 150  # no of time steps
time = linspace(0, T, N+1)


from matplotlib.pyplot import *
#change some default values to make plots more readable on the screen
LNWDT=5; FNT=25
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT
figure()
legends=[]
linet=['r-',':','.','-.','--']

i=0
for solver in solvers:
    solver.set_initial_condition([2.0, 0.0])
    u, t = solver.solve(time)
    plot(t,u[:,0],linet[i])
    legends.append(str(solver))
    i+=1

xlabel('time')
legend(legends)

show()








=======
mu = 3.0
def f(u, t):
    """2x2 system for a van der Pool oscillator."""
    return [u[1], mu*(1. - u[0]*u[0])*u[1] - u[0]]

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
solvers.append(odespy.RK4(f))
solvers.append(odespy.RK2(f))
solvers.append(odespy.RK3(f))
solvers.append(Kutta4(f))

from numpy import linspace, exp
T = 30  # end of simulation
N = 5000  # no of time steps
time = linspace(0, T, N+1)


from matplotlib.pyplot import *
#change some default values to make plots more readable on the screen
LNWDT=3; FNT=20
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT
figure()
legends=[]
linet=['r-',':','.','-.','--']

i=0
for solver in solvers:
    solver.set_initial_condition([2.0, 0.0])
    u, t = solver.solve(time)
    plot(t,u[:,0],linet[i])
    legends.append(str(solver))
    i+=1

xlabel('time')
legend(legends)

show()








>>>>>>> master:allfiles/Kap2/avsnitt23/blasius/vdPool.py
