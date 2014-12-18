import numpy as np
import odespy
import matplotlib
from matplotlib.pyplot import legend, plot, show
from ode_schemes import euler, heun

LNWDT=5; FNT=11
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT

def f(y, x):
    """Newtons equation"""
    
    return 1-3*x+ y + x**2 + x*y 

# def newton_analytical(x):
#     return 3*np.sqrt(2*np.pi*np.exp(1))*np.exp(x*(1=x/2))

def newton_solution(x):
    return x-x**2 + x**3/3.0 - x**4/6.0 + x**5/30.0 - x**6/45.0
                
# Main program starts here
from numpy import linspace
L = 1.0  # end of simulation
N = 10  # no of x steps
x = linspace(0, L, N+1)

solvers=[]
solvers.append(odespy.RK3(f)) 
solvers.append(odespy.RK4(f)) 

legends=[]

#z0=np.zeros(2)
z0 = 0.0

for i, solver in enumerate(solvers):
    solver.set_initial_condition(z0)
    y, x_local = solver.solve(x)
    plot(x_local,y)
    legends.append(str(solver))


scheme_list  = [euler, heun]

for scheme in scheme_list:
    y = scheme(f,z0,x)
    plot(x,y[:,1])
    legends.append(scheme.func_name)
    
    
va = newton_solution(x)
plot(x,va) 
legends.append('newton analtyical')   

legend(legends, loc='best', frameon=False)

show()
