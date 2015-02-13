from ODEschemes import euler, heun, rk4
from matplotlib.pyplot import *
# Change some default values to make plots more readable on the screen
LNWDT=3; FNT=20
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT
from numpy import linspace

def fblasius(y, x):
    """ODE-system for the Blasius-equation"""
    return [y[1],y[2],-y[0]*y[2]]

solvers = [euler, heun, rk4] #list of solvers
solver=solvers[2] # select specific solver

xmin = 0
xmax = 5.75

N = 50  # no x-values
x = linspace(xmin, xmax, N+1)

# Guessed values
#s=[0.1,0.8]
s_guesses=np.linspace(0.01,5.0)

z0=np.zeros(3)

beta=1.0 #Boundary value for eta=infty

phi = []
for s_guess in s_guesses:
    z0[2] = s_guess
    u = solver(fblasius, z0, x)
    phi.append(u[-1,1]-beta)

plot(s_guesses,phi)

title('Phi-function for the Blasius equation')
ylabel('phi')
xlabel('s')
grid(b=True, which='both')
show()




