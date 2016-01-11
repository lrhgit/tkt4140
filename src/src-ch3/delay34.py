# src/src-ch3/delay34.py;TRIdiagonalSolvers.py @ git@lrhgit/tkt4140/src/src-ch3/TRIdiagonalSolvers.py;


import numpy as np
from matplotlib.pyplot import *
from TRIdiagonalSolvers import tdma
from TRIdiagonalSolvers import tripiv
# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)

h = 0.05 # steplength
n = round(1./h)-1 # number of equations
fac = (3./2.)*h**2
nitr = 12 # number of iterations

#initilizing:
ym = np.zeros(n) # initial guess of y. ym = np.zeros(n) will converge to y1. ym = -20*x*(1-x) for instance will converge to y2


for k in range(nitr):
    """iteration process of linearized system of equations using delay method"""
    # need to set a, b, c and d vector inside for loop since indices are changed in tdma solver
    a = np.ones(n) 
    c = np.ones(n)
    b = -(np.ones(n)*2. + fac*ym)
    d = np.zeros(n)
    d[n-1] = -1.
    d[0] = -4.
    
    ym1 = tdma(a,b,c,d) # solution
    dymax  = np.max(np.abs((ym1 - ym)/ym))
    ym = ym1

    print 'k = {},  dymax = {} '.format(k, dymax)

xv = np.linspace(h, 1 - h, n)
ya = 4./(1 + xv)**2
feil = np.abs((ym1 - ya)/ya)
print "\n"

for l in range(len(xv)):
    print 'x = {},  y = {},  ya = {}'.format(xv[l], ym1[l], ya[l])

# plotting:
legends=[] # empty list to append legends as plots are generated
# Add the labels
plot(xv,ya,"r") # plot analytical solution
legends.append('analytic')
plot(xv,ym1,"b--") # plot numerical solution
legends.append('delay linearization')

legend(legends,loc='best',frameon=False) # Add the legends
ylabel('y')
xlabel('x')
#grid(b=True, which='both', axis='both',linestyle='-')
grid(b=True, which='both', color='0.65',linestyle='-')

show()
