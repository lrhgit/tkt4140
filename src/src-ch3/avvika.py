# src-ch3/avvika.py;TRIdiagonalSolvers.py @ git@lrhgit/tkt4140/src/src-ch3/TRIdiagonalSolvers.py;


import numpy as np
from matplotlib.pyplot import *
from TRIdiagonalSolvers import tdma
from TRIdiagonalSolvers import tripiv
# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)

# numerical parameters
h = 0.05 # steplength. 
n = round(1./h)-1 # number of equations. h should be set so that n is an integer
fac = (3.)*h**2

#initilizing:
x = np.linspace(h,1-h,n)
ym = -20*x*(1-x) # initial guess of y. ym = np.zeros(n) will converge to y1. ym = -20*x*(1-x) for instance will converge to y2

itmax = 15

for it in range(itmax):
    """iteration process of linearized system of equations"""
    # need to set a, b, c and d vector inside for loop since indices are changed in tdma solver
    a = np.ones(n) 
    c = np.ones(n)
    b = -(np.ones(n)*2. + fac*ym)
    d = -(np.ones(n)*fac*0.5*ym**2)
    d[n-1] = -1.
    d[0] = -4.
    
    ym1 = tdma(a,b,c,d) # solution
    dy = np.abs(ym1-ym)
    ta1 = np.max(dy) # abs stop critertia 1
    ta2 = np.sum(dy)/n # abs stop critertia 2
    ta3 = np.sqrt(np.dot(dy,dy))/n # abs stop critertia 3
    
    ym = ym1

    print 'it = {},  ta2 = {},  ta3 = {} '.format(it, ta2, ta3)


xv = np.linspace(h,1-h,n)
ya = 4./(1+xv)**2

# plotting:
legends=[] # empty list to append legends as plots are generated
# Add the labels
plot(xv,ym1,"b") 
legends.append('y2 numerical')
#plot(xv,ya,"r") 
#legends.append('y1 analytic')
legend(legends,loc='best',frameon=False) # Add the legends
ylabel('y')
xlabel('x')
#grid(b=True, which='both', axis='both',linestyle='-')
grid(b=True, which='both', color='0.65',linestyle='-')

show()
