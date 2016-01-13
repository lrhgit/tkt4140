# src-ch3/delta34.py;ODEschemes.py @ git@lrhgit/tkt4140/src/src-ch3/ODEschemes.py;


import numpy as np
from matplotlib.pyplot import *
from TRIdiagonalSolvers import tdma
from TRIdiagonalSolvers import tripiv
# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)

h = 0.05 # steplength
n = int(round(1./h)-1) # number of equations
fac = (3.)*h**2
nitr = 6 # number of iterations

#initilizing:

ym = np.zeros(n) # initial guess of y. ym = np.zeros(n) will converge to y1. ym = -20*x*(1-x) for instance will converge to y2

it, itmax, dymax, relTol = 0, 10, 1., 10**-5 # numerical tollerance limits

while (dymax > relTol) and (it < itmax):
    """iteration process of linearized system of equations"""
    it = it + 1
    a = np.ones(n) 
    c = np.ones(n)
    b = -(np.ones(n)*2. + fac*ym)
    d = (np.ones(n)*fac*0.5*ym**2)
    for j in range(1,n-1):
        d[j] = d[j] - (ym[j+1] - 2*ym[j] +ym[j-1])
    
    d[n-1] = d[n-1]-(1 - 2*ym[n-1] + ym[n-2])
    d[0] = d[0] -(ym[1] - 2*ym[0] + 4)

    dy = tdma(a,b,c,d) # solution
    ym = ym + dy
    dymax  = np.max(np.abs((dy)/ym))
    

    print 'it = {},  dymax = {} '.format(it, dymax)


xv = np.linspace(h,1-h,n)
ya = 4./(1+xv)**2

print "\n"

for l in range(len(xv)):
    print 'x = {},  y = {},  ya = {}'.format(xv[l], ym[l], ya[l])


legends=[] # empty list to append legends as plots are generated

plot(xv,ya,"r") # plot analytical solution
legends.append('y1 analytic')
plot(xv,ym,"b--") # plot numerical solution
legends.append('delta linearization')
# Add the labels
legend(legends,loc='best',frameon=False) # Add the legends
ylabel('y')
xlabel('x')
#grid(b=True, which='both', axis='both',linestyle='-')
grid(b=True, which='both', color='0.65',linestyle='-')

show()
