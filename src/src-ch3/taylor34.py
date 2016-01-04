# chapter3/src-ch3/taylor34.py;TRIdiagonalSolvers.py @ git@lrhgit/tkt4140/allfiles/digital_compendium/chapter3/src-ch3/TRIdiagonalSolvers.py;


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
print n
fac = (3.)*h**2

#initilizing:
x = np.linspace(h,1-h,n)
print x
ym = -20*x*(1-x) # initial guess of y. ym = np.zeros(n) will converge to y1. ym = -20*x*(1-x) for instance will converge to y2


it, itmax, dymax, relTol = 0, 15, 1., 10**-10
legends=[]
while (dymax > relTol) and (it < itmax):
    """iteration process of linearized system of equations using taylor 
    """
    plot(x,ym) # plot ym for iteration No. it
    legends.append(str(it))
    it = it + 1
    a = np.ones(n) 
    c = np.ones(n)
    b = -(np.ones(n)*2. + fac*ym)
    d = -(fac*0.5*ym**2)
    d[n-1] = d[n-1]-1
    d[0] = d[0]-4

    ym1 = tdma(a,b,c,d) # solution
    dymax  = np.max(np.abs((ym1-ym)/ym))
    ym = ym1

    print 'it = {},  dymax = {} '.format(it, dymax)
    
legend(legends,loc='best',frameon=False)
show()

ya = 4./(1+x)**2
feil = np.abs((ym1-ya)/ya)
print "\n"

for l in range(len(x)):
    print 'x = {},  y = {},  ya = {}'.format(x[l], ym1[l], ya[l])



# plotting:
legends=[] # empty list to append legends as plots are generated
# Add the labels
plot(x,ya,"r") # plot analytical solution
legends.append('y1 analytic')
plot(x,ym1,"b--") # plot numerical solution
legends.append('y2 taylor linearization')

legend(legends,loc='best',frameon=False) # Add the legends
ylabel('y')
xlabel('x')
#grid(b=True, which='both', axis='both',linestyle='-')
grid(b=True, which='both', color='0.65',linestyle='-')

show()
