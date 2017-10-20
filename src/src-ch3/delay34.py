# src-ch3/delay34.py;TRIdiagonalSolvers.py @ git@lrhgit/tkt4140/src/src-ch3/TRIdiagonalSolvers.py;


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
fac = (3./2.)*h**2
Nmax = 30 # max number of iterations

# iterative solution with Nmax iterations
#initilizing:
ym = np.zeros(n) # initial guess of y. ym = np.zeros(n) will converge to y1. ym = -20*x*(1-x) for instance will converge to y2

for m in range(Nmax):
    """iteration process of linearized system of equations using delay method"""
    # need to set a, b, c and d vector inside for loop since indices are changed in tdma solver
    a = np.ones(n) 
    c = np.ones(n)
    b = -(np.ones(n)*2. + fac*ym)
    d = np.zeros(n)
    d[n-1] = -1.
    d[0] = -4.
    
    ym1 = tdma(a,b,c,d) # solution
        
    if(m>0): 
        max_dym  = np.max(np.abs((ym1 - ym)/ym))
        print('m = {} \t  max dym = {:6.4g} '.format(m, max_dym))
        
    ym = ym1


# compute analytical solution
xv = np.linspace(h, 1 - h, n)
ya = 4./(1 + xv)**2
error = np.abs((ym1 - ya)/ya)

print('\nThe max relative error after {} iterations is: {:4.3g}'.format(Nmax,np.max(error)))


# iterative solution with max_dym stop criterion

m=0; RelTol=1.0e-8
max_dym=1.0

ym = np.zeros(n) # initial guess of y. ym = np.zeros(n) will converge to y1. ym = -20*x*(1-x) for instance will converge to y2

a = np.ones(n) 
c = np.ones(n)

while (m<Nmax and max_dym>RelTol):

    d = np.zeros(n)
    d[n-1] = -1.
    d[0] = -4.

    b = -(np.ones(n)*2. + fac*ym)
    
    ym1 = tdma(a,b,c,d) # solution

    if(m>0): 
        max_dym  = np.max(np.abs((ym1 - ym)/ym))
        print('m = {} \t  max dym = {:6.4g} '.format(m, max_dym))
        
    m+=1        
    ym = ym1
    
    
error = np.abs((ym1 - ya)/ya)
print('\nThe max relative after {} iterations is: {:4.3g}'.format(m,np.max(error)))    

# print results nicely with pandas
print_results=0
if (print_results):
    import pandas as pd 
    data=np.column_stack((xv,ym1,ya))
    data_labeled = pd.DataFrame(data, columns=['xv','ym','ya'])  
    print(data_labeled.round(3).to_string(index=False))

# plotting:
legends=[] # empty list to append legends as plots are generated
# Add the labels
plot(xv,ya,"r") # plot analytical solution
legends.append('analytic')
plot(xv,ym1,"b--") # plot numerical solution
legends.append('delayed linearization')

legend(legends,loc='best',frameon=False) # Add the legends
ylabel('y')
xlabel('x')
#grid(b=True, which='both', axis='both',linestyle='-')
grid(b=True, which='both', color='0.65',linestyle='-')

show()
