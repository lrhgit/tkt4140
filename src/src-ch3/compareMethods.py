# src-ch3/compareMethods.py;TRIdiagonalSolvers.py @ git@lrhgit/tkt4140/src/src-ch3/TRIdiagonalSolvers.py;

import numpy as np
from matplotlib.pyplot import *
from TRIdiagonalSolvers import tdma
from TRIdiagonalSolvers import tripiv
"""A more compact script comparing solutions using three different 
    linearization methods: delay, taylor and deltavalues
"""

# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)



def taylor(ym, h, n):
    """Method caluclating input diagonals a,b,c(sub, main and top), and rhs vecor d
        of the tridiagonal system of equations based on taylor linearization
    
        Args:
            ym(array): solutian array from last iteration
            h(float): step length
            n(integer): number of equations/unknowns

    
        Returns:
            a(array): sub diagonal
            b(array): main diagonal
            c(array): top diagonal
            d(array): rhs
    """
    fac = (3.)*h**2 # nonlinear coefficient
    
    a = np.ones(n) 
    c = np.ones(n)
    b = -(np.ones(n)*2. + fac*ym)
    d = -(fac*0.5*ym**2)
    d[n-1] = d[n-1] - 1
    d[0] = d[0] - 4
    
    return [a, b ,c ,d]

def delay(ym, h, n):
    """Method caluclating input diagonals a,b,c(sub, main and top), and rhs vecor d
        of the tridiagonal system of equations based on delay linearization
    
        Args:
            ym(array): solutian array from last iteration
            h(float): step length
            n(integer): number of equations/unknowns

    
        Returns:
            a(array): sub diagonal
            b(array): main diagonal
            c(array): top diagonal
            d(array): rhs
    """
    fac = (3./2)*h**2 # nonlinear coefficient
    
    a = np.ones(n) 
    c = np.ones(n)
    b = -(np.ones(n)*2. + fac*ym)
    d = np.zeros(n)
    d[n-1] = -1.
    d[0] = -4.
    
    return [a, b ,c ,d]

def delta(ym, h, n):
    """Method caluclating input diagonals a,b,c(sub, main and top), and rhs vecor d
        of the tridiagonal system of equations based on delta value linearization
    
        Args:
            ym(array): solutian array from last iteration
            h(float): step length
            n(integer): number of equations/unknowns

    
        Returns:
            a(array): sub diagonal
            b(array): main diagonal
            c(array): top diagonal
            d(array): rhs
    """
    fac = (3.)*h**2 # nonlinear coefficient
    
    a = np.ones(n) 
    c = np.ones(n)
    b = -(np.ones(n)*2. + fac*ym)
    d = (np.ones(n)*fac*0.5*ym**2)
    for j in range(1,n-1):
        d[j] = d[j] - (ym[j+1] - 2*ym[j] +ym[j-1])
    
    d[n-1] = d[n-1]-(1 - 2*ym[n-1] + ym[n-2])
    d[0] = d[0] -(ym[1] - 2*ym[0] + 4)
    
    return [a, b ,c ,d]

"""main program:"""

#initilizing:
h = 0.05 # steplength
n = int(round(1./h)-1) # number of equations
print n
x = np.linspace(h,1-h,n)
print x
ym = -20*x*(1-x) # initial guess of y. ym = np.zeros(n) will converge to y1. ym = -20*x*(1-x) for instance will converge to y2

it, itmax, dymax, relTol = 0, 15, 1., 10**-10 # iteration parameters
methods = [taylor, delay, delta] #different linearization methods to use
solutions = {"taylor":None, "delay":None, "delta":None} # dictionary containing solutions
loop = 0

for method in methods:
    methodname = solutions.keys()[loop] # string variable equal to the name of the method
    print methodname
    figure() # new figure
    legends=[] # empty legend list
    while (dymax > relTol) and (it < itmax):
        """iteration process of linearized system of equations
        """
        
        it = it + 1
        [a, b, c, d] = method(ym, h, n) # get updated matrix and rhs arrays
        if methodname != "delta":
            """taylor and delay method both use y rather than dy
            """
            #plot ym from last iteration
            plot(x,ym)
            legends.append(str(it))
                
            ym1 = tdma(a,b,c,d) # solution
            dymax  = np.max(np.abs((ym1-ym)/ym)) # tollerance limit
            ym = ym1
            
        else:
            """delta method use dy rather than y
            """
            #plot ym from last iteration
            plot(x,ym)
            legends.append(str(it))
            
            dy = tdma(a,b,c,d) # solution
            ym = ym + dy 
            dymax  = np.max(np.abs((dy)/ym)) # tollerance limit
    
        print 'it = {},  dymax = {} '.format(it, dymax)
    
    solutions[methodname] = ym # update solutions dict with iterated solution
    #legend, title and label
    legend(legends,loc='best',frameon=False) # Add the legends
    title(methodname+ " method")
    ylabel('y')
    xlabel('x')
    grid(b=True, which='both', color='0.65',linestyle='-')
    legends=[]
    ym = -20*x*(1-x) # reset initial guess of y for next method
    it = 0 # reset number of iterations for next method
    dymax = 1 # reset dymax for next method
    loop = loop + 1
    print"\n"


xv = np.linspace(h,1-h,n)
ya = 4./(1+xv)**2

ytaylor = solutions["taylor"]
ydelay = solutions["delay"]
ydelta = solutions["delta"]


figure() # new figure
# plotting:
legends=[] # empty list to append legends as plots are generated
# Add the labels
plot(xv,ytaylor,"b") 
legends.append('taylor')
plot(xv,ydelay,"g") 
legends.append('delay')
plot(xv,ydelta,"y--") 
legends.append('delta')
#plot(xv,ya,"r") # plot analytical solution
#legends.append('analytic')
legend(legends,loc='best',frameon=False) # Add the legends
ylabel('y')
xlabel('x')
#grid(b=True, which='both', axis='both',linestyle='-')
grid(b=True, which='both', color='0.65',linestyle='-')

show()
