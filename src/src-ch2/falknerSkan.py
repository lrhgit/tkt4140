# chapter2/src-ch2/falknerSkan.py;ODEschemes.py @ git@lrhgit/tkt4140/allfiles/digital_compendium/chapter2/src-ch2/ODEschemes.py;

from ODEschemes import euler, heun, rk4
from numpy import cos, sin
import numpy as np
from matplotlib.pyplot import *

# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)

def svalue(Beta):
    """method that calculates etainf (length of computational domain), and appropriate starting values of S to start iterations.
        See Appendix C.3
    
    Args:
        Beta(float): angle
        
    Returns:
        S0(float): first guess of S
        S1(float): second guess of S
        etainf(float): approximation of eta infinity (length of computational domain)
    
    Raises:
        ValueError: if Beta is out of bounds
    """
    
    if Beta >= Betasep and Beta <= 1.0:
        s0 = (1.27*(0.2+Beta))**0.56
        s1 = (1.23*(Beta-Betasep))**0.54
        
    elif Beta > 1.0 and Beta <= 1.999:
        s0 = 1.23*Beta**0,454
        s1 = -0.0693*Beta**2 + 0.661*Beta+0.642
        
    else: print "Error Beta out of value"
    
    if Beta >= Betasep and Beta <= 0.0:
        etainf= 36.76*Beta**2 + 2.0*Beta + 5.87
    
    elif Beta > 0.0 and Beta <= 1.0:
        etainf = 0.633*Beta**2 -1.68*Beta + 5.76
        
    elif Beta > 1.0 and Beta <= 1.999:
        etainf = 0.125*Beta**2 - 0.9*Beta + 5.436
        
    else:
        print "Error Beta out of value"
    
    return [s0, s1, etainf]

def dsfunction(phi0, phi1, s0, s1):
    """Method caluclating ds, based on the secant method
    
        Args:
            phi0(float): phi(s0)
            phi1(float): phi(s1)
            s0(floa): first value of s
            s1(floa): second value of s
    
        Returns:
            ds(float): s-s1
    """
    
    if (abs(phi1-phi0)>0.0):   
        return    -phi1 *(s1 - s0)/float(phi1 - phi0)
    else:
        return 0.0

def falknerSkan(f, eta):
    """The Falkner Skan equation reduced to a system of first order equations
    
    Args:
        f(array): an array containg f and its derivatives up to second order. (RHS)
        eta: independent variable
    Returns:
        dfdeta(array): df/d(eta)
    """
    
    dfdeta = np.zeros_like(f)
    dfdeta[0] = f[1]
    dfdeta[1] = f[2]
    dfdeta[2] = -(f[0]*f[2] + Beta*(1 - f[1]**2))
    
    return dfdeta

Betasep = -0.19883768 # value of Beta at which sepparation occurs
Beta = 1 # Beta angle randomly chosen to be one. Check with different values
eta0 = 0 # start of computational domain
N = 100 # number of elements
    
solverList = [euler, heun, rk4] 
solver = solverList[2] 
    
[s0, s1, etainf] = svalue(Beta) 
print 's0 = {}  s1 = {} and etainf = {} \n'.format(s0, s1, etainf)
eta = np.linspace(eta0, etainf, N + 1) # allocate space

Y0 = np.array([0, 0, s0]) # initial values
f1 = solver(falknerSkan, Y0, eta) 
phi0 = f1[-1, 1] - 1 # function that should equal zero for correct value of S

itmax, epsi, ds = 10, 10**-10, 1 # numerical tollerance limits in iteration process
it = 0

while abs(ds) > epsi and it < itmax:
    
    """Iteration loop: A new value of s is calculated with the secant method based on
        phi0, phi1, s0 and s1.
    """
    
    it = it + 1
    Y0 = np.array([0, 0, s1])
    f1 = solver(falknerSkan, Y0 , eta)
    phi1 = f1[-1,1]-1 
    ds = dsfunction(phi0, phi1, s0 , s1) 
    s = s1 + ds
    s0 = s1
    s1 = s
    phi0 = phi1
    print 'it = {}  s = {} and ds = {} \n'.format(it, s, ds)
    
legendList=[] 

plot(f1[:,1],eta)
legendList.append('dimensionless shear stress ')
plot(f1[:,2],eta)
legendList.append('dimensionless velocity')
#Add the labels
legend(legendList,loc='best',frameon=False) # Add the legends
ylabel('$\eta$')
xlabel("$f ',  f ''$")
grid(b=True, which='both', color='0.65',linestyle='-')

show()
