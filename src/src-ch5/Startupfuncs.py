# chapter5/src-ch6/startup.py;TRIdiagonalSolvers.py @ git@lrhgit/tkt4140/allfiles/digital_compendium/chapter5/src-ch6/TRIdiagonalSolvers.py;


import numpy as np
from scipy.special import jv as besselj


def analyticSolution(r,t,eps):
    
    """ Method that calculates the analytic solution to the differential equation:
            dw/dt = w'' + w/r, w = w(r,t), 0 < r < 1
            where t is the time.
            The accuracy can be controlled by changing the value of sumtol1, which depicts the
            relative accuracy.
            
        Args:
            r(float): radial coordinat
            t(float): time
            eps(float): numerical tollerance

    
        Returns:
            w(float): velocity, us - ur

    """
    
    if r == 1:
        w = 0
    elif r < eps and t < eps:
        w = 1
    else:
        # calculate first element seperately
        n = 1
        lam1 = j0zero(n)
        arg1 = r*lam1
        term1 = np.exp(-t*lam1**2)*besselj(0,arg1)/(besselj(1,lam1)*lam1**3)
        sum1 = term1
        sumtol = 1*10**-8
        test = 1
        
        while test > sumtol:
            n = n + 1
            lamn = j0zero(n)
            arg = r*lamn
            term = np.exp(-t*lamn**2)*besselj(0,arg)/(besselj(1,lamn)*lamn**3)
            sum1 = sum1 + term
            test = np.abs(term/term1)
        
        w = 8 * sum1
    
    return w
        

def j0zero(s):
    
    """ method that computes root number s, s = 1,2,...,
        of the Besselfunction J0 where zj0 is the root,
        using table-values, asymptotic formulae and Newton-Raphsons method.
        Rellative error app. 1.0e-14 to 1.0e-15
        Args:
            s(integer): solutian array from last iteration

        Returns:
            zj0(float): root No. s of besselfunction
    """
    
    sz = np.array([2.4048225557695773, 5.520078110286311, 8.653727912911012, 11.79153443901428, 14.93091770848779])
    b0 = (s - 0.25)*np.pi
    b08 = 0.125/b0
    b082 = b08**2
    z0 = b0 + b08*(1. - b082*(124./3 -  b082*120928./15))
    
    if s <= 5:
        zj0 = sz[s - 1]
    elif s > 30:
        zj0 = z0
    else:
        dz0 = besselj(0, z0)/besselj(1,  z0)
        z0 = z0 + dz0
        zj0 = z0
    return zj0