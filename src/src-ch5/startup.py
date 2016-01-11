# src/src-ch5/startup.py;TRIdiagonalSolvers.py @ git@lrhgit/tkt4140/src/src-ch5/TRIdiagonalSolvers.py;


import numpy as np
from matplotlib.pyplot import *
from scipy.special import jv as besselj
from TRIdiagonalSolvers import tdma
from Startupfuncs import analyticSolution


# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)



def startup():
    """Method setting up the startup of a flow in a pipe (Szymanski's problem)
        the equation of w is given by:
            dw/dt = w'' + w/r, w = w(r,t), 0 < r < 1
        velocity field: u(r,t) = us(r) -w(r,t)
        where us = 1 - r^2
        BCs : w(r,0) = 1- r^2 = us
        Timestep : dt = D*h^2, h = dr
        
        using nested functions:
            theta = 1: FTCS-scheme
            theta = 1/2: Crank-Nicolsen
            theta = 1: Laasonen
    """
    N = 50 # No. of parts
    h = 1./N # length of t-step
    D = 0.4 # numerical diffusion number
    dt = D*h**2
    nmax = 500 # No. of time-steps
    tid = nmax*dt
    
    print "*********************************"
    print "*   theta = 0 FTSC-Scheme       *"
    print "*   theta = 1 Laasonen          *"
    print "*   theta = 1/2 Crank-Nicolsen  *"
    print "*********************************"
    print "\n No. of time-steps..........", nmax
    print " r-step length..........", h
    print " Diffusuion-number D..........", D
    print " Timestep....................", dt
    print " Elapsed time..........", tid
    
    wnew = np.zeros(N)
    wa = np.zeros(N+1)
    r = np.linspace(0, 1, N +1)
    
    # === theta = 0 FTCS scheme ===
    theta = 0
    wold = 1 - r**2
    w1 = solut(theta, D, nmax, N, wold)
    
    # === theta = 1  FTCS Laasonen ===
    theta = 1
    wold = 1 - r**2
    w2 = solut(theta, D, nmax, N, wold) # Laason
    
    # === theta =  0.5 Crank-Nicolsen ===
    theta = 0.5
    wold = 1 - r**2
    w3 = solut(theta, D, nmax, N, wold) # Crank Nicholson
    
    eps = 2.2204*10**-16
    for k in range(len(r)):
        x = r[k]
        wa[k] = analyticSolution(x, tid, eps)
    
    w1 = np.append(w1, 0)
    w2 = np.append(w2, 0)
    w3 = np.append(w3, 0)
    w1 = 1 - r**2 - w1
    w2 = 1 - r**2 - w2
    w3 = 1 - r**2 - w3
    wa = 1 - r**2 - wa
    print 'r ,  u(ftcs) , u(ftcs) , u(laasonen) , u(crannk) , u(analytic) '
    for k1 in range(len(r)):
        
        print '{},   {},  {},  {}, {}'.format(r[k1], w1[k1], w2[k1], w3[k1], wa[k1])
    
    

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

def solut(theta, D, nmax, N, wold):
    """ Method that sets up the tridiagonal implicitt scheme to be solved with tdma using an iteration process
        Args:
            theta(float): parameter between 0 and 1 seperating different schemes (FTCS, Laasonen, Crank-Nicolson...)
            D(float): Numerical diffusion number
            n(int): number of iterations
            N(int): number of equations
            wold(array): solution from previous iteration

    
        Returns:
            wnew(array): iterated solution 
    """
    
    tmp1 = D*(1. - theta)
    
    for n in range(nmax):
        a = np.zeros(N)
        b = np.zeros(N)
        c = np.zeros(N)
        d = np.zeros(N)

        for j in range(1, N):
            fac = 0.5 / (j)
            a[j] = -D*theta*(1 - fac)
            b[j] = (1 + 2*D*theta)
            c[j] = -D*theta*(1 + fac)
            d[j] = tmp1*((1 - fac)*wold[j - 1] + (1 + fac)*wold[j + 1]) + (1 - 2*tmp1)*wold[j]
            
        b[0] = 1 + 4*D*theta
        c[0] = -4*D*theta
        d[0] = wold[0] + 4*tmp1*(wold[1] - wold[0])
        wnew = tdma(a, b, c, d)
        wold = np.append(wnew, [0])
    
    return wnew
            

startup()

