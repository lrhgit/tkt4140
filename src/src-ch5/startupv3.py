# src/src-ch5/startupv3.py


import numpy as np
from matplotlib.pyplot import *
from scipy.special import jv as besselj
from Startupfuncs import analyticSolution


# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)

def solut(theta, D, nmax, N, wold):
    """ Method that sets up the tridiagonal implicitt scheme to be solved with modified version of tdma,
         using an iteration process
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
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    wnew = np.zeros(N)

    for j in range(1, N):
        fac = 0.5 / (j)
        a[j] = -D*theta*(1 - fac)
        b[j] = (1 + 2*D*theta)
        c[j] = -D*theta*(1 + fac)
        
    b[0] = 1 + 4*D*theta
    c[0] = -4*D*theta
    # === elimination ===
    c[0] = c[0]/b[0]
    
    for j in range(1, N):
        b[j] = b[j] - a[j]*c[j - 1]
        c[j] = c[j]/b[j]
    
    for n in range(nmax):
        # === compute right hand side ===
        for j in range(1, N):
            fac = 0.5 / (j)
            d[j] = tmp1 *((1 - fac)*wold[j - 1] + (1 + fac)*wold[j+1]) + (1 - 2*tmp1)*wold[j]
        
        
        d[0] = wold[0] + 4*tmp1*(wold[1] - wold[0])

        # === backsubstitution ===
        d[0] = d[0]/b[0]
        for j in range(1, N):
            d[j] = (d[j] - a[j]*d[j - 1])/b[j]
        
        wnew[N-1] = d[N-1]
        
        for j in range(N - 2, -1, -1):
            wnew[j] = d[j] - c[j]*wnew[j + 1]
        
        wold = np.append(wnew,[0])
    
    return wnew



"""Main program setting up the startup of a flow in a pipe (Szymanski's problem)
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

print "\n *********************************"
print " *   theta = 0 FTSC-Scheme       *"
print " *   theta = 1 Laasonen          *"
print " *   theta = 1/2 Crank-Nicolsen  *"
print " *********************************"
print "\n No. of time-steps.......", nmax
print " r-step length.............", h
print " Diffusuion-number D.......", D
print " Timestep..................", dt
print " Elapsed time..............", tid
print "\n"

wnew = np.zeros(N)
wanalytic = np.zeros(N + 1)
r = np.linspace(0, 1, N +1)

thetaFTSC, thetaLaasonen, thetaCrank = 0, 1, 0.5

wold = 1 - r**2
wFTSC = solut(thetaFTSC, D, nmax, N, wold) 

wold = 1 - r**2
wLaasonen = solut(thetaLaasonen, D, nmax, N, wold) 

wold = 1 - r**2
wCrank = solut(thetaCrank, D, nmax, N, wold) 

eps = 2.2204*10**-16
for k in range(len(r)):
    x = r[k]
    wanalytic[k] = analyticSolution(x,tid,eps)

wFTSC = np.append(wFTSC, 0)
wLaasonen = np.append(wLaasonen, 0)
wCrank = np.append(wCrank, 0)
wFTSC = 1 - r**2 - wFTSC
wLaasonen = 1 - r**2 - wLaasonen
wCrank = 1 - r**2 - wCrank
wanalytic = 1 - r**2 - wanalytic

print ' r        u(FTSC)         u(Laasonen)         u(Crank)      u(analytic) '
for k1 in range(len(r)):
    
    print ' {},   {},  {},  {}, {}'.format(r[k1], wFTSC[k1], wLaasonen[k1], wCrank[k1], wanalytic[k1])



lstyle = ['r-', ':', '.', '-.']
legendList = ['wanalytic', 'wFTSC', 'wLaasonen', 'wCrank']
solutionList =[wanalytic, wFTSC, wLaasonen, wCrank]
i = 0
for solution in solutionList:
    plot(r, solution, lstyle[i])
    i = i +1

legend(legendList)
title('Flow in pipe')
xlabel('radius')
ylabel('velocity')
#savefig('startupv3.pdf')
show()

        

