# src/src-ch5/startupv2.py


import numpy as np
from matplotlib.pyplot import *
from scipy.special import jv as besselj
from TRIdiagonalSolvers import tdma
from Startupfuncs import analyticSolution


# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)


def solut(theta, D, nmax, N, wold):
    """ Method that sets up the tridiagonal implicitt scheme to be solved with tdma using an iteration process
        Args:
            theta(float): parameter between 0 and 1 seperating different schemes (FTCS, Laasonen, Crank-Nicolson...)
            D(float): Numerical diffusion number
            nmax(int): number of timesteps
            N(int): number of equations
            wold(array): solution from previous iteration

    
        Returns:
            wnew(array): iterated solution at timestep nmax
    """
    
    tmp1 = D*(1. - theta)
    
    for n in range(nmax):
        a = np.zeros(N)
        b = np.zeros(N)
        c = np.zeros(N)
        d = np.zeros(N)

        for j in range(1,N):
            fac = 0.5 / (j + 1)
            a[j] = -D*theta*(1 - fac)
            b[j] = (1 + 2*D*theta)
            c[j] = -D*theta*(1 + fac)
            d[j] = tmp1*((1 - fac)*wold[j - 1] + (1 + fac)*wold[j + 1]) + (1 - 2*tmp1)*wold[j]
        b[0] = 1 + 4*D*theta/3
        c[0] = -4*D*theta/3
        d[0] = (1-4*tmp1/3.)*wold[0] + 4*tmp1*wold[1]/3.
        wnew = tdma(a, b, c, d)
        wold = np.append(wnew, [0])
    
    return wnew


"""script setting up the startup of a flow in a pipe (Szymanski's problem
    the equation of w is given by:
        dw/dt = w'' + w'/r, w = w(r,t), 0 < r < 1
    velocity field: u(r,t) = us(r) -w(r,t)
    where us = 1 - r^2
    Boundary conditions: w(0) = (4w(1)-w(2))/3)
    Initial condition: w(r,0) = 1- r^2 = us
    Timestep : dt = D*h^2, h = dr
    
    using nested functions:
        theta = 1: FTCS-scheme
        theta = 1/2: Crank-Nicolsen
        theta = 1: Laasonen
"""
Np1 = 50 # No. of parts
h = 1./Np1 # length of r-step
N = Np1 - 1 # No. of equations
D = 0.4 # numerical diffusion number
dt = D*h**2
nmax = 500 # No. of time-steps
tid = nmax*dt

print "*********************************"
print "*   theta = 0 FTSC-Scheme       *"
print "*   theta = 1 Laasonen          *"
print "*   theta = 1/2 Crank-Nicolsen  *"
print "*********************************"
print "\n No. of time-steps..........\n", nmax
print " r-step length..........\n", h
print "\n Diffusuion-number D..........\n", D
print "\n Timestep....................\n", dt
print "\n Elapsed time..........\n", tid

wnew = np.zeros(N)
wanalytic = np.zeros(Np1 + 1)
r = np.linspace(h, 1, Np1)

thetaFTSC, thetaLaasonen, thetaCrank = 0., 1., 0.5

wold = 1 - r**2 # initial Pouseille profile
wFTSC = solut(thetaFTSC, D, nmax, N, wold) 

wold = 1 - r**2
wLaasonen = solut(thetaLaasonen, D, nmax, N, wold) 

wold = 1 - r**2
wCrank = solut(thetaCrank, D, nmax, N, wold)

eps = 2.2204*10**-16
r = np.append([0], r)
for k in range(len(r)):
    x = r[k]
    wanalytic[k] = analyticSolution(x, tid, eps)
    

# === add BC's ===
w0 = (4*wFTSC[0]-wFTSC[1])/3.
wFTSC = np.append(np.append(w0, wFTSC), [0])
w0 = (4*wLaasonen[0]-wLaasonen[1])/3.
wLaasonen = np.append(np.append(w0, wLaasonen), [0])
w0 = (4*wCrank[0]-wCrank[1])/3.
wCrank = np.append(np.append(w0, wCrank), [0])
wFTSC = 1 - r**2 - wFTSC
wLaasonen = 1 - r**2 - wLaasonen
wCrank = 1 - r**2 - wCrank
wanalytic = 1 - r**2 - wanalytic

print 'r ,  u(ftcs) , u(ftcs) , u(laasonen) , u(crannk) , u(analytic) '
for k1 in range(len(r)):
    
    print '{},   {},  {},  {}, {}'.format(r[k1], wFTSC[k1], wLaasonen[k1], wCrank[k1], wanalytic[k1])

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
#savefig('startupv2.pdf')
show()