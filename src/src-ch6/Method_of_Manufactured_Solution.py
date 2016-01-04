# ../Kap6/advection_schemes.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy import interpolate
from numpy import where
from math import sqrt, pi



LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT



# macCormack for burgers quation
def macCormack(u, t):
    """method that solves u(n+1) based on the modified Burgers equaion (RHS != 0):
        du/dt + dF/dx = RHS,
        where F = 0.5u^2
        with use of the MacCormack scheme
        
        Args:
            u(array): an array containg the previous solution of u, u(n). (RHS)
            t(float): an array 
        Returns:
            u[1:-1](array): the solution of the interior nodes for the next timestep, u(n+1).
    """
    up = u.copy()
    up[:-1] = u[:-1] - (dt/dx)*(F(u[1:]) - F(u[:-1])) + dt*RHS(t-0.5*dt, x[:-1])
    u[1:] = .5*(u[1:] + up[1:] -  (dt/dx)*(F(up[1:]) - F(up[:-1])) + dt*RHS(t-0.5*dt, x[1:]))
    return u[1:-1]

def lax_friedrich_Flux(u, t):
    u[1:-1] = (u[:-2] +u[2:])/2.0 -  dt*(F(u[2:])-F(u[:-2]))/(2.0*dx) + dt*(RHS(t, x[:-2]) + RHS(t, x[2:]))/2.0
    return u[1:-1]

def macCormack2(u, t):
    up = u.copy()
    up[:-1] = u[:-1] - c*(u[1:]-u[:-1]) + dt*RHS(t - 0.5*dt, x[:-1])
    u[1:] = .5*(u[1:]+up[1:] -  c*(up[1:]-up[:-1])) + 0.5*dt*RHS(t-0.5*dt, x[1:])
    return u[1:-1]

def F(u):
    """method that calculates the Flux term F = 0.5u^2
    
        Args:
            u(array): an array containg the solution of u
        Returns:
            Flux(array): The flux 0.5u^2.
    """
    Flux = u.copy()
    Flux[:] = (u**2)/2.0
    
    return Flux

def F2(u):
    """method that calculates the Flux term F = 0.5u^2
    
        Args:
            u(array): an array containg the solution of u
        Returns:
            Flux(array): The flux 0.5u^2.
    """
    Flux = u.copy()
    Flux[:] = a*u
    
    return Flux

def Analytic_RHS(Expression):
    """method that initialize the analytical solution (coefficient function) that is wanted to be used in the MMS.
        With use of the sympy module the analytical solution of u is created, and RHS is calculated based on the 
        modified Burgers equation:
        du/dt + dF/dx = RHS,
        where F = 0.5u^2
        Args:
            Expression(str): a function Analytic(t,x) readable by sympy, that is chosen to be the analytical solution
            meanU(float): mean velocity
            xmin(float): start of xdomain
            xmax(float): end of xdomain
            tmin(float): start of timedomain
            tmax(float): end of timedomain
        Returns:
            Analytic(method): a python function Analytic(t, x) calculated from uAnalytic = meanU + uPerturbation.
            RHS(method): a python function RHS(t, x) 
    """        
    from sympy import symbols, diff, integrate, Rational, lambdify
    from sympy import latex, pi, sqrt, exp, sin, cos, tan, Abs
    print "hey"
    t, x = symbols('t x')
    #Expression = str(meanU) + 
    Analytic = eval(Expression) # expression in term of variables t and x
    F = 0.5*Analytic**2
    dAnalyticdt = diff(Analytic, t)
    dFdx = diff(F, x)
    RHS = (dAnalyticdt + dFdx)
    print "Source term aka. RHS: ", RHS
    #minValue, maxValue = findExtremalvaluesAnalytic(Analytic, xmin, xmax, tmin, tmax)
    #print minValue, maxValue
    titlestring = latex(Analytic, mode='inline')
    
    Analytic = lambdify([t, x], Analytic, np)
    RHS = lambdify([t, x], RHS, np)
    return Analytic, RHS, titlestring

def Analytic_RHS2(Expression):
    """method that initialize the analytical solution (coefficient function) that is wanted to be used in the MMS.
        With use of the sympy module the analytical solution of u is created, and RHS is calculated based on the 
        modified Burgers equation:
        du/dt + dF/dx = RHS,
        where F = 0.5u^2
        Args:
            Expression(str): a function Analytic(t,x) readable by sympy, that is chosen to be the analytical solution
            meanU(float): mean velocity
            xmin(float): start of xdomain
            xmax(float): end of xdomain
            tmin(float): start of timedomain
            tmax(float): end of timedomain
        Returns:
            Analytic(method): a python function Analytic(t, x) calculated from uAnalytic = meanU + uPerturbation.
            RHS(method): a python function RHS(t, x) 
    """        
    from sympy import symbols, diff, integrate, Rational, lambdify
    from sympy import latex, pi, sqrt, exp, sin, cos, tan, Abs

    t, x = symbols('t x')
    #Expression = str(meanU) + 
    Analytic = eval(Expression) # expression in term of variables t and x
    F = a*Analytic
    dAnalyticdt = diff(Analytic, t)
    dFdx = diff(F, x)
    RHS = (dAnalyticdt + dFdx)
    print "RHS", RHS
    #minValue, maxValue = findExtremalvaluesAnalytic(Analytic, xmin, xmax, tmin, tmax)
    #print minValue, maxValue
    titlestring = latex(Analytic, mode='inline')

    Analytic = lambdify([t, x], Analytic, np)
    RHS = lambdify([t, x], RHS, np)
    return Analytic, RHS, titlestring

def RMS_error(fnumeric, fanalytic):
    N = len(fnumeric)
    error = np.abs((fnumeric - fanalytic)/np.max(fanalytic))
    error = np.sum(error)/N
    
    return error

def abs_error(fnumeric, fanalytic):
    
    error = np.abs((fnumeric - fanalytic))
    error = np.max(error)
    
    return error

if __name__ == '__main__':
    
    def convergence_test():
        from matplotlib.pyplot import plot, show, legend, hold,rcParams,rc, figure, axhline, close,\
        xticks, title, xlabel, ylabel, savefig, axis, grid
        from numpy import log2, log10, log
        global a, tmin, tmax, xmin, xmax, Nx, dx, c, x, dx, dt, Nt, time, Analytic, RHS
        a = 1.0 # wave speed
        tmin, tmax = 0.0, 1 # start and stop time of simulation
        xmin, xmax = 0.0, 1 # start and end of spatial domain
        Nx = 10 # number of spatial points
        CFL = 0.8#0.8#0.015625 # courant number, need c<=1 for stability
        c = CFL
        
        
        Expression = '5 + sin(2*pi*x-t)**2'
        meanU = 5 
        Analytic, RHS, titlestring = Analytic_RHS(Expression)
        
        #minValue, maxValue = findExtremalvaluesNumeric(Analytic, xmin, xmax, tmin, tmax)
        #print Analytic
        UmeanPlusC = max([abs(a + meanU), abs(a + meanU)])
        
                # Discretize
        x = np.linspace(xmin, xmax, Nx + 1) # discretization of space
        dx = float((xmax-xmin)/Nx) # spatial step size
        #meanU_plus_A = c + meanU # calculate wavespeed + Meanvelocity
        dt = dx*CFL/UmeanPlusC # stable time step calculated from stability requirement
        #tmax = dt
        time = np.arange(tmin, tmax + dt, dt) # discretization of time
        ActualCFL = UmeanPlusC*dt/dx
    
        # solve from tmin to tmax
    
        solvers = [macCormack, lax_friedrich_Flux]   
        Ntds = 5
        spaceErrorList = [] #empty list of space Error
        spaceOrderList = []
        temporalErrorList = [] #empty list of temporal Error
        temporalOrderList = []
        
        for p in range(len(solvers)): 
            # create empty lists for all schemes for order and spaceError
            spaceErrorList.append([]), spaceOrderList.append([])
            temporalErrorList.append([]), temporalOrderList.append([])
        
        #subplot ploting solutions for last timestep
        fig, ax = plt.subplots(Ntds, 1, sharex = True, squeeze=False)
        lstyle = ['b', 'r--', 'g', 'm']
        lstyleAnalytic = 'c--'
        legendList = []
        LNWDT=2; FNT=8
        plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT
        print "\n"
        for n in range(Ntds):
            jump = 2**n

            for k, solver in enumerate(solvers): # Solve for all solvers in list
                u = Analytic(0, x)
                temporalError = 0
                temporalElements = 0
                temporalList = [0]
                temporaltimes = [0]
                temporalAnalyticlist = [0]
                for i, t in enumerate(time[1:]):
#                     if i ==0:
#                         print "i, t", i, t
#                         print "time[jump+1::jump]", time[jump::jump]
                    
                    #u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
                    
                    u[1:-1] = solver(u[:], t) # calculate numerical solution of interior
                    #u[-1] = u_bc(x[-1] - a*dt) # interpolate along a characteristic to find the boundary value
                    u[0] = Analytic(t, x[0])
                    u[-1] = Analytic(t, x[-1])
                    if float(i + 1)/jump ==temporalElements + 1.0: #Only use samplepoints from every timestep of first iteration
                        #print t
                        temporalError += (u[Nx/2] - Analytic(t, x[Nx/2]))**2
                        temporalElements += 1
                        
                    temporalList.append(u[Nx/2])
                    temporalAnalyticlist.append(Analytic(t, x[Nx/2]))
                    temporaltimes.append(t)
                uanalytical = Analytic(t, x)
                #spaceError = np.max(np.abs(u[0::jump]-uanalytical[0::jump]))
                spaceError = np.sqrt(np.sum((u-uanalytical)**2)/(len(u)))
                
                if k==0:
                    temporalError = np.sqrt(temporalError/temporalElements)
                    #print "temporalElements: ", temporalElements
                    spaceError = np.sqrt(np.sum((u-uanalytical)**2)/(len(u)))
                if k==1:
                    temporalError = np.sqrt(np.sum((np.asanyarray(temporalList)[jump::jump]-np.asanyarray(temporalAnalyticlist)[jump::jump])**2)/(len(temporalList[jump::jump])))
                    #print "len(temporalList[jump::jump]): ", len(temporalList[jump::jump])
                    spaceError = np.sqrt(np.sum((u[0::jump]-uanalytical[0::jump])**2)/(len(u[0::jump])))
                spaceErrorList[k].append(spaceError)
                temporalErrorList[k].append(temporalError)
                #print "temporaltimes[jump::jump]", temporaltimes[jump::jump]
                print " finished iteration {0} of {1}, dx = {2}, dt = {3}, tsample = {4}".format(n+1, Ntds, dx, dt, t)
                if n>0:
                    #calculate order approximation
                    spaceOrderList[k].append(log(spaceErrorList[k][n-1]/spaceErrorList[k][n])/log(2))
                    temporalOrderList[k].append(log(temporalErrorList[k][n-1]/temporalErrorList[k][n])/log(2))
                    
                ax[n][0].plot(x, u, lstyle[k])
                
                if n==Ntds-1:
                    legendList.append(solver.func_name)
                    
            if n==Ntds-1:
                legendList.append('Analytic')
            ax[n][0].plot(x,uanalytical,lstyleAnalytic)
            ax[n][0].set_title('Iteration {0}'.format(n + 1))
            #ax[n][0].set_xlim(1.2,1.8)
            ax[0][0].legend(legendList, frameon=False)            
            
            #print " finished iteration {0} of {1}, dx = {2}, dt = {3}, tsample = {4}".format(n+1, Ntds, dx, dt, t)
            #print " finished iteration {0} of {1}, dxa = {2}, dta = {3}, tsamplea = {4}".format(n+1, Ntds, x[1]-x[0], time[1]-time[0], t)
            # refine grid and dt:   
            Nx *= 2
            x = np.linspace(xmin, xmax, Nx + 1) # discretization of space
            dx = float((xmax-xmin)/Nx) # spatial step size
            dt = dx*CFL/UmeanPlusC # stable time step calculated from stability requirement
            time = np.arange(tmin, tmax + dt, dt) # discretization of time
            #print "\n"
        
        
        # Plot calculated rms spaceErrors and order approximations
        
        legendList = []
        fig , axarr = plt.subplots(2, 2, squeeze=False)
        for k, solver in enumerate(solvers):
            
            axarr[0][0].plot(spaceErrorList[k],lstyle[k])
            axarr[0][1].plot(temporalErrorList[k],lstyle[k])
            axarr[1][0].plot(spaceOrderList[k],lstyle[k])
            axarr[1][1].plot(temporalOrderList[k],lstyle[k])
            legendList.append(solver.func_name)
        plt.suptitle('test_MES_convergence(): Results from convergence test using Method of Exact Solution')

        axarr[1][0].axhline(1.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
        axarr[1][0].axhline(2.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
        axarr[1][1].axhline(1.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
        axarr[1][1].axhline(2.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
        axarr[1][0].set_ylim(0, 3)
        axarr[1][1].set_ylim(0, 3)
        axarr[0][0].set_ylabel('rms Error')
        axarr[0][0].set_title('space Error')
        axarr[1][0].set_ylabel('rms Error')
        axarr[0][1].set_title('temporal Error')
        axarr[1][0].set_ylabel('order')
        axarr[1][0].set_title('space order')
        axarr[1][1].set_ylabel('order')
        axarr[1][1].set_title('temporal order')
        axarr[0][1].legend(legendList, frameon=False)
        print "\n"
    
    convergence_test()
    plt.show()
