# ../Kap6/advection_schemes.py

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from numpy import where

# function defining the initial condition

def init_step(x):
    """Assigning a value of 1.0 for values less than 0.1"""
    f = np.zeros_like(x)
    f[np.where(x <= 0.1)] = 1.0
    return f

def init_sine4(x):
    """A smooth sin^ function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = 0.4
    x_right = 0.6
    xm = (x_right-x_left)/2.0
    f = where((x>x_left) & (x<x_right), np.sin(np.pi*(x-x_left)/(x_right-x_left))**4,f) 
    return f

# Forward in time, backward in space
def ftbs(u): # forward time backward space
    u[1:-1] = (1-c)*u[1:-1] + c*u[:-2]
    return u[1:-1]

# Lax-Wendroff
def lax_wendroff(u): 
    u[1:-1] = c/2.0*(1+c)*u[:-2] + (1-c**2)*u[1:-1] - c/2.0*(1-c)*u[2:]
    return u[1:-1]

# Lax-Friedrich Flux formulation
def lax_friedrich_Flux(u):
    u[1:-1] = (u[:-2] +u[2:])/2.0 -  dt*(F(u[2:])-F(u[:-2]))/(2.0*dx)
    return u[1:-1] 

# Lax-Friedrich Advection
def lax_friedrich(u):
    u[1:-1] = (u[:-2] +u[2:])/2.0 -  c*(u[2:] - u[:-2])/2.0
    return u[1:-1] 

# macCormack for advection quation
def macCormack(u):
    up = u.copy()
    up[:-1] = u[:-1] - c*(u[1:]-u[:-1])
    u[1:] = .5*(u[1:]+up[1:] -  c*(up[1:]-up[:-1]))
    return u[1:-1] 

def F(u):
    return a*u


if __name__=='__main__':
    # Tests:
    # Convergence test
    def test_convergence_MES():
        from numpy import log 
        from math import sqrt
        global c, dt, dx, a
        
        # Change default values on plots
        LNWDT=2; FNT=13
        plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT
        
        a = 1.0 # wave speed
        tmin, tmax = 0.0, 1 # start and stop time of simulation
        xmin, xmax = 0.0, 2.0 # start and end of spatial domain
        Nx = 160 # number of spatial points
        c = 0.8 # courant number, need c<=1 for stability
    
        # Discretize
        x = np.linspace(xmin, xmax, Nx + 1) # discretization of space
        dx = float((xmax-xmin)/Nx) # spatial step size
        dt = c/a*dx # stable time step calculated from stability requirement
        time = np.arange(tmin, tmax + dt, dt) # discretization of time
        
        init_funcs = [init_step, init_sine4] # Select stair case function (0) or sin^4 function (1)
        f = init_funcs[1]
        
        solvers = [ftbs, lax_friedrich, lax_wendroff]
        
        errorDict = {} # empty dictionary to be filled in with lists of errors
        orderDict = {}
        
        for solver in solvers:
            errorDict[solver.func_name] = [] # for each  solver(key) assign it with value=[], an empty list (e.g: 'ftbs'=[])
            orderDict[solver.func_name] = []
            
        hx = [] # empty list of spatial step-length
        ht = [] # empty list of temporal step-length
        
        Ntds = 5 # number of grid/dt refinements
        # iterate Ntds times:
        for n in range(Ntds):
            hx.append(dx)
            ht.append(dt)
            
            for k, solver in enumerate(solvers): # Solve for all solvers in list
                u = f(x) # initial value of u is init_func(x)
                error = 0 
                samplePoints = 0
                
                for i, t in enumerate(time[1:]):
                    u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
                    u[1:-1] = solver(u[:]) # calculate numerical solution of interior
                    u[-1] = u_bc(x[-1] - a*dt) # interpolate along a characteristic to find the boundary value

                    error += np.sum((u - f(x-a*t))**2) # add error from this timestep
                    samplePoints += len(u)

                error = sqrt(error/samplePoints) # calculate rms-error
                errorDict[solver.func_name].append(error)
                
                if n>0:
                    previousError = errorDict[solver.func_name][n-1]
                    orderDict[solver.func_name].append(log(previousError/error)/log(2))            
            
            print " finished iteration {0} of {1}, dx = {2}, dt = {3}, tsample = {4}".format(n+1, Ntds, dx, dt, t)
            # refine grid and dt:   
            Nx *= 2
            x = np.linspace(xmin, xmax, Nx+1) # new x-array, twice as big as the previous
            dx = float((xmax-xmin)/Nx) # new spatial step size, half the value of the previous
            dt = c/a*dx # new stable time step 
            time = np.arange(tmin, tmax + dt, dt) # discretization of time
        
        # Plot error-values and corresponding order-estimates:
        fig , axarr = plt.subplots(2, 1, squeeze=False)
        lstyle = ['b', 'r', 'g', 'm']
        legendList = []
        
        N = Nx/2**(Ntds + 1)
        N_list = [N*2**i for i in range(1, Ntds+1)]
        N_list = np.asarray(N_list)
        
        epsilonN = [i for i in range(1, Ntds)]
        epsilonlist = ['$\epsilon_{0} , \epsilon_{1}$'.format(str(i), str(i+1)) for i in range(1, Ntds)]
        
        for k, solver in enumerate(solvers):
            axarr[0][0].plot(N_list, np.log10(np.asarray(errorDict[solver.func_name])),lstyle[k])
            axarr[1][0].plot(epsilonN, orderDict[solver.func_name],lstyle[k])
            
            legendList.append(solver.func_name)
        
        
        
        axarr[1][0].axhline(1.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
        axarr[1][0].axhline(2.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
        
        #legends and labels
        #plt.suptitle('test_MES_convergence(): Results from convergence test using MES')
        axarr[1][0].set_ylim(0, 3)
        axarr[0][0].set_ylabel(r'$\log_{10}\left(\sqrt{\frac{1}{N}\sum\left(\hat{f}-f\right)^2}\right)$')
        axarr[0][0].set_xlabel(r'$N_x$')
        axarr[0][0].set_xticks(N_list)
        axarr[0][0].set_yticks([-1, -2, -3, -4, -5])
        axarr[1][0].set_xticks(epsilonN)
        axarr[1][0].set_xticklabels(epsilonlist)
        axarr[1][0].set_yticks([1,2])
        axarr[1][0].set_xlabel(r'$\left(\epsilon_{n-1} , \epsilon_{n}\right)$')
        print epsilonN, epsilonlist

        axarr[1][0].set_ylabel(r'$log_{2}( \frac{{\epsilon}_{n-1}}{{\epsilon}_n} )$', fontsize=25)
        axarr[0][0].legend(legendList, frameon=False, loc='best')
        fig.tight_layout()
#        plt.savefig('../figs/Advection_schemes_convergence.png', transparent=True, frameon=False)
    
    # CFL test
    def test_CFL():
        
        # Change default values on plots
        LNWDT=2; FNT=10
        plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT
        plt.rcParams.update({'figure.autolayout': True})
        
        global c
        a = 1.0 # wave speed
        tmin, tmax = 0.0, 0.5 # start and stop time of simulation
        xmin, xmax = 0.0, 1 # start and end of spatial domain
        Nx = 80 # number of spatial points
        
        x = np.linspace(xmin, xmax, Nx+1) # discretization of space
        dx = x[1] - x[0] # spatial step size
        # solve from tmin to tmax
    
        solvers = [ftbs,lax_friedrich, macCormack]
        init_funcs = [init_step, init_sine4]
        CFLnumbers = np.linspace(0.25, 1, 4)
        
        lstyle = ['b', 'r', 'g', 'c']
        lstyleAnalytic = 'k--'

        for f in init_funcs:
            fig, ax = plt.subplots(len(CFLnumbers), len(solvers), sharey = True, sharex = True, squeeze=False)
            
            for n, c in enumerate(CFLnumbers):
                # adjust dt for new cfl number
                dt = c/a*dx # stable time step calculated from stability requirement
                Nt = int((tmax-tmin)/dt) # number of time steps
                time = np.linspace(tmin, tmax, Nt) # discretization of time
                drawing = 0
                totalDrawings = 4
                jump = int(0.5/totalDrawings/dt) # only draw certain amount of drawings for each solver
                
                for k, solver in enumerate(solvers): # Solve for all solvers in list
                    u = f(x)
    
                    for i, t in enumerate(time[1:]):
                        
                        
                        uanalytical = f(x-a*t) # compute analytical solution for this time step
                        
                        u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
                        
                        u[1:-1] = solver(u[:]) # calculate numerical solution of interior
                        u[-1] = u_bc(x[-1] - a*dt) # interpolate along a characteristic to find the boundary value
                        #u[-1] = f(x[-1]-a*t)
                        #u[0] = 1
                        if (round(i/jump)-drawing)==0 and drawing<totalDrawings:
                            drawing += 1
                            ax[n][k].plot(x, u, lstyle[k])
                            ax[n][k].plot(x, uanalytical, lstyleAnalytic)
                            
                            if n==0:
                                ax[n][k].set_title(solver.func_name)
                                
                            if n ==len(CFLnumbers)-1:
                                ax[n][k].set_xlabel('x')
                                
                            if k==0:
                                if f.func_name=='init_step':
                                    ax[n][k].text(-0.55, 0.5, 'CFL = {0}'.format(c))
                                else:
                                    ax[n][k].text(-0.55, 0.5, 'CFL = {0}'.format(c))
                                ax[n][k].set_ylabel('\n \n \n y')
                                
                    drawing = 0
            fig.tight_layout()
#            plt.savefig('../figs/Advection_Schemes_CFL_' + f.func_name + '.png', transparent=True, frameon=False)
        print "\n"
                            



        
    
    test_convergence_MES()
    test_CFL()
    plt.show()
