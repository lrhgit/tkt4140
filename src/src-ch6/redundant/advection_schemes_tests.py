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

def init_box(x):
    """A smooth sin^2 function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = dx
    x_right = 0.2
    f = where((x>x_left) & (x<x_right), 1,f) 
    return f

def init_box_sine(x):
    """A smooth sin^2 function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = dx
    x_right = 0.2
    f = where((x>x_left) & (x<x_right), 1,f)
    x_left2 = 0.5
    x_right2 = 0.75
    f = where((x>x_left2) & (x<x_right2), np.sin(np.pi*(x-x_left2)/(x_right2-x_left2))**2,f) 
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
    
    def test_convergemce_MES():
        from numpy import log, log2, log10
        global c, dt, dx, a
        import time as timeModule
        
        a = 1.0 # wave speed
        tmin, tmax = 0.0, 1 # start and stop time of simulation
        xmin, xmax = 0.0, 2.0 # start and end of spatial domain
        Nx = 40 # number of spatial points
        c = 0.8 # courant number, need c<=1 for stability
    
        # Discretize
        x = np.linspace(xmin, xmax, Nx + 1) # discretization of space
        dx = float((xmax-xmin)/Nx) # spatial step size
        dt = c/a*dx # stable time step calculated from stability requirement
        time = np.arange(tmin, tmax + dt, dt) # discretization of time
        Ntds = 14 # number of grid/dt refinements
        # solve from tmin to tmax
        
        init_funcs = [init_step, init_sine4] # Select stair case function (0) or sin^4 function (1)
        f = init_funcs[1]
        #f = f3
        solvers = [ftbs, lax_friedrich, lax_wendroff, macCormack]
        spaceErrorList = [] #empty list of space Error
        spaceOrderList = []
        temporalErrorList = [] #empty list of temporal Error
        temporalOrderList = []
        maxErrorList = []
        globalErrorList = []
        hx = []
        ht = []
        
        for p in range(len(solvers)): 
            # create empty lists for all schemes for order and Error
            spaceErrorList.append([]), spaceOrderList.append([])
            temporalErrorList.append([]), temporalOrderList.append([])
            maxErrorList.append([]), globalErrorList.append([])
        #subplot ploting solutions for last timestep
        fig, ax = plt.subplots(Ntds, 1, sharex = True, squeeze=False)
        lstyle = ['b', 'r', 'g', 'm']
        lstyleAnalytic = 'c--'
        legendList = []

        print "\n"
        for n in range(Ntds):
            jump = 2**n
            hx.append(dx)
            ht.append(dt)
            
            cpustartT = timeModule.time()
            for k, solver in enumerate(solvers): # Solve for all solvers in list
                u = f(x)
                temporalError = 0
                temporalElements = 0
                maxerror = 0
                globalerror = 0
                globalElements = 0
                for i, t in enumerate(time[1:]):
                    
                    if k==0:
                        uanalytical = f(x-a*t) # compute analytical solution for this time step
                    
                    #u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
                    
                    u[1:-1] = solver(u[:]) # calculate numerical solution of interior
                    #u[-1] = u_bc(x[-1] - a*dt) # interpolate along a characteristic to find the boundary value
                    u[-1] = f(x[-1]-a*t)
                    u[0] = f(x[0]-a*t)
                    error = max(abs(u - f(x-a*t)))
                    globalerror += np.sum((u[1:-1] - f(x[1:-1]-a*t))**2)
                    globalElements += len(u[1:-1])
                    if error>maxerror:
                        maxerror = error
                    if float(i + 1)/jump == temporalElements + 1.0: #Only use samplepoints from every timestep of first iteration
                        #print t
                        temporalError += (u[Nx/2] - f(x[Nx/2]-a*t))**2
                        temporalElements += 1

                globalerror = np.sqrt(globalerror/globalElements)
                spaceError = np.sqrt(np.sum((u-uanalytical)**2)/(len(u)))
                spaceErrorList[k].append(spaceError)
                temporalError = np.sqrt(temporalError/temporalElements)
                temporalErrorList[k].append(temporalError)
                maxErrorList[k].append(maxerror)
                globalErrorList[k].append(globalerror)
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
            
            cpuendT = timeModule.time()
            
            print " finished iteration {0} of {1}, dx = {2}, dt = {3}, tsample = {4}, CPU-time = {5}".format(n+1, Ntds, dx, dt, t, cpuendT- cpustartT)
            # refine grid and dt:   
            Nx *= 2
            x = np.linspace(xmin, xmax, Nx+1) # discretization of space
            dx = float((xmax-xmin)/Nx) # spatial step size
            dt = c/a*dx # stable time step calculated from stability requirement
            time = np.arange(tmin, tmax + dt, dt) # discretization of time
        
        for k, solver in enumerate(solvers):
            with open(solver.func_name + "_errors.txt", 'w') as filename:
                filename.write('hx: '+ str(hx))
                filename.write('\n')
                filename.write('\n')
                filename.write('ht: '+ str(ht))
                filename.write('\n')
                filename.write('\n')
                filename.write('spaceError: '+ str(spaceErrorList[k]))
                filename.write('\n')
                filename.write('\n')
                filename.write('temporalError: '+ str(temporalErrorList[k]))
                filename.write('\n')
                filename.write('\n')
                filename.write('maxError: '+ str(maxErrorList[k]))
                filename.write('\n')
                filename.write('\n')
                filename.write('globalError: '+ str(globalErrorList[k]))
            filename.close()
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
        
    def test_CFL():
        
        # Change default values on plots
        LNWDT=2; FNT=10
        plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT
        
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
                                
                            if n ==len(solvers)-1:
                                ax[n][k].set_xlabel('x')
                                
                            if k==0:
                                ax[n][k].text(-0.5, 0.5, 'CFL = {0}'.format(c))
                                ax[n][k].set_ylabel('y')
                                
                    drawing = 0
        print "\n"
                            



        
    
    test_convergemce_MES()
    #test_CFL()
    plt.show()
