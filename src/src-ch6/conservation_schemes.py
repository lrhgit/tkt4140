# ../Kap6/advection_schemes.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy import interpolate
from numpy import where
from math import sqrt, pi



LNWDT=2; FNT=10
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT



# ftbs for conservation equation
def ftbs(u, t):
    """method that solves u(n+1), for the scalar conservation equation with source term:
        du/dt + dF/dx = RHS,
        where F = 0.5u^2 for the burger equation
        with use of the forward in time backward in space (upwind) scheme
        
        Args:
            u(array): an array containg the previous solution of u, u(n). (RHS)
            t(float): an array 
        Returns:
            u[1:-1](array): the solution of the interior nodes for the next timestep, u(n+1).
    """
    u[1:-1] = u[1:-1] -  (dt/dx)*(F(u[1:-1])-F(u[:-2])) + dt*RHS(t-0.5*dt, x[1:-1])
    return u[1:-1]

# lax_friedrich_Flux for conservation equation
def lax_friedrich_Flux(u, t):
    """method that solves u(n+1), for the scalar conservation equation with source term:
        du/dt + dF/dx = RHS,
        where F = 0.5u^2 for the burger equation
        with use of the lax-friedrich scheme
        
        Args:
            u(array): an array containg the previous solution of u, u(n). (RHS)
            t(float): an array 
        Returns:
            u[1:-1](array): the solution of the interior nodes for the next timestep, u(n+1).
    """
    u[1:-1] = (u[:-2] +u[2:])/2.0 -  dt*(F(u[2:])-F(u[:-2]))/(2.0*dx) + dt*(RHS(t, x[:-2]) + RHS(t, x[2:]))/2.0
    return u[1:-1]


# Lax_W_Two_Step for conservation equation
def Lax_W_Two_Step(u, t):
    """method that solves u(n+1), for the scalar conservation equation with source term:
        du/dt + dF/dx = RHS,
        where F = 0.5u^2 for the burger equation
        with use of the Two-step Lax-Wendroff scheme
        
        Args:
            u(array): an array containg the previous solution of u, u(n).
            t(float): time at t(n+1) 
        Returns:
            u[1:-1](array): the solution of the interior nodes for the next timestep, u(n+1).
    """
    ujm = u[:-2].copy() #u(j-1)
    uj = u[1:-1].copy() #u(j)
    ujp = u[2:].copy() #u(j+1)
    up_m = 0.5*(ujm + uj) - 0.5*(dt/dx)*(F(uj)-F(ujm)) + 0.5*dt*RHS(t-0.5*dt, x[1:-1] - 0.5*dx) #u(n+0.5dt,j-0.5dx)
    up_p = 0.5*(uj + ujp) - 0.5*(dt/dx)*(F(ujp)-F(uj)) + 0.5*dt*RHS(t-0.5*dt, x[1:-1] + 0.5*dx)#u(n+0.5dt,j+0.5dx)
    
    u[1:-1] = uj -(dt/dx)*(F(up_p) - F(up_m)) + dt*RHS(t-0.5*dt, x[1:-1])
    return u[1:-1]

# macCormack for conservation equation
def macCormack(u, t):
    """method that solves u(n+1), for the scalar conservation equation with source term:
        du/dt + dF/dx = RHS,
        where F = 0.5u^2 for the burger equation
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

# Lax_W_One_Step_J for conservation equation
def Lax_W_One_Step_J(u, t):
    """method that solves u(n+1), for the scalar conservation equation with source term:
        du/dt + dF/dx = RHS,
        where F = 0.5u^2 for the burger equation
        with use of a general Lax-Wendroff one step method
        
        Args:
            u(array): an array containg the previous solution of u, u(n). (RHS)
            t(float): an array 
        Returns:
            u[1:-1](array): the solution of the interior nodes for the next timestep, u(n+1).
    """
    ujm = u[:-2].copy() #u(j-1)
    uj = u[1:-1].copy() #u(j)
    ujp = u[2:].copy() #u(j+1)
    
    vjm = (dt/dx)*wavespeed(ujm, uj)
    vjp = (dt/dx)*wavespeed(uj, ujp)#A()
    
    Fjm = F(ujm)
    Fj = F(uj)
    Fjp = F(ujp)
    
    Qjm = RHS(t-dt, x[:-2])
    Qj = RHS(t-dt, x[1:-1])
    Qjp = RHS(t-dt, x[2:])
    Qjpnp = RHS(t, x[1:-1])
    
    Source = 0.5*dt*(Qj+Qjpnp)-0.25*dt*(vjp*(Qjp+Qj)-vjm*(Qj+Qjm))
    
    u[1:-1] = uj -0.5*(dt/dx)*(Fjp - Fjm) + 0.5*(dt/dx)*(vjp*(Fjp-Fj)-vjm*(Fj-Fjm)) + Source
    return u[1:-1]

# Lax_W_One_Step for conservation equation
def Lax_W_One_Step(u, t):
    """method that solves u(n+1), for the scalar conservation equation with source term:
        du/dt + dF/dx = RHS,
        where F = 0.5u^2 for the burger equation
        with use of a general Lax-Wendroff one-step method
        
        Args:
            u(array): an array containg the previous solution of u, u(n). (RHS)
            t(float): an array 
        Returns:
            u[1:-1](array): the solution of the interior nodes for the next timestep, u(n+1).
    """
    ujm = u[:-2].copy() #u(j-1)
    uj = u[1:-1].copy() #u(j)
    ujp = u[2:].copy() #u(j+1)
    ajp = wavespeed(uj, ujp)
    ajm = wavespeed(ujm, uj)
    
    vjp = ajp*dt/dx
    vjm = ajm*dt/dx

    Qjm = RHS(t-dt, x[:-2])
    Qj = RHS(t-dt, x[1:-1])
    Qjp = RHS(t-dt, x[2:])
    Qjpnp = RHS(t, x[1:-1])
    
    Source = 0.5*dt*(Qj+Qjpnp)-0.25*dt*(vjp*(Qjp+Qj)-vjm*(Qj+Qjm))
    
    u[1:-1] = uj -(dt/dx)*(F_LW(uj, ujp, ajp) - F_LW(ujm, uj, ajm)) + Source
    return u[1:-1]




def F(u):
    """method that calculates the Flux term for the burger equation F = 0.5u^2, for
    
        Args:
            u(array): an array containg the solution of u
        Returns:
            Flux(array): The flux 0.5u^2.
    """
    
    Flux = (u**2)/2.0
    
    return Flux

def F_LW(uj, ujp, a):
    """method that calculates the general Flux term F_LW for the Lax-Wendroff Flux
        F(j+0.5,n+0.5) = F(uj) + 0.5*a*(1-a*dt/dx)*(u{j+1}-u{j})
        where an aproximation of the numerical wavespeed at j+0.5
    
        Args:
            u(array): an array containg the solution of u
        Returns:
            Flux(array): The Lax-Wendroff Flux .
    """
    
    Flux = F(uj) + 0.5*a*(1-a*dt/dx)*(ujp-uj)
    
    return Flux


def wavespeed(uj, ujp, tol=1e-10):
    """method that approximates the numerical wavespeed at j+0.5:
                    { F(u{j+1})-F(u{j})/(u{j+1}-u{j}) if u{j+1} != u{1}
        a{j+0.5} =  {
                    {F'(u{j}                          if u{j+1} == u{1}
        
        the method assures conservation
        Args:
            uj(array): an array containg the solution of u{j}
            ujp(array): an array containg the solution of u{j+1}
        Returns:
            A(array): the numerical wavespeed at j+0.5
    """
    A = np.ones_like(uj)
    for n in range(len(uj)):
        if abs(ujp[n]-uj[n])<tol:
            A[n]=uj[n]
        else:
            A[n]=(F(ujp[n])-F(uj[n]))/(ujp[n]-uj[n])

    return A

def wavespeed_non_conservative(uj, ujp):
    """method that approximates the numerical wavespeed at j+0.5:
                    
        a{j+0.5} = 0.5*(u{j+1} + u{j})
        
        the method does not assure conservation
        Args:
            uj(array): an array containg the solution of u{j}
            ujp(array): an array containg the solution of u{j+1}
        Returns:
            A(array): the numerical wavespeed at j+0.5
    """
    
    A = 0.5*(uj+ujp)

    return A



if __name__ == '__main__':
    

    def example():
        
        def init_sine2(x):
            """A smooth sin^2 function between x_left and x_right"""
            f = np.zeros_like(x)
            x_left = 0.25
            x_right = 0.75
            xm = (x_right-x_left)/2.0
            f = where((x>x_left) & (x<x_right), np.sin(np.pi*(x-x_left)/(x_right-x_left))**2,f)
            f = f + 2 
            return f
        
        def init_sine3(x):
            """A smooth sin^2 function between x_left and x_right"""
            f = np.zeros_like(x)
            x_left = 0.25
            x_right = 0.75
            xm = (x_right-x_left)/2.0
            f = where((x>x_left) & (x<x_right), np.sin(2*np.pi*(x-x_left)/(x_right-x_left)),f)
            f = f + 1 
            return f
        
        def init_box(x):
            """A smooth sin^2 function between x_left and x_right"""
            f = np.zeros_like(x)
            x_left = dx
            x_right = 0.2
            f = where((x>x_left) & (x<x_right), 1,f) 
            return f
        
        def RHS(x,t):
            return 0
        
        
        # Discretize
        global dx, dt, RHS, x
        a = 2.0 # wave speed
        tmin, tmax = 0.0, 1.0 # start and stop time of simulation
        xmin, xmax = 0.0, 2.0 # start and end of spatial domain
        Nx = 401 # number of spatial points
        c = 0.9 # courant number, need c<=1 for stability
        
        
        x = np.linspace(xmin, xmax, Nx+1) # discretization of space
        dx = float((xmax-xmin)/Nx) # spatial step size
        dt = c/a*dx # stable time step calculated from stability requirement
        Nt = int((tmax-tmin)/dt) # number of time steps
        time = np.arange(tmin, tmax, dt) # discretization of time
        
        f = init_sine3
        # solve from tmin to tmax
        
        solvers = [macCormack, lax_friedrich_Flux, Lax_W_One_Step, Lax_W_Two_Step, ftbs]

        u_solutions=np.zeros((len(solvers),len(time),len(x)))
        
        
        
            
        for k, solver in enumerate(solvers): # Solve for all solvers in list
            u = f(x)
            un = np.zeros((len(time), len(x))) # holds the numerical solution
        
            for i, t in enumerate(time[1:]):
                
                    
                u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
                
                u[1:-1] = solver(u[:],t) # calculate numerical solution of interior
                u[-1] = u_bc(x[-1] - a*dt) # interpolate along a characteristic to find the boundary value
                
                un[i,:] = u[:] # storing the solution for plotting
            
            u_solutions[k,:,:] = un
        
        
        
        ### Animation 
         
        # First set up the figure, the axis, and the plot element we want to animate
        fig = plt.figure()
        ax = plt.axes(xlim=(xmin,xmax), ylim=(np.min(un), np.amax(un)*1.2))
        
        lines=[]     # list for plot lines for solvers and analytical solutions
        legends=[]   # list for legends for solvers and analytical solutions
        
        for solver in solvers:
            line, = ax.plot([], [])
            lines.append(line)
            legends.append(solver.func_name)
        
        
        plt.xlabel('x-coordinate [-]')
        plt.ylabel('Amplitude [-]')
        plt.legend(legends, loc='best', frameon=False)
         
        # initialization function: plot the background of each frame
        def init():
            for line in lines:
                line.set_data([], [])
            return lines,
        
        # animation function.  This is called sequentially

        
        def animate_alt(i):
            for k, line in enumerate(lines):

                line.set_data(x, u_solutions[k,i,:])
            return lines,
        
#        Writer = animation.writers['ffmpeg']
#        writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate_alt, init_func=init, frames=Nt, interval=100, blit=False)
#        anim.save('../figs/burger.mov', writer=writer)
         
        plt.show()
    
    def convergence_test():
        from sympy import symbols, diff, latex, lambdify, sin
        from numpy import log
        global  x, dx, dt, RHS
        
        tmin, tmax = 0.0, 1 # start and stop time of simulation
        xmin, xmax = 0.0, 1 # start and end of spatial domain
        Nx = 40 # number of spatial points
        CFL = 0.8#0.8#0.015625 # courant number, need c<=1 for stability
        
        # symbolic computations:
        a = 1
        meanU = 5 
        t, x = symbols('t x')
         
        Analytic = meanU + sin(2*pi*x-a*t)**4 # expression in term of variables t and x
        F = 0.5*Analytic**2
        dAnalyticdt = diff(Analytic, t)
        dFdx = diff(F, x)
        RHS = (dAnalyticdt + dFdx)

        
        Analytic = lambdify([t, x], Analytic, np)
        RHS = lambdify([t, x], RHS, np)
            

        UmeanPlusC = max([abs(a + meanU), abs(a + meanU)])
        
        # Discretize
        x = np.linspace(xmin, xmax, Nx + 1) # discretization of space
        dx = float((xmax-xmin)/Nx) # spatial step size
        #meanU_plus_A = c + meanU # calculate wavespeed + Meanvelocity
        dt = dx*CFL/UmeanPlusC # stable time step calculated from stability requirement
        #tmax = dt
        time = np.arange(tmin, tmax + dt, dt) # discretization of time

        solvers = [ftbs, lax_friedrich_Flux, Lax_W_Two_Step, Lax_W_One_Step, macCormack]   
        Ntds = 3

        ErrorList = [] #empty list of temporal Error
        OrderList = []
        
        for p in range(len(solvers)): 
            # create empty lists for all schemes for order and spaceError
            ErrorList.append([]), OrderList.append([])

        for n in range(Ntds):

            for k, solver in enumerate(solvers): # Solve for all solvers in list
                u = Analytic(0, x)
                error = 0
                elements = 0

                for i, t in enumerate(time[1:]):
                    
                    #u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
                    
                    u[1:-1] = solver(u[:], t) # calculate numerical solution of interior
                    #u[-1] = u_bc(x[-1] - A*dt) # interpolate along A characteristic to find the boundary value
                    u[0] = Analytic(t, x[0])
                    u[-1] = Analytic(t, x[-1])

                    
                    error += np.sum((u-Analytic(t, x))**2)
                    elements += len(u)  

                error = np.sqrt(error/elements)
                ErrorList[k].append(error)

                if n>0:
                    #calculate order approximation
                    OrderList[k].append(log(ErrorList[k][n-1]/ErrorList[k][n])/log(2))
                            
            
            print " finished iteration {0} of {1}, dx = {2}, dt = {3}, tsample = {4}".format(n+1, Ntds, dx, dt, t)
            Nx *= 2
            x = np.linspace(xmin, xmax, Nx + 1) # discretization of space
            dx = float((xmax-xmin)/Nx) # spatial step size
            dt = dx*CFL/UmeanPlusC # stable time step calculated from stability requirement
            time = np.arange(tmin, tmax + dt, dt) # discretization of time
            #print "\n"
        
        
        # Plot calculated rms spaceErrors and order approximations
        
        fig , axarr = plt.subplots(2, 1, squeeze=False)
        lstyle = ['b', 'r', 'g', 'c', 'm--']
        legendList = []
        
        N = Nx/2**(Ntds + 1)
        N_list = [N*2**i for i in range(1, Ntds+1)]
        N_list = np.asarray(N_list)
        
        epsilonN = [i for i in range(1, Ntds)]
        print epsilonN
        epsilonlist = ['$\epsilon_{0} , \epsilon_{1}$'.format(str(i), str(i+1)) for i in range(1, Ntds)]
        
        for k, solver in enumerate(solvers):
            axarr[0][0].plot(N_list, np.log10(np.asarray(ErrorList[k])),lstyle[k])
            axarr[1][0].plot(epsilonN, OrderList[k],lstyle[k])
            
            legendList.append(solver.func_name)
        
        
        
        axarr[1][0].axhline(1.0, xmin=0, xmax=epsilonN[-1], linestyle=':', color='k')
        axarr[1][0].axhline(2.0, xmin=0, xmax=epsilonN[-1], linestyle=':', color='k')
        
        #legends and labels
        #plt.suptitle('test_MES_convergence(): Results from convergence test using MES')
        axarr[1][0].set_ylim(0, 3)
        axarr[1][0].set_xlim(epsilonN[0], epsilonN[-1])
        axarr[0][0].set_ylabel(r'$\log_{10}(E)$')
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
        
    example()
    #convergence_test()
    plt.show()
