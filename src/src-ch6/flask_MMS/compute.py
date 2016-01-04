import os, glob, math
import time as timemodule
import numpy as np

from math import sqrt, pi
import matplotlib.pyplot as plt
from matplotlib import animation

LNWDT=4; FNT=20
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT


def compute_MMS_PDE(Expression='exp(cos(x - t))*(cos(t))/exp(1)', c=1.0, meanU=4,
                    tdomain='[0,2*pi]', xdomain='[0,2*pi]', Nx=160., CFL=0.9 ):


    tmin, tmax = eval(tdomain)[0], eval(tdomain)[1] # start and stop time of simulation
    xmin, xmax = eval(xdomain)[0], eval(xdomain)[1] # start and end of spatial domain

    
    Analytic, RHS, titlestring = Analytic_RHS(Expression, meanU, xmin, xmax, tmin, tmax) #initilize manufactured solution and RHS
    minValue, maxValue = findExtremalvaluesNumeric(Analytic, xmin, xmax, tmin, tmax)
    
    UmaxPlusC = max([abs(c + meanU), abs(c + meanU)]) # find critical velocity for calculating dt from CFL (CFL = max(U+c)*dt/dx
    # Discretize
    x = np.linspace(xmin, xmax, Nx + 1) # discretization of space
    dx = float((xmax-xmin)/Nx) # spatial step size
    #meanU_plus_A = c + meanU # calculate wavespeed + Meanvelocity
    dt = dx*CFL/UmaxPlusC # stable time step calculated from stability requirement
    Nt = int((tmax-tmin)/dt) # number of time steps
    time = np.linspace(tmin, tmax, Nt) # discretization of time
    ActualCFL = UmaxPlusC*dt/dx
        
    # solve from tmin to tmax
    solvers = [macCormack]
    #initialize solution matrices
    u_solutions=np.zeros((len(solvers), len(time), len(x)))
    uanalytical = np.zeros((len(time), len(x))) # holds the analytical solution
    un = np.zeros((len(time), len(x)))
        
    for k, solver in enumerate(solvers): # Solve for all solvers in list
        u = Analytic(0, x)
        un = np.zeros((len(time), len(x))) # holds the numerical solution
        un[0, :] = u
        uanalytical[0, :] = Analytic(0, x)
        for i, t in enumerate(time[1:]):
            #print 'i, t: ', i, t
             
            if k==0:
                uanalytical[i+1,:] = Analytic(t, x) # compute analytical solution for this time step
                 
            #u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interpolate at right boundary
             
            u[1:-1] = solver(u[:], t, x, dt, dx, RHS) # calculate numerical solution of interior
            u[-1] = Analytic(t, x[-1]) # interpolate along a characteristic to find the boundary value
            u[0] = Analytic(t, x[0])
            un[i+1,:] = u[:] # storing the solution for plotting
         
        u_solutions[k,:,:] = un
    
    
    animationfile = createAnimation(uanalytical, u_solutions, dt, Nt,
                                   x, xmin, xmax,  meanU, solvers,
                                   titlestring)
    
    return animationfile, ActualCFL, dt, dx

def Analytic_RHS(Expression, meanU, xmin, xmax, tmin, tmax):
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
    Analytic = meanU + eval(Expression) # expression in term of variables t and x
    F = 0.5*Analytic**2
    dAnalyticdt = diff(Analytic, t)
    dFdx = diff(F, x)
    RHS = (dAnalyticdt + dFdx)
    #minValue, maxValue = findExtremalvaluesAnalytic(Analytic, xmin, xmax, tmin, tmax)
    #print minValue, maxValue
    titlestring = latex(Analytic, mode='inline')

    Analytic = lambdify([t, x], Analytic, np)
    RHS = lambdify([t, x], RHS, np)
    return Analytic, RHS, titlestring

def findExtremalvaluesAnalytic(f, xmin, xmax, tmin, tmax): #f,xmin,xmax,tmin,tmax
    """Method that find extreme values of a function f(t, x) on the domain limited by xmin, xmax, tmin, tmax.
        In the interior both dfdt(t,x) and dfdx(t,x) has to be zero for f to have an extreme minimum or maximum at (t,x)
        In addition the maximums could be on the boundaries.
        Args:
            f(funciton): a sympy function
            xmin(float): start of xdomain
            xmax(float): end of xdomain
            tmin(float): start of timedomain
            tmax(float): end of timedomain
        Returns:
            minvalue(float): the maximum value of f in the domain
            maxvalue(float): the maximum value of f in the domain
    """
    from sympy import symbols, diff, solve, lambdify, pi, sqrt, exp, sin, cos, tan, Abs
    # differentiate and solve to locate extreme values of f:
    t, x = symbols('t x')
    dfdt = diff(f, t)
    dfdx = diff(f,x)
    print "finished differentiating"
    print "dfdt", dfdt
    roots1 = solve(dfdt, t, x)
    roots2 = solve(dfdx, t, x)
    roots1dict = {"time":[], "space":[]}
    roots2dict = {"time":[], "space":[]}
    print "found roots: ", roots1
    print "found roots: ", roots2
    # Save roots in dictionaries 
    for root in roots1:
        try:
            troot = float(root[t])
            if  tmin < troot < tmax:
                roots1dict["time"].append(troot)
        except:
            xroot = float(root[x])
            if  xmin < troot < xmax:
                roots1dict["space"].append(xroot)
    
    for root in roots2:
        try:
            troot = float(root[t])
            if  tmin < troot < tmax:
                roots2dict["time"].append(troot)#, root.value()
        except:
            xroot = float(root[x])
            if  xmin<xroot<xmax:
                roots2dict["space"].append(xroot)
    
            
    f = lambdify([t, x], f, np)
    #extract maximum values from interior:
    fextremevalues = np.array([])
    for time in roots1dict["time"]:
        tempvalue = f(time, np.asarray(roots2dict["space"]))
        fextremevalues = np.append(fextremevalues, tempvalue)

    for time in roots2dict["time"]:
        tempvalue = f(time, np.asarray(roots1dict["space"]))
        fextremevalues = np.append(fextremevalues, tempvalue)   
   

    maxvalue = np.max(fextremevalues)
    minvalue = np.min(fextremevalues)
    # check values on the border numerically:
    time = np.linspace(tmin,tmax, 100)
    space = np.linspace(xmin,xmax,100)
    
    bcmax = max([np.max(f(xmin,time)), np.max(f(xmax,time)), np.max(f(space,tmin)), np.max(f(space,tmax))])
    bcmin = min([np.min(f(xmin,time)), np.min(f(xmax,time)), np.min(f(space,tmin)), np.min(f(space,tmax))])
    
    maxvalue = max([maxvalue, bcmax])
    minvalue = min([minvalue, bcmin])
    
    return minvalue, maxvalue

def findExtremalvaluesNumeric(f, xmin, xmax, tmin, tmax): #f,xmin,xmax,tmin,tmax
    """Method that find extreme values of a function f(t, x) on the domain limited by xmin, xmax, tmin, tmax.
        by evaluating the values of the domain numerically on a coarse grid.
        Args:
            f(funciton): a python function
            xmin(float): start of xdomain
            xmax(float): end of xdomain
            tmin(float): start of timedomain
            tmax(float): end of timedomain
        Returns:
            minvalue(float): the maximum value of f in the domain
            maxvalue(float): the maximum value of f in the domain
    """
    # discretize with a course grid
    maxvalue = -float("inf")
    minvalue = float("inf")
    time = np.linspace(tmin,tmax, 100)
    space = np.linspace(xmin,xmax,100)
    for t in time:
        tempmaxvalue = np.max(f(t,space))
        tempminvalue = np.min(f(t,space))
        
        if tempmaxvalue > maxvalue:
            maxvalue = tempmaxvalue
        if tempminvalue < minvalue:
            minvalue = tempminvalue
    
    return minvalue, maxvalue

def F(u):
    """method that calculates the Flux term F = 0.5u^2
    
        Args:
            u(array): an array containg the solution of u
        Returns:
            Flux(array): The flux 0.5u^2.
    """
    Flux = u.copy()
    Flux[:] = (u**2)/2
    
    return Flux

def macCormack(u, t, x, dt, dx, RHS):
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
    up[:-1] = u[:-1] - (dt/dx)*(F(u[1:]) - F(u[:-1])) + dt*RHS(t, x[:-1])
    u[1:] = .5*(u[1:] + up[1:] -  (dt/dx)*(F(up[1:]) - F(up[:-1])) + dt*RHS(t, x[1:]))
    return u[1:-1]

    ### Animation 
    
    # First set up the figure, the axis, and the plot element we want to animate
    
def createAnimation(uanalytical, u_solutions, dt, Nt,
                   x, xmin, xmax,  meanU, solvers,
                   titlestring):
    """Method that creates and saves the animation of uNumeric and uAnalytic. The .mp4 file is stored in a directory "static". 
    """
    fig = plt.figure()
    ax = plt.axes(xlim=(xmin,xmax), ylim=((np.min(uanalytical)-meanU)*1.1 + meanU, (np.max(uanalytical)-meanU)*1.1 + meanU))
    
    lines=[]     # list for plot lines for solvers and analytical solutions
    legends=[]   # list for legends for solvers and analytical solutions
    
    line, = ax.plot([], [], 'b') #add extra plot line for analytical solution
    lines.append(line)
    legends.append('Analytical')
    
    for solver in solvers:
        line, = ax.plot([], [],'m--')
        lines.append(line)
        legends.append(solver.func_name)
    
    plt.title(titlestring)
    plt.xlabel('x-coordinate [-]')
    plt.ylabel('Amplitude [-]')
    plt.legend(legends, loc=3, frameon=False)
     
    # initialization function: plot the background of each frame
    def init():
        #line.set_data([], [])
        for line in lines:
            line.set_data([], [])
        return lines,
    
    # animation function.  This is called sequentially
    
    def animate_alt(i):
        i = i*jump
#         for k, line in enumerate(lines):
#             if (k==len(lines)-1):
#                 line.set_data(x, uanalytical[i,:])
#             else:
#                 line.set_data(x, u_solutions[k,i,:])
        lines[0].set_data(x, uanalytical[i,:])
        lines[1].set_data(x, u_solutions[0,i,:])

        return lines,
    #FFwriter = animation.FFMpegWriter()
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    
    jump = int(0.06/dt) + 1 # reasonable number of timesteps visualized.
    if jump<1:
        jump = 1
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate_alt, init_func=init, frames=Nt/jump, interval=100, blit=False)
    if not os.path.isdir('static'):
        os.mkdir('static')
    else:
        #remove old animations
        for filename in glob.glob(os.path.join('static','*mp4')):
            os.remove(filename)
#     print time.time()
    print timemodule.time()
    animationfile = os.path.join('static', str(timemodule.time()) + '.mp4')
    print animationfile
    anim.save(animationfile, writer=writer)
    
    return animationfile

def convergence_test(Expression, meanU, xmin, xmax, tmin, tmax,  UmaxPlusC, Nxstart=40, CFLstart=9, test='both'):
    from matplotlib.pyplot import plot, show, legend, hold,rcParams,rc, figure, axhline, close,\
    xticks, title, xlabel, ylabel, savefig, axis, grid
    from numpy import log2, log10
#     global a, tmin, tmax, xmin, xmax, Nx, dx, c, x, dx, dt, Nt, time, Analytic, RHS
#     a = 5.0 # wave speed
#     tmin, tmax = 0.0, np.pi # start and stop time of simulation
#     xmin, xmax = 0.0, 2*np.pi # start and end of spatial domain
#     Nx = 40*8 # number of spatial points
#     c = 1#0.015625 # courant number, need c<=1 for stability
    Analytic, RHS, titlestring = Analytic_RHS(Expression, meanU, xmin, xmax, tmin, tmax)
    x = np.linspace(xmin, xmax, Nx+1) # discretization of space
    dx = float((xmax-xmin)/Nx) # spatial step size
    dt = c/a*dx # stable time step calculated from stability requirement
    Nt = int((tmax-tmin)/dt) # number of time steps
    time = np.linspace(tmin, tmax, Nt) # discretization of time
    
    x = np.linspace(xmin, xmax, Nxstart + 1) # discretization of space
    dx = float((xmax-xmin)/Nx) # spatial step size
    #meanU_plus_A = c + meanU # calculate wavespeed + Meanvelocity
    dt = dx*CFLstart/UmaxPlusC # stable time step calculated from stability requirement
    Nt = int((tmax-tmin)/dt) # number of time steps
    time = np.linspace(tmin, tmax, Nt) # discretization of time
    ActualCFL = UmaxPlusC*dt/dx



    # solve from tmin to tmax

    solvers = [macCormack]
    solver_order={}
    Ndts = 5 # Number of times to refine delta x in convergence test
    for k, solver in enumerate(solvers): # Solve for all solvers in list
        solver_error = [] 
        for n in range(Ndts+1):
            u = Analytic(0, x)
            error = 0
            Measurementpoints = 0
            print "len(x): ", len(x)
            print "\n"
            print "------ it = ", n,  "------"
            print "------ dt = ", dt, "------"
            print "------ dx = ", dx, "------"
            print "------ CFL = ", a*dt/dx, "------"
            print "\n"
            
            for i, t in enumerate(time[1:]):
                #print 'i, t: ', i, t
                uanalytical = Analytic(t, x) # compute analytical solution for this time step   
                #u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
                u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
                u[1:-1] = solver(u[:], t) # calculate numerical solution of interior
                
                u[-1] = u_bc(x[-1] - a*dt) # interpolate along a characteristic to find the boundary value
                #u[-1] = Analytic(t, x[-1]) # interpolate along a characteristic to find the boundary value
                u[0] = Analytic(t, x[0])
                temperror = (u-uanalytical)**2
                
                error = error + np.sum(temperror)
                Measurementpoints = Measurementpoints + len(uanalytical)
            figure()
            plot(x, uanalytical,'b')
            plot(x, u,'r--')
            print "max u: ", np.amax(u)
            error = sqrt(error/(Measurementpoints))
            solver_error.append(error)
            test = 'dt'
            if test=='both':
                Nx *=2
                x = np.linspace(xmin, xmax, Nx+1)
                u = Analytic(0, x)
                dx = float((xmax-xmin)/Nx) # spatial step size
                dt = dx*c/a # stable time step calculated from stability requirement
                Nt = int((tmax-tmin)/dt) # number of time steps
                time = np.linspace(tmin, tmax, Nt) # discretization of time
            elif test == 'dx':
                Nx *=2
                x = np.linspace(xmin, xmax, Nx+1)
                u = Analytic(0, x)
                dx = float((xmax-xmin)/Nx) # spatial step size
                #dt = dx*c/a # stable time step calculated from stability requirement
                #Nt = int((tmax-tmin)/dt) # number of time steps
                #time = np.linspace(tmin, tmax, Nt) # discretization of time
            elif test=='dt':
#                 Nx *=2
#                 x = np.linspace(xmin, xmax, Nx+1)
#                 u = Analytic(0, x)
#                 dx = float((xmax-xmin)/Nx) # spatial step size
                dt = 0.5*dt # stable time step calculated from stability requirement
                Nt = int((tmax-tmin)/dt) # number of time steps
                time = np.linspace(tmin, tmax, Nt) # discretization of time
        solver_order[solver.func_name] = solver_error
        solver_error = np.asarray(solver_error)
        print solver_error
        order_approx = log10(solver_error[:-1]/solver_error[1:])/log10(2)
        print order_approx
        figure()
        plot(solver_error)
        title('global RMS error')
        xlabel('iteration')
        figure()
        plot(order_approx)
        title('order approximation')
        xlabel('iteration')
        show()
def compute_vib(A, b, w, T, resolution=500):
    """Return filename of plot of the damped_vibration function."""
    t = linspace(0, T, resolution+1)
    y = damped_vibrations(t, A, b, w)
    plt.figure()  # needed to avoid adding curves in plot
    plt.plot(t, y)
    plt.title('A=%g, b=%g, w=%g' % (A, b, w))

    # Make Matplotlib write to BytesIO file object and grab
    # return the object's string
    from io import BytesIO
    figfile = BytesIO()
    plt.savefig(figfile, format='png')
    figfile.seek(0)  # rewind to beginning of file
    import base64
    figdata_png = base64.b64encode(figfile.getvalue())
    figfile = BytesIO()
    plt.savefig(figfile, format='svg')
    figfile.seek(0)
    figdata_svg = '<svg' + figfile.getvalue().split('<svg')[1]
    figdata_svg = unicode(figdata_svg,'utf-8')
    return figdata_png, figdata_svg

def gamma_density(x, a, h, A):
    # http://en.wikipedia.org/wiki/Gamma_distribution
    xA = x/float(A)
    return abs(h)/(math.gamma(a)*A)*(xA)**(a*h-1)*exp(-xA**h)

def gamma_cumulative(x, a, h, A):
    # Integrate gamma_density using the Trapezoidal rule.
    # Assume x is array.
    g = gamma_density(x, a, h, A)
    r = zeros_like(x)
    for i in range(len(r)-1):
        r[i+1] = r[i] + 0.5*(g[i] + g[i+1])*(x[i+1] - x[i])
    return r

def compute_gamma(a=0.5, h=2.0, A=math.sqrt(2), resolution=500):
    """Return plot and mean/st.dev. value of the gamma density."""
    gah = math.gamma(a + 1./h)
    mean = A*gah/math.gamma(a)
    stdev = A/math.gamma(a)*math.sqrt(
        math.gamma(a + 2./h)*math.gamma(a) - gah**2)
    x = linspace(0, 7*stdev, resolution+1)
    y = gamma_density(x, a, h, A)
    plt.figure()  # needed to avoid adding curves in plot
    plt.plot(x, y)
    plt.title('a=%g, h=%g, A=%g' % (a, h, A))
    # Make Matplotlib write to BytesIO file object and grab
    # return the object's string
    from io import BytesIO
    figfile = BytesIO()
    plt.savefig(figfile, format='png')
    figfile.seek(0)  # rewind to beginning of file
    import base64
    figdata_density_png = base64.b64encode(figfile.getvalue())
    figfile = BytesIO()
    plt.savefig(figfile, format='svg')
    figfile.seek(0)
    figdata_density_svg = '<svg' + figfile.getvalue().split('<svg')[1]
    figdata_density_svg = unicode(figdata_density_svg,'utf-8')

    y = gamma_cumulative(x, a, h, A)
    plt.figure()
    plt.plot(x, y)
    plt.grid(True)
    figfile = BytesIO()
    plt.savefig(figfile, format='png')
    figfile.seek(0)
    figdata_cumulative_png = base64.b64encode(figfile.getvalue())
    figfile = BytesIO()
    plt.savefig(figfile, format='svg')
    figfile.seek(0)
    figdata_cumulative_svg = '<svg' + figfile.getvalue().split('<svg')[1]
    figdata_cumulative_svg = unicode(figdata_cumulative_svg,'utf-8')
    return figdata_density_png, figdata_cumulative_png, \
           figdata_density_svg, figdata_cumulative_svg, \
           '%.2f' % mean, '%.2f' % stdev
