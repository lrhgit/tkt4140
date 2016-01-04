# chapter2/src-ch1/ODEschemes.py

import numpy as np
from matplotlib.pyplot import plot, show, legend, hold,rcParams,rc, figure, axhline, close,\
    xticks, xlabel, ylabel, savefig, axis, grid
    
from ODEschemes import euler, heun, rk4
from sympy import symbols, diff, integrate, Rational, lambdify
from math import sqrt, pi

# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)

def Manufactured_solution():
    
    from sympy import exp, sin, cos
    t, A, b = symbols('t A b')
    
    f = A*exp(b*t)
    dfdt = diff(f, t)
    d2fdt = diff(dfdt, t)
    print 'f, dfdt, d2fdt: ', f, dfdt, d2fdt
    
    RHS = f+d2fdt
    print 'RHS', RHS
    f = lambdify([t, A, b], f)
    dfdt = lambdify([t, A, b], dfdt)
    RHS = lambdify([t, A, b], RHS)
    
    T = 1
    Nsteps = 50
    time = np.linspace(0, T, Nsteps + 1)
    
    
    def ODEset(y,t):
        yout = np.zeros_like(y)
        #print 't, RHS(t, 2, 3)',t, RHS(t, 2, 0.5)
        
        #print 't, RHS(t, 2, 3)-y[0]: ', RHS(t, 2, 2)-f(t, 2, 0.5)
        yout[:] = [y[1], RHS(t, 2, 0.5)-y[0]]
        
        return yout
        
    solver = euler
    
    y0 = np.array([f(0, 2, 0.5), dfdt(0, 2, 0.5)])
    print y0
    
    u = solver(ODEset, y0, time)
    fnumeric = u[:,0]
    fanalytic = 2*np.exp(0.5*time)
    legends = []
    plot(time, fnumeric)
    legends.append('fnumeric')
    plot(time, fanalytic)
    legends.append('fanalytic')
    legend(legends)
    show()
    
def Manufactured_solution2():
    
    from sympy import exp, sin, cos
    t, sigma, mu = symbols('t sigma mu')
    
    f = (1/(sigma*sqrt(2*pi)))*exp(-((t-mu)**2)/(2*sigma**2))
    dfdt = diff(f, t)
    d2fdt = diff(dfdt, t)
    print 'f, dfdt, d2fdt: ', f, dfdt, d2fdt
    
    RHS = f+d2fdt
    print 'RHS', RHS
    f = lambdify([t, sigma, mu], f)
    dfdt = lambdify([t, sigma, mu], dfdt)
    RHS = lambdify([t, sigma, mu], RHS)
    
    t0 = -1.5
    Nsteps = 50
    time = np.linspace(-1.5, 2.5, Nsteps + 1)
    
    
    def ODEset(y,t):
        yout = np.zeros_like(y)
        #print 't, RHS(t, 2, 3)',t, RHS(t, 2, 0.5)
        
        #print 't, RHS(t, 2, 3)-y[0]: ', RHS(t, 2, 2)-f(t, 2, 0.5)
        yout[:] = [y[1], RHS(t, 2, 0.5)-y[0]]
        
        return yout
        
    solver = heun
    
    y0 = np.array([f(t0, 2, 0.5), dfdt(t0, 2, 0.5)])
    print y0
    
    u = solver(ODEset, y0, time)
    fnumeric = u[:,0]
    fanalytic = (1/(2*sqrt(2*pi)))*np.exp(-((time-0.5)**2)/(2*2**2))
    legends = []
    plot(time, fnumeric)
    legends.append('fnumeric')
    plot(time, fanalytic)
    legends.append('fanalytic')
    legend(legends)
    show()
    
def Manufactured_solution3(Expression='(1/(sigma*sqrt(2*pi)))*exp(-((t-mu)**2)/(2*sigma**2))', 
                           variable='t', sigma=2, mu=0.5, Domain=[-1.5,2.5], solver='heun', Nsteps = 50 
                           ):
    
    
    #print Expression
    from sympy import exp, sin, cos
    commaSeperatedVariables = ''

 
    #print command
    print "solving equation f''' + f''*f + f' = RHS"
    print "which lead to f''' = RHS - f''*f - f"
    t = symbols(variable)
    
    f = eval(Expression)
    dfdt = diff(f, t)
    d2fdt = diff(dfdt, t)
    d3fdt = diff(d2fdt, t)
    
    
    RHS = d3fdt + dfdt*d2fdt + f

    f = lambdify([t], f)
    dfdt = lambdify([t], dfdt)
    d2fdt = lambdify([t], d2fdt)
    RHS = lambdify([t], RHS)
    
    t0 = Domain[0]
    tend = Domain[1]
    
    time = np.linspace(t0, tend, Nsteps + 1)
    
    
    def func(y,t):
        yout = np.zeros_like(y)
        yout[:] = [y[1], y[2], RHS(t) -y[0]- y[1]*y[2]]
        
        return yout
        
    
    
    y0 = np.array([f(t0), dfdt(t0), d2fdt(t0)])
    
    solver = eval(solver)
    u = solver(func, y0, time)
    fnumeric = u[:,0]
    fanalytic = np.zeros_like(fnumeric)
    i = 0
    for t in time:
        fanalytic[i] = f(t)
        i = i + 1
        
    legends = []
    plot(time, fnumeric)
    legends.append('fnumeric')
    plot(time, fanalytic)
    legends.append('fanalytic')
    legend(legends)
    show()



Manufactured_solution3()

# a = 0.2
# b = 3.0
# u_exact = lambda t: a*t   +  b


    

def f_local(u,t):
    """A function which returns an np.array but less easy to read
    than f(z,t) below. """
    return np.asarray([a + (u - u_exact(t))**5])

def f(z, t):
    """Simple to read function implementation """
    return [a + (z - u_exact(t))**5]


def test_ODEschemes():
    """Use knowledge of an exact numerical solution for testing."""
    from numpy import linspace, size

    tol = 1E-15
    T = 2.0  # end of simulation
    N = 20  # no of time steps
    time = linspace(0, T, N+1)


    z0 = np.zeros(1)
    z0[0] = u_exact(0.0)

    schemes  = [euler, heun, rk4]

    for scheme in schemes:
        z = scheme(f, z0, time)
        max_error = np.max(u_exact(time) - z[:,0])
        msg = '%s failed with error = %g' % (scheme.func_name, max_error)
        assert max_error < tol, msg

# f3 defines an ODE with ananlytical solution in u_nonlin_analytical
def f3(z, t, a=2.0, b=-1.0):
    """ """
    return a*z + b

def u_nonlin_analytical(u0, t, a=2.0, b=-1.0):
    from numpy import exp
    TOL = 1E-14
    if (abs(a)>TOL):
        return (u0 + b/a)*exp(a*t)-b/a
    else:
        return u0 + b*t
        
     
# Function for convergence test
def convergence_test():
    """ Test convergence rate of the methods """
    from numpy import linspace, size, abs, log10, mean, log2
    figure()
    tol = 1E-15
    T = 8.0   # end of simulation
    Ndts = 5 # Number of times to refine timestep in convergence test

    z0 = 2

    schemes =[euler, heun, rk4]
    legends=[]
    schemes_order={}
    
    colors = ['r', 'g', 'b', 'm', 'k', 'y', 'c']
    linestyles = ['-', '--', '-.', ':', 'v--', '*-.']
    iclr = 0
    for scheme in schemes:
        N = 30    # no of time steps
        time = linspace(0, T, N+1)

        order_approx = []
        
        for i in range(Ndts+1):
            z = scheme(f3, z0, time)   
            abs_error = abs(u_nonlin_analytical(z0, time)-z[:,0])
            log_error = log2(abs_error[1:]) # Drop 1st elt to avoid log2-problems (1st elt is zero)
            max_log_err = max(log_error)
            plot(time[1:], log_error, linestyles[i]+colors[iclr], markevery=N/5)
            legends.append(scheme.func_name +': N = ' + str(N))
            hold('on')
            
            if i > 0: # Compute the log2 error difference
                order_approx.append(previous_max_log_err - max_log_err) 
            previous_max_log_err = max_log_err

            N *=2
            time = linspace(0, T, N+1)
        
        schemes_order[scheme.func_name] = order_approx
        iclr += 1

    legend(legends, loc='best')
    xlabel('Time')
    ylabel('log(error)')
    grid()
    
    N = N/2**Ndts
    N_list = [N*2**i for i in range(1, Ndts+1)]
    N_list = np.asarray(N_list)
    
    figure()
    for key in schemes_order:
        plot(N_list, (np.asarray(schemes_order[key])))
    
    # Plot theoretical n for 1st, 2nd and 4th order schemes
    axhline(1.0, xmin=0, xmax=N, linestyle=':', color='k')
    axhline(2.0, xmin=0, xmax=N, linestyle=':', color='k')
    axhline(4.0, xmin=0, xmax=N, linestyle=':', color='k')
    xticks(N_list, rotation=-70)
    legends = schemes_order.keys()
    legends.append('theoretical') 
    legend(legends, loc='best', frameon=False)
    xlabel('Number of unknowns')
    ylabel('Scheme order approximation')
    axis([0, max(N_list), 0, 5])
    savefig('ConvergenceODEschemes.png', transparent=True)
    
def plot_ODEschemes_solutions():
    """Plot the solutions for the test schemes in schemes"""
    from numpy import linspace
    figure()
    T = 1.5  # end of simulation
    N = 50  # no of time steps
    time = linspace(0, T, N+1)

    z0 = 2.0

    schemes  = [euler, heun, rk4]
    legends = []

    for scheme in schemes:
        z = scheme(f3, z0, time)
        plot(time, z[:,-1])
        legends.append(scheme.func_name)

    plot(time, u_nonlin_analytical(z0, time))
    legends.append('analytical')
    legend(legends, loc='best', frameon=False)

