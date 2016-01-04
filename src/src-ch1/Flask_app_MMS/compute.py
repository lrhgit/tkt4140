from numpy import exp, cos, linspace, zeros_like
from ODEschemes import euler, heun, rk4
from math import sqrt, pi

import matplotlib.pyplot as plt
import os, time, glob, math
import numpy as np
from numpy import sin, pi
from numpy import exp as e

def compute_Manufactured_solution(Expression='(1/(sigma*sqrt(2*pi)))*exp(-((t-mu)**2)/(2*sigma**2))',
                                  sigma=0.5, mu=0.5, Domain='[-1.5,2.5]',
                                  solver='euler', Nsteps = 100, erase='No' 
                                  ):
    
    """Function that solves the equation f + f''*f + f''' = RHS numerically with different ODEsolvers with 
        use of the method of manufactured solution. The necesary Right hand side (RHS) is calculated from
        the differential equation by choosing/manufacturing a wanted solution. The differentiations are
        performed with use of the sympolic python module sympy. 
        In this case the chosen solution f=Expression is set default to be a normal distribution.
        In order to be able to handle various functions f/Ekspression on various domains and with different
        solvers, filled in on a form on a webpage  the arguments "Expression", "Domain" and "solver" are of 
        type str. and are evaluated into python expressions with use of the eval() function.
        The default values of the arguments are automatically filled into the form of "view.html" file.
        
        Args:
            Expression(str): the manufactured solution
            sigma(float): standard deviation
            mu(float): mean value
            Domain(str): the computational domain 
            solver(str): the chosen ODEsolver
            Nsteps(int): the number of timesteps
            erase(str): erase previous plot (Yes/No)
            
        Returns:
            figdata_png(string): the plot of (t, fnumeric) and (t, fanalytic) encoded with BytesIo to the fed into html file
        """
        
    from sympy import symbols, diff, integrate, Rational, lambdify, exp, sin, cos
    t = symbols('t')
    #### Turn strings into python expressions ####
    f = eval(Expression)
    [t0, tend] = eval(Domain)
    solver = eval(solver)
    #### Perform needed differentiations based on the differential equation ####
    dfdt = diff(f, t)
    d2fdt = diff(dfdt, t)
    d3fdt = diff(d2fdt, t)
    RHS = d3fdt + dfdt*d2fdt + f
    #### Create Python functions of f, RHS and needed differentiations of f ####
    f = lambdify([t], f)
    dfdt = lambdify([t], dfdt)
    d2fdt = lambdify([t], d2fdt)
    RHS = lambdify([t], RHS)
    
    #### Discretize time ####
    time = np.linspace(t0, tend, Nsteps + 1)
    
    def func(y, t):
        """ Function that returns the dfn/dt of the differential equation f + f''*f + f''' = RHS
            as a system of 1st order equation; f = f1 
                    f1' = f2 
                    f2' = f3
                    f3' = RHS - f1 - f2*f3
        
        Args:
            y(array): solutian array [f1, f2, f3] at time t
            t(time): current time

        Returns:
            yout(array): differantiation array [f1', f2', f3'] at time t
        """
        yout = np.zeros_like(y)
        yout[:] = [y[1], y[2], RHS(t) - y[0]- y[1]*y[2]]
        
        return yout
    
    #### Solve for fnumeric and evaluate fanalytic on the domain: ####
    y0 = np.array([f(t0), dfdt(t0), d2fdt(t0)]) #initial values calculated from the differentiations of f at t0
    u = solver(func, y0, time)
    fnumeric = u[:,0]
    fanalytic = np.zeros_like(fnumeric)
    i = 0
    for t in time:
        fanalytic[i] = f(t)
        i = i + 1
        
    if erase=='Yes':
        plt.figure()
        print "hello"
        
    #### Create plots and generate BytesIo plot data  ####
    plt.figure()
    legends = []
    plt.plot(time, fanalytic,'g')
    
    
    plt.plot(time, fnumeric, 'r--')
    analyticlegend = 'fanalytic'#, sigma = ' + str(sigma) + ', mu = ' + str(mu)
    numericlegend = 'fnumeric'#, solver = ' + solvername + ', Nsteps = ' + NstepsStr
    legends.append(analyticlegend)
    legends.append(numericlegend)
    plt.title('Normal distribution')
    plt.xlabel('t')
    plt.ylabel('f(t)')
    plt.legend(legends,loc='best',frameon=False)
    
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
    return figdata_png#, figdata_svg


