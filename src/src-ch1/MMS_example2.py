# src-ch1/MMS_example2.py; ODEschemes.py @ git@lrhgit/tkt4140/src/src-ch1/ODEschemes.py;
 
import numpy as np
import matplotlib.pylab as plt
from ODEschemes import euler, heun, rk4
from sympy import exp, symbols, diff, lambdify
from math import sqrt, pi

# Newton solver
def Newton_solver_sympy(error, h, x0):
    """ Function that solves for the nonlinear set of equations
        error1 = C*h1^p --> f1 = C*h1^p - error1 = 0
        error2 = C*h2^p --> f2 = C h2^p - error 2 = 0
        where C is a constant h is the step length and p is the order,
        with use of a newton rhapson solver. In this case C and p are
        the unknowns, whereas h and error are knowns. The newton rhapson 
        method is an iterative solver which take the form:
        xnew = xold - (J^-1)*F, where J is the Jacobi matrix and F is the 
        residual funcion.
            x = [C, p]^T
            J = [[df1/dx1  df2/dx2],
                 [df2/dx1  df2/dx2]]
            F = [f1, f2]
        This is very neatly done with use of the sympy module
        
        Args:
            error(list): list of calculated errors [error(h1), error(h2)]
            h(list): list of steplengths corresponding to the list of errors
            x0(list): list of starting (guessed) values for x
        
        Returns:
            x(array): iterated solution of x = [C, p]

    """
    from sympy import Matrix
    #### Symbolic computiations: ####
    C, p = symbols('C p')
    f1 = C*h[-2]**p - error[-2]
    f2 = C*h[-1]**p - error[-1]
    F = [f1, f2]
    x = [C, p]
    
    def jacobiElement(i,j):
        return diff(F[i], x[j])
        
    Jacobi = Matrix(2, 2, jacobiElement) # neat way of computing the Jacobi Matrix
    JacobiInv = Jacobi.inv()
    #### Numerical computations: ####
    JacobiInvfunc = lambdify([x], JacobiInv)
    Ffunc = lambdify([x], F)
    x = x0
    
    for n in range(8): #perform 8 iterations
        F = np.asarray(Ffunc(x))
        Jinv = np.asarray(JacobiInvfunc(x))
        xnew = x - np.dot(Jinv, F)
        x = xnew
#        print "n, x: ", n, x, F
    x[0] = round(x[0], 2)
    x[1] = round(x[1], 3)
 
    return x

# Differential function
def func(y,t):
    """ Function that returns the dfn/dt of the differential equation f''' + f''*f + f' = RHS
        as a system of 1st order equations; f = f0 
                f0' = f1 
                f1' = f2
                f2' = RHS - f2*f0 - f1 
    
        Args:
            y(array): solutian array [f0, f1, f2] at time t
            t(float): current time

        Returns:
            yout(array): differantiation array [f0', f1', f2'] at time t
    """
    yout = np.zeros_like(y)
    yout[:] = [y[1], y[2], RHS(t) -y[1]- y[0]*y[2]]
    
    return yout

#### Main Program starts here ####

t = symbols('t')
sigma=  0.5 # standard deviation
mu = 0.5 # mean value

#### Perform needed differentiations based on the differential equation ####
f = (1/(sigma*sqrt(2*pi)))*exp(-((t-mu)**2)/(2*sigma**2))
dfdt = diff(f, t)
d2fdt = diff(dfdt, t)
d3fdt = diff(d2fdt, t)
RHS = d3fdt + f*d2fdt + dfdt

#### Create Python functions of f, RHS and needed differentiations of f ####
f = lambdify([t], f, np)
dfdt = lambdify([t], dfdt, np)
d2fdt = lambdify([t], d2fdt)
RHS = lambdify([t], RHS)

t0, tend = -1.5, 2.5
z0 = np.array([f(t0), dfdt(t0), d2fdt(t0)]) # initial values

schemes = [euler, heun, rk4] # list of schemes; each of which is a function
schemes_error = {} # empty dictionary. 
h = [] # empty list of time step

Ntds = 8 # number of times to refine dt
 
for k, scheme in enumerate(schemes):
    N = 20    # initial number of time steps
    error = [] # start of with empty list of errors for all schemes
    legendList = [] 
    
    for i in range(Ntds + 1):
        time = np.linspace(t0, tend, N+1)
        
        if k==0:
            h.append(time[1] - time[0]) # add this iteration's dt to list h
        z = scheme(func, z0, time) # e.g: euler(func, z0, time) 
        fanalytic = f(time) # call analytic function f 
        
        abs_error = np.abs(z[:,0]- fanalytic) # calculate infinity norm
        error.append(max(abs_error))
        
        N *=2 # refine dt
    
    schemes_error[scheme.func_name] = error # Add a key:value pair to the dictionary


ht = np.asarray(h) 
eulerError = np.asarray(schemes_error["euler"])
heunError = np.asarray(schemes_error["heun"])
rk4Error = np.asarray(schemes_error["rk4"])

[C_euler, p_euler] = Newton_solver_sympy(eulerError, ht, [1,1])
[C_heun, p_heun] = Newton_solver_sympy(heunError, ht, [1,2])
[C_rk4, p_rk4] = Newton_solver_sympy(rk4Error, ht, [1,4])

print C_euler, p_euler
print C_heun, p_heun
print C_rk4, p_rk4

from sympy import latex
h = symbols('h')
epsilon_euler = C_euler*h**p_euler
epsilon_euler_latex = '$' + latex(epsilon_euler) + '$'
epsilon_heun = C_heun*h**p_heun
epsilon_heun_latex = '$' + latex(epsilon_heun) + '$'
epsilon_rk4 = C_rk4*h**p_rk4
epsilon_rk4_latex = '$' + latex(epsilon_rk4) + '$'

epsilon_euler = lambdify(h, epsilon_euler, np)
epsilon_heun = lambdify(h, epsilon_heun, np)
epsilon_rk4 = lambdify(h, epsilon_rk4, np)

N = N/2**(Ntds + 2)
N_list = [N*2**i for i in range(1, Ntds + 2)]
N_list = np.asarray(N_list)

plt.figure()

# plt.plot(N_list, np.log2(eulerError), 'b')
# plt.plot(N_list, np.log2(epsilon_euler(ht)), 'b--')
# plt.plot(N_list, np.log2(heunError), 'g')
# plt.plot(N_list, np.log2(epsilon_heun(ht)), 'g--')
# plt.plot(N_list, np.log2(rk4Error), 'r')
# plt.plot(N_list, np.log2(epsilon_rk4(ht)), 'r--')


plt.plot(np.log2(N_list), np.log2(eulerError), 'b')
plt.plot(np.log2(N_list), np.log2(epsilon_euler(ht)), 'b--')
plt.plot(np.log2(N_list), np.log2(heunError), 'g')
plt.plot(np.log2(N_list), np.log2(epsilon_heun(ht)), 'g--')
plt.plot(np.log2(N_list), np.log2(rk4Error), 'r')
plt.plot(np.log2(N_list), np.log2(epsilon_rk4(ht)), 'r--')


LegendList = ['${\epsilon}_{Euler}$', epsilon_euler_latex, '${\epsilon}_{Heun}$', epsilon_heun_latex, '${\epsilon}_{RK4}$', epsilon_rk4_latex]
plt.legend(LegendList, loc='best', frameon=False,fontsize='20')
plt.xlabel('log2(N)')
plt.ylabel('log2($\epsilon$)')
#plt.savefig('../fig-ch1/MMS_example2.png')
plt.show()
