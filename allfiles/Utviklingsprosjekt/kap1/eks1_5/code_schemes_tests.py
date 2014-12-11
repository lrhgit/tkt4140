import numpy as np
from matplotlib.pyplot import plot, show, legend

def f(z, t):
    """2x2 syst for sphere with constant drag."""
    zout = np.zeros_like(z)
    CD = 0.5
    alpha = 3.0*rho_f/(4.0*rho_s*d)*CD
    zout[:] = [z[1], g - alpha*z[1]**2]
    return zout 


def euler(func,z0, time):
    """The Euler scheme for solution of systems of of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""
    
    dt = time[1]-time[0]
    z = np.zeros((np.size(time),2))
    z[0,:] = z0

    for i, t in enumerate(time[1:]):
        z[i+1,:]=z[i,:] + func(z[i,:],t)*dt

    return z

def heun(func,z0, time):
    """The Heun scheme for solution of systems of of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""
    
    dt = time[1]-time[0]
    z = np.zeros((np.size(time),2))
    z[0,:] = z0
    zp = np.zeros_like(z0)
    
    for i, t in enumerate(time[1:]):
        zp = z[i,:] + func(z[i,:],t)*dt   # Predictor step
        z[i+1,:] = z[i,:] + (func(z[i,:],t) + func(zp,t+dt))*dt/2.0 # Corrector step

    return z


class ForwardEuler:
    """
    Class for solving an ODE,

      du/dt = f(u, t)

    by the ForwardEuler solver.

    Class attributes:
    t: array of time values
    u: array of solution values (at time points t)
    k: step number of the most recently computed solution
    f: callable object implementing f(u, t)
    dt: time step (assumed constant)
    """
    def __init__(self, f):
        if not callable(f):
            raise TypeError('f is %s, not a function' % type(f))
        self.f = f

    def set_initial_condition(self, U0):
        self.U0 = float(U0)

    
    def solve(self, time_points):
        """Compute u for t values in time_points list."""
        self.t = np.asarray(time_points)
        self.u = np.zeros(len(time_points))
        # Assume self.t[0] corresponds to self.U0
        self.u[0] = self.U0

        for k in range(len(self.t)-1):
            self.k = k
            self.u[k+1] = self.advance()
        return self.u, self.t

    def advance(self):
        """Advance the solution one time step."""
        # Load attributes into local variables to
        # obtain a formula that is as close as possible
        # to the mathematical notation.
        u, f, k, t = self.u, self.f, self.k, self.t

        dt = t[k+1] - t[k]
        u_new = u[k] + dt*f(u[k], t[k])
        return u_new


if __name__ == '__main__':              
#Check whether this file is executed (name==main) or imported as module


    
    def test_ode_schemes():
        """Use knowledge of an     exact numerical solution for testing."""
        from numpy import linspace
        T = 2.0  # end of simulation
        N = 20  # no of time steps
        time = linspace(0, T, N+1)

        a = 0.2
        b = 3.0 
        u_exact = lambda t: a*t   +  b 
         
#         def u_exact(t):
#            return 0.2*t + 3
#         
        def f_local(u,t):
            return np.asarray([0.2 + (u - 0.2*t -3)**5])
            # return np.asarray([a])


        def f(z, t):
            zout = np.zeros_like(z)
            zout[0] = np.asarray([0.2 + (z[0] - u_exact(t))**5])
            
            return zout 

        
        z0=np.zeros(1)
        z0[0] = u_exact(0.0)
        
#        euler(f_local,z0, time)

        # scheme_list  = [euler, heun]
        scheme_list  = [heun]
         
        for scheme in scheme_list:
             z = scheme(f_local,z0,time)
             
        #solver = ForwardEuler(lambda u, t: 0.2 + (u - u_exact(t))**4)
        solver = ForwardEuler(f)

        # Solve for first time interval [0, 1.2]
        solver.set_initial_condition(u_exact(0))

        zc, tc = solver.solve(time)
        
        #print 'success', np.abs(np.asarray(u_exact(time))-z[:,-1]).max(), np.abs(u_exact(time)-zc).max
        print 'success', np.max(u_exact(time)-zc[:])
        
    def plot_ode_schemes_solutions():
        """Use knowledge of an     exact numerical solution for testing."""
        from numpy import linspace
        T = 1.5  # end of simulation
        N = 5  # no of time steps
        time = linspace(0, T, N+1)

        a = 0.2
        b = 3.0 
        u_exact = lambda t: a*t   +  b 
         
#         def u_exact(t):
#            return 0.2*t + 3
#         
        def f_local(u,t):
            return np.asarray([0.2 + (u - 0.2*t -3)**5])
         

        def f(z, t):
            zout = np.zeros_like(z,dtype=np.float64)
            zout[0] = np.asarray([0.2 + (z[0] - u_exact(t))**5])
            
            return zout 

        
        z0=np.zeros(1)
        z0[0] = u_exact(0.0)
        
        euler(f_local,z0, time)

        scheme_list  = [euler, heun]
        legends = []
        
         
        for scheme in scheme_list:
             z = scheme(f_local,z0,time)
             plot(time,z[:,-1])
             legends.append(scheme.func_name)
             
             
        legend(legends) 
        show()
        
test_ode_schemes()
plot_ode_schemes_solutions()               

        
        