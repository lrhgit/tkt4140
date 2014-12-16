import numpy as np
from matplotlib.pyplot import plot, show, legend


def euler(func,z0, time):
    """The Euler scheme for solution of systems of of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""
    
    z = np.zeros((np.size(time),2))
    z[0,:] = z0

    for i in range(len(time)-1):
        dt = time[i+1]-time[i] 
        z[i+1,:]=z[i,:] + np.asarray(func(z[i,:],time[i]))*dt

    return z

def euler2(func,z0, time):
    """The Euler scheme for solution of systems of of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""

    def f_np(z,t):
        return np.asarray(func(z,t))

    z = np.zeros((len(time),2))
    z[0,:] = z0

    for i, t in enumerate(time[0:-1]):
        dt = time[i+1] - t
        z[i+1,:]=z[i,:] + f_np(z[i,:],t)*dt

    return z

def euler3(func,z0, time):
    """The Euler scheme for solution of systems of of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""
    
    time_local = np.asarray(time)
    z = np.zeros((np.size(time_local),2))
    z[0,:] = float(z0)

    for i, t in enumerate(time_local[0:-1]):
        dt = time_local[i+1] - t
        z[i+1,:]=z[i,:] + np.asarray(func(z[i,:],t))*dt

    return z

def euler4(func,z0, time):
    """The Euler scheme for solution of systems of of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""
    
    time_local = np.asarray(time)
    z = np.zeros((np.size(time_local),2))
    z[0,:] = float(z0)

    for i, t in enumerate(time_local[1:]):
        dt = t-time_local[i] 
        z[i+1,:]=z[i,:] + np.asarray(func(z[i,:],time_local[i]))*dt

    return z


def heun(func,z0, time):
    """The Heun scheme for solution of systems of of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""
    
    def f_np(z,t):
        """A local function to ensure that the return of func is an np array 
        and to avoid lengthy code for implementation of the Heun algorithm"""
        return np.asarray(func(z,t))
    
    z = np.zeros((np.size(time),2))
    z[0,:] = z0
    zp = np.zeros_like(z0)
    
    for i, t in enumerate(time[0:-1]):
        dt = time[i+1]-time[i]
        zp = z[i,:] + f_np(z[i,:],t)*dt   # Predictor step
        z[i+1,:] = z[i,:] + (f_np(z[i,:],t) + f_np(zp,t+dt))*dt/2.0 # Corrector step

    return z

def heun2(func,z0, time):
    """The Heun scheme for solution of systems of of ODEs. 
    z0 is a vector for the initial conditions, 
    the right hand side of the system is represented by func which returns 
    a vector with the same size as z0 ."""
    
    z = np.zeros((np.size(time),2))
    z[0,:] = z0
    zp = np.zeros_like(z0)
    
    for i, t in enumerate(time[0:-1]):
        dt = time[i+1]-time[i]
        zp = z[i,:] + np.asarray(func(z[i,:],t))*dt   # Predictor step
        z[i+1,:] = z[i,:] + (np.asarray(func(z[i,:],t)) + np.asarray(func(zp,t+dt)))*dt/2.0 # Corrector step

    return z




if __name__ == '__main__':              
#Check whether this file is executed (name==main) or imported as module

    a = 0.2
    b = 3.0 
    u_exact = lambda t: a*t   +  b 
     
    def f_local(u,t):
        """A function which returns an np.array but less easy to read
        than f(z,t) below. """
        return np.asarray([a + (u - u_exact(t))**5])

    def f(z, t):
        """Simple to read function implementation """
        return [a + (z - u_exact(t))**5]
    
    
    def test_ode_schemes():
        """Use knowledge of an exact numerical solution for testing."""
        from numpy import linspace
        
        tol = 1E-15
        T = 2.0  # end of simulation
        N = 20  # no of time steps
        time = linspace(0, T, N+1)
    
        
        z0=np.zeros(1)
        z0[0] = u_exact(0.0)
        
        scheme_list  = [euler, euler2, euler3, euler4, heun, heun2]
         
        for scheme in scheme_list:
            z = scheme(f,z0,time)
            max_error = np.max(u_exact(time)-z[:,1])
            msg = '%s failed with error = %g' % (scheme.func_name, max_error)
            assert max_error < tol, msg
             
    def plot_ode_schemes_solutions():
        """Plot the linear solutions for the test schemes in scheme_list"""
        from numpy import linspace
        T = 1.5  # end of simulation
        N = 5  # no of time steps
        time = linspace(0, T, N+1)

        z0 = np.zeros(1)
        z0[0] = u_exact(0.0)
        
        scheme_list  = [euler, heun]
        legends = []
        
         
        for scheme in scheme_list:
            z = scheme(f_local,z0,time)
            plot(time,z[:,-1])
            legends.append(scheme.func_name)
             
             
        legend(legends) 
        show()
        
    test_ode_schemes()
    #plot_ode_schemes_solutions()               

        
        