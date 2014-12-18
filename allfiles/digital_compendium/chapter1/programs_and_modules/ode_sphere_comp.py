import numpy as np

from ode_schemes import euler, heun

if __name__ == '__main__':              
#Check whether this file is executed (name==main) or imported as module

    import odespy
    import matplotlib
    from matplotlib.pyplot import legend, plot, show

    LNWDT=5; FNT=11
    matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT
    
    g = 9.81      # Gravity m/s^2
    d = 41.0e-3     # Diameter of the sphere
    rho_f = 1.22  # Density of fluid [kg/m^3]
    rho_s = 1275  # Density of sphere [kg/m^3]
    nu = 1.5e-5   # Kinematical viscosity [m^2/s]

    def f(z, t):
        """2x2 syst for sphere with constant drag."""
        zout = np.zeros_like(z)
        CD = 0.5
        alpha = 3.0*rho_f/(4.0*rho_s*d)*CD
        zout[:] = [z[1], g - alpha*z[1]**2]
        return zout 

    
                
    # Main program starts here
    from numpy import linspace
    T = 30  # end of simulation
    N = 25  # no of time steps
    time = linspace(0, T, N+1)
    
    solvers=[]
    solvers.append(odespy.RK3(f)) 
    solvers.append(odespy.RK4(f)) 
    
    legends=[]
    
    z0=np.zeros(2)
    z0[0] = 2.0
    
    for i, solver in enumerate(solvers):
        solver.set_initial_condition(z0)
        z, t = solver.solve(time)
        plot(t,z[:,1])
        legends.append(str(solver))
    
    
    scheme_list  = [euler, heun]
    
    for scheme in scheme_list:
        z = scheme(f,z0,time)
        plot(time,z[:,1])
        legends.append(scheme.func_name)
        
    
    legend(legends, loc='best', frameon=False)
    
    show()
