# src-ch5/couette_Flow_FTCS.py;Visualization.py @ git@lrhgit/tkt4140/src/src-ch5/Visualization.py;


import numpy as np
from math import exp, sin, pi


def analyticSolution(y, t, N=100):
    
    """ Method that calculates the analytical solution to the differential equation:
        du/dt = d^2(u)/dx^2 , u = u(y,t), 0 < y < 1
        Boundary conditions: u(0, t) = 1, u(1, t) = 0
        Initial condition: u(t, 0) = 0 t<0,  u(t, 0) = 1 t>0
            
        Args:
            y(float): radial coordinat
            t(float): time
            N(int): truncation integer. Truncate sumation after N elements

    
        Returns:
            w(float): velocity, us - ur
    """
    sumValue = 0
    for n in range(1,N+1):
        temp = exp(-t*(n*pi)**2)*sin(n*pi*y)/n
        sumValue += temp
    u = 1 - y - (2/pi)*sumValue
    return u


def solveNextTimestepFTCS(Uold, D):
    """ Method that solves the transient couetteflow using the FTCS-scheme..
        At time t=t0 the plate starts moving at y=0
        The method solves only for the next time-step.
        The Governing equation is:
        
        du/dt = d^2(u)/dx^2 , u = u(y,t), 0 < y < 1
        
        Boundary conditions: u(0, t) = 1, u(1, t) = 0
        
        Initial condition: u(t, 0) = 0 t<0,  u(t, 0) = 1 t>0
        
        Args:
            uold(array): solution from previous iteration
            D(float): Numerical diffusion number
            
        Returns:
            unew(array): solution at time t^n+1
    """
    Unew = np.zeros_like(Uold)
    
    Uold_plus = Uold[2:]
    Uold_minus = Uold[:-2]
    Uold_mid = Uold[1:-1]
    
    Unew[1:-1] = D*(Uold_plus + Uold_minus) + (1 - 2*D)*Uold_mid
    Unew[0] = 1
    
    return Unew



            
if __name__ == '__main__':
    
    import numpy as np
    from Visualization import createAnimation
    
    D = 0.4
    
    N = 20
    Y = np.linspace(0, 1, N + 1)
    h = Y[1] - Y[0]
    dt = D*h**2 # numerical diffusion number
    T = 0.2 # simulation time
    time = np.arange(0, T + dt, dt)
    
    
    # solution matrices:
    Usolutions = np.zeros((len(time), N + 1))
    Usolutions[0, 0] = 1 # no slip condition at the plate boundary
    
    Uanalytic = np.zeros((len(time), N + 1))
    Uanalytic[0, 0] = 1
    
    
    for n, t in enumerate(time[1:]):
        
        Uold = Usolutions[n, :]
        Unew = solveNextTimestepFTCS(Uold, D)
        
        Unew_analytic = np.zeros_like(Y)
        
        for i, y in enumerate(Y):
            Unew_analytic[i] = analyticSolution(y, t)
            
        Usolutions[n + 1, :] = Unew
        Uanalytic[n + 1, :] = Unew_analytic
        
    Usolutions_Visualization = np.zeros((1, len(time), N + 1))
    Usolutions_Visualization[0, :, :] = Usolutions
    
    createAnimation(Usolutions_Visualization, Uanalytic, ["FTCS"], Y, time, symmetric=False)

    
    

