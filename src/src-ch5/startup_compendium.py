# src-ch5/startup.py;Visualization.py @ git@lrhgit/tkt4140/src/src-ch5/Visualization.py;Startupfuncs.py @ git@lrhgit/tkt4140/src/src-ch5/Startupfuncs.py;

import scipy 
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

import numpy as np

# Theta-scheme and using L'hopital for r=0
def thetaSchemeNumpyV1(theta, D, N, wOld):
    """ Algorithm for solving w^(n+1) for the startup of pipeflow 
        using the theta-schemes. L'hopitals method is used on the 
        governing differential equation for r=0.
        
        Args:
            theta(float): number between 0 and 1. 0->FTCS, 1/2->Crank, 1->Laasonen
            D(float): Numerical diffusion number [dt/(dr**2)]
            N(int): number of parts, or dr-spaces. In this case equal to the number of unknowns
            wOld(array): The entire solution vector for the previous timestep, n.
            
        Returns:
            wNew(array): solution at timestep n+1
    """
    
    superDiag = np.zeros(N - 1)
    subDiag = np.zeros(N - 1)
    mainDiag = np.zeros(N)
    
    RHS = np.zeros(N)
    
    j_array = np.linspace(0, N, N + 1)
    tmp = D*(1. - theta)
    
    superDiag[1:] = -D*theta*(1 + 0.5/j_array[1:-2])
    mainDiag[1:] = np.ones(N - 1)*(1 + 2*D*theta)
    subDiag[:] = -D*theta*(1 - 0.5/j_array[1:-1])
    
    a = tmp*(1 - 1./(2*j_array[1:-1]))*wOld[0:-2]
    b = (1 - 2*tmp)*wOld[1:-1]
    c = tmp*(1 + 1/(2*j_array[1:-1]))*wOld[2:]
    
    RHS[1:] = a + b + c
    
    superDiag[0] = -4*D*theta
    mainDiag[0] = 1 + 4*D*theta
    RHS[0] = (1 - 4*tmp)*wOld[0] + 4*tmp*wOld[1] 
    
    A = scipy.sparse.diags([subDiag, mainDiag, superDiag], [-1, 0, 1], format='csc')
    
    wNew = scipy.sparse.linalg.spsolve(A, RHS)
    wNew = np.append(wNew, 0)

    
    return wNew

# Theta-scheme and using 2nd order forward difference for r=0
def thetaScheme_numpy_V2(theta, D, N, wOld):
    """ Algorithm for solving w^(n+1) for the startup of pipeflow 
        using the theta-schemes. 2nd order forward difference is used  
        on the von-Neumann bc at r=0.
        
        Args:
            theta(float): number between 0 and 1. 0->FTCS, 1/2->Crank, 1->Laasonen
            D(float): Numerical diffusion number [dt/(dr**2)]
            N(int): number of parts, or dr-spaces.
            wOld(array): The entire solution vector for the previous timestep, n.
            
        Returns:
            wNew(array): solution at timestep n+1
    """
    superDiag = np.zeros(N - 2)
    subDiag = np.zeros(N - 2)
    mainDiag = np.zeros(N-1)
    
    RHS = np.zeros(N - 1)
    
    j_array = np.linspace(0, N, N + 1)
    tmp = D*(1. - theta)
    
    superDiag[1:] = -D*theta*(1 + 0.5/j_array[2:-2])
    mainDiag[1:] = np.ones(N - 2)*(1 + 2*D*theta)
    subDiag[:] = -D*theta*(1 - 0.5/j_array[2:-1])
    
    a = tmp*(1 - 1./(2*j_array[2:-1]))*wOld[1:-2]
    b = (1 - 2*tmp)*wOld[2:-1]
    c = tmp*(1 + 1/(2*j_array[2:-1]))*wOld[3:]
    
    RHS[1:] = a + b + c
    
    superDiag[0] = -(4./3)*D*theta
    mainDiag[0] = 1 + (4./3)*D*theta
    RHS[0] = (1 - (4./3)*tmp)*wOld[1] + (4./3)*tmp*wOld[2] 
    
    A = scipy.sparse.diags([subDiag, mainDiag, superDiag], [-1, 0, 1], format='csc')
    
    wNew = scipy.sparse.linalg.spsolve(A, RHS)
    w_0 = (1./3)*(4*wNew[0] - wNew[1])
    
    wNew = np.append(w_0, wNew)
    wNew = np.append(wNew, 0)

    
    return wNew



if __name__ == '__main__':
    
    import matplotlib.pylab as plt
    import time as timeModule
    from Visualization import createAnimation
    from Startupfuncs import analyticSolution

    # change some default values to make plots more readable 
    LNWDT=3; FNT=9
    plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT
        
    #### main program starts here ####
    
    N = 25 # No. of parts
    h = 1./N # length of r-step
    D = 0.4 # numerical diffusion number
    dt = D*h**2
    nmax = 1500 # No. of time-steps
    T = nmax*dt
    time = np.linspace(0, T, nmax + 1)
    r = np.linspace(0, 1, N + 1)
    eps = 2.2204*10**-16 #tollerance to be used in calculations of analytical solution

    theta = 0
    scheme_name = 'FTCS'
    
    wStart = 1 - r**2 # initial profile is the pouseille profile 
    wOld = wStart
    #initialize solution matrices
    uSolutions = np.zeros((1, len(time), N + 1))
    uAnalytic = np.zeros((len(time), N + 1))
    
    # solve numerically with vectorized solver
    
    tic = timeModule.time()
    for n in range(nmax):
        
        wNew = thetaSchemeNumpyV1(theta, D, N, wOld)
        uNew = 1 - r**2 - wNew
        
        
        uSolutions[0, n + 1, :] = uNew
        wOld = wNew
    
    toc = timeModule.time()
    print "vectorized scipy sparse solver CPU-time: ", toc - tic
    

    animate = True
    if animate:
        # solve analytically:
        tic = timeModule.time()
        for i in range(nmax):
             
            t = time[i + 1]
            for k, r_variable in enumerate(r):
                uAnalytic[i + 1, k] = 1 - r_variable**2 - analyticSolution(r_variable, t, eps)
     
        toc = timeModule.time()
        print "analytic solution CPU-time: ", toc - tic

        createAnimation(uSolutions, uAnalytic, [scheme_name], r, time, jump = int(round(0.002/dt)))
        

    
        

    
    

