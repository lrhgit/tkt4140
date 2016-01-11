# src/src-ch5/startupPythonic.py;TRIdiagonalSolvers.py @ git@lrhgit/tkt4140/src/src-ch5/TRIdiagonalSolvers.py;


import numpy as np
from scipy.special import jv as besselj
from TRIdiagonalSolvers import tdma


def analyticSolution(r,t,eps):
    
    """ Method that calculates the analytical solution to the differential equation:
            dw/dt = w'' + w'/r, w = w(r,t), 0 < r < 1.
            The accuracy can be controlled by changing the value of eps, which depicts the
            relative accuracy.
            
        Args:
            r(float): radial coordinat
            t(float): time
            eps(float): numerical tollerance

    
        Returns:
            w(float): velocity, us - ur
    """
    
    if r == 1:
        w = 0
    elif r < eps and t < eps:
        w = 1
    else:
        # calculate first element seperately
        n = 1
        lam1 = j0zero(n)
        arg1 = r*lam1
        term1 = np.exp(-t*lam1**2)*besselj(0,arg1)/(besselj(1,lam1)*lam1**3)
        sum1 = term1
        sumtol = 1*10**-8
        test = 1
        
        while test > sumtol:
            n = n + 1
            lamn = j0zero(n)
            arg = r*lamn
            term = np.exp(-t*lamn**2)*besselj(0,arg)/(besselj(1,lamn)*lamn**3)
            sum1 = sum1 + term
            test = np.abs(term/term1)
        
        w = 8 * sum1
    
    return w
        
        
def j0zero(s):
    
    """ method that computes root number s, s = 1,2,...,
        of the Besselfunction J0 where zj0 is the root,
        using table-values, asymptotic formulae and Newton-Raphsons method.
        Rellative error app. 1.0e-14 to 1.0e-15
        Args:
            s(integer): solutian array from last iteration

        Returns:
            zj0(float): root No. s of besselfunction
    """
    
    sz = np.array([2.4048225557695773, 5.520078110286311, 8.653727912911012, 11.79153443901428, 14.93091770848779])
    b0 = (s - 0.25)*np.pi
    b08 = 0.125/b0
    b082 = b08**2
    z0 = b0 + b08*(1. - b082*(124./3 -  b082*120928./15))
    
    if s <= 5:
        zj0 = sz[s - 1]
    elif s > 30:
        zj0 = z0
    else:
        dz0 = besselj(0, z0)/besselj(1,  z0)
        z0 = z0 + dz0
        zj0 = z0
    return zj0


def solveNextTimestepStartup(theta, D, N, wold):
    """ Method that sets up the tridiagonal theta-scheme for the startup of flow in a pipe (Szymanski's problem).
        The method solves only for the next time-step, using standard tdma solver.
        The Governing equation is:
        dw/dt = w'' + w'/r, w = w(r,t), 0 < r < 1
        where w' denotes dw/dr.
        velocity field: u(r,t) = us(r) -w(r,t)
        where us = 1 - r^2
        Boundary conditions: w^n+1(0) = (1-2D)w^n(1) + 4Du^n(1), w(r) = 0, w(r) = 0
        Initial condition: w(r,0) = 1- r^2 = us
        Args:
            theta(float): parameter between 0 and 1 seperating different schemes (FTCS, Laasonen, Crank-Nicolson...)
            D(float): Numerical diffusion number
            N(int): number of equations
            wold(array): solution from previous iteration

    
        Returns:
            wnew(array): solution at time t^n+1
    """
    tmp1 = D*(1. - theta)
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    for j in range(1, N):
        fac = 0.5 / (j)
        a[j] = -D*theta*(1 - fac)
        b[j] = (1 + 2*D*theta)
        c[j] = -D*theta*(1 + fac)
        d[j] = tmp1*((1 - fac)*wold[j - 1] + (1 + fac)*wold[j + 1]) + (1 - 2*tmp1)*wold[j]
        
    b[0] = 1 + 4*D*theta
    c[0] = -4*D*theta
    d[0] = wold[0] + 4*tmp1*(wold[1] - wold[0])
    wnew = tdma(a, b, c, d)
    wnew = np.append(wnew, [0]) # add boundary
    
    return wnew

def solveAllTimestepsStartup(theta, D, nmax, N, wold):
    """ Method that sets up the tridiagonal theta-scheme for the startup of flow in a pipe (Szymanski's problem).
        The method solves without storing solutions from t0 to t = dt*nmax, using standard tdma solver.
        The Governing equation is:
        dw/dt = w'' + w'/r, w = w(r,t), 0 < r < 1
        where w' denotes dw/dr.
        velocity field: u(r,t) = us(r) -w(r,t)
        where us = 1 - r^2
        Boundary conditions: w^n+1(0) = (1-2D)w^n(1) + 4Du^n(1), w(r) = 0
        Initial condition: w(r,0) = 1- r^2 = us
        Args:
            theta(float): parameter between 0 and 1 seperating different schemes (FTCS, Laasonen, Crank-Nicolson...)
            D(float): Numerical diffusion number
            nmax(int): number of iterations
            N(int): number of equations
            wold(array): solution from previous iteration

    
        Returns:
            wnew(array): solution at time t = dt*nmax
    """
    
    tmp1 = D*(1. - theta)
    
    for n in range(nmax):
        a = np.zeros(N)
        b = np.zeros(N)
        c = np.zeros(N)
        d = np.zeros(N)

        for j in range(1, N):
            fac = 0.5 / (j)
            a[j] = -D*theta*(1 - fac)
            b[j] = (1 + 2*D*theta)
            c[j] = -D*theta*(1 + fac)
            d[j] = tmp1*((1 - fac)*wold[j - 1] + (1 + fac)*wold[j + 1]) + (1 - 2*tmp1)*wold[j]
            
        b[0] = 1 + 4*D*theta
        c[0] = -4*D*theta
        d[0] = wold[0] + 4*tmp1*(wold[1] - wold[0])
        wnew = tdma(a, b, c, d)
        wold = np.append(wnew, [0])
    
    return wold

def solveNextTimestepStartupV2(theta, D, N, wold):
    """ Method that sets up the tridiagonal theta-scheme for the startup of flow in a pipe (Szymanski's problem).
        The method solves only for the next time-step, using standard tdma solver.
        The Governing equation is:
        dw/dt = w'' + w'/r, w = w(r,t), 0 < r < 1
        where w' denotes dw/dr.
        velocity field: u(r,t) = us(r) -w(r,t)
        where us = 1 - r^2
        Boundary conditions: w(0) = (4w(1)-w(2))/3), w(r) = 0
        Initial condition: w(r,0) = 1- r^2 = us
        Args:
            theta(float): parameter between 0 and 1 seperating different schemes (FTCS, Laasonen, Crank-Nicolson...)
            D(float): Numerical diffusion number
            N(int): number of equations
            wold(array): solution from previous iteration

    
        Returns:
            wnew(array): solution at time t^n+1
    """
    wold = wold[1:] # in this version of startup w0 is not solved with tdma
    N = N - 1 # this version has one less unknown to be solved with tdma
    tmp1 = D*(1. - theta)
    
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)

    for j in range(1,N):
        fac = 0.5 / (j + 1)
        a[j] = -D*theta*(1 - fac)
        b[j] = (1 + 2*D*theta)
        c[j] = -D*theta*(1 + fac)
        d[j] = tmp1*((1 - fac)*wold[j - 1] + (1 + fac)*wold[j + 1]) + (1 - 2*tmp1)*wold[j]
    b[0] = 1 + 4*D*theta/3
    c[0] = -4*D*theta/3
    d[0] = (1-4*tmp1/3.)*wold[0] + 4*tmp1*wold[1]/3.
    wnew = tdma(a, b, c, d) # solve with tdma
    w0 = (4*wnew[0]-wnew[1])/3
    wnew = np.append(np.append(w0, wnew), [0]) # add BCs
    
    return wnew

def solveAllTimestepsStartupV3(theta, D, nmax, N, wold):
    """ Method that sets up the tridiagonal theta-scheme for the startup of flow in a pipe (Szymanski's problem).
        The method solves without storing solutions from t0 to t = dt*nmax, using a modified tdma solver
        (possible since only the RHS is time dependent).
        The Governing equation is:
        dw/dt = w'' + w'/r, w = w(r,t), 0 < r < 1
        where w' denotes dw/dr.
        velocity field: u(r,t) = us(r) -w(r,t)
        where us = 1 - r^2
        Boundary conditions: w^n+1(0) = (1-2D)w^n(1) + 4Du^n(1), w(r) = 0
        Initial condition: w(r,0) = 1- r^2 = us
        Args:
            theta(float): parameter between 0 and 1 seperating different schemes (FTCS, Laasonen, Crank-Nicolson...)
            D(float): Numerical diffusion number
            nmax(int): number of iterations
            N(int): number of equations
            wold(array): solution from previous iteration

    
        Returns:
            wnew(array): solution at time t = dt*nmax
    """    
    
    tmp1 = D*(1. - theta)
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    wnew = np.zeros(N)

    for j in range(1, N):
        fac = 0.5 / (j)
        a[j] = -D*theta*(1 - fac)
        b[j] = (1 + 2*D*theta)
        c[j] = -D*theta*(1 + fac)
        
    b[0] = 1 + 4*D*theta
    c[0] = -4*D*theta
    # === elimination ===
    c[0] = c[0]/b[0]
    
    for j in range(1, N):
        b[j] = b[j] - a[j]*c[j - 1]
        c[j] = c[j]/b[j]
    
    for n in range(nmax):
        # === compute right hand side ===
        for j in range(1, N):
            fac = 0.5 / (j)
            d[j] = tmp1 *((1 - fac)*wold[j - 1] + (1 + fac)*wold[j+1]) + (1 - 2*tmp1)*wold[j]
        
        
        d[0] = wold[0] + 4*tmp1*(wold[1] - wold[0])

        # === backsubstitution ===
        d[0] = d[0]/b[0]
        for j in range(1, N):
            d[j] = (d[j] - a[j]*d[j - 1])/b[j]
        
        wnew[N-1] = d[N-1]
        
        for j in range(N - 2, -1, -1):
            wnew[j] = d[j] - c[j]*wnew[j + 1]
        
        wold = np.append(wnew,[0])
    
    return wold
            
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import time as timeModule
    from math import sqrt
    from Visualization import createAnimation, plotErrorAndOrder
    
    # change some default values to make plots more readable 
    LNWDT=3; FNT=9
    plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT
    
    def test_MES():
        """Code verification using the Method of Exact solution"""
        N = 25 # No. of parts
        h = 1./N # length of r-step
        D = 0.4 # numerical diffusion number
        dt = D*h**2
        nmax = 1500 # No. of time-steps
        T = nmax*dt
        time = np.linspace(0, T, nmax + 1)
        r = np.linspace(0, 1, N +1)
        
        schemes_theta = [0, 0.5, 1]
        schemes_name = ['FTCS', 'Crank-Nicolson', 'Laasonen']
        schemes_global_rms_error = []
        wstart = 1 - r**2 # initial profile is the pouseille profil
        #initialize solution matrices
        wSolutions = np.zeros((len(schemes_theta),len(time), N + 1))
        wAnalytic = np.zeros((len(time), N + 1))
        print """\n######## Started test 'test_MES()': #########
         
        Solving startup of a flow in a pipe (Szymanski's problem) using theta schemes:
        {0}, {1}, {2}""".format(schemes_name[0], schemes_name[1], schemes_name[2])
        print """        Numerical discretisation values:
        h = {0}, dt = {1}, D = {2}, t0 = 0 to T = {3}""".format(h, dt, D, T)
        print """        Acceptance criteria: Percent error, Visual inspection\n"""
        startTime = timeModule.time()
        
        #iterate all theta-schemes, solve them with solveNextTimestep tdma solver, and store solutions
        eps = 2.2204*10**-16 #tollerance to be used in calculations of analytical solution
        for scheme, theta in enumerate(schemes_theta):
            wold = wstart
            for i, t in enumerate(time[1:]):
                wnew = solveNextTimestepStartup(theta, D, N, wold)
                wold = wnew
                wSolutions[scheme, i + 1, :] = 1 - r**2 - wnew
                
                if scheme==0:
                    for k, r_variable in enumerate(r):
                        wAnalytic[i + 1, k] = 1 - r_variable**2 - analyticSolution(r_variable, t, eps)
            
            Nelements = (len(time)-1)*(N)
            error = sqrt(np.sum(((wSolutions[scheme, 1:, :-1] - wAnalytic[1:,:-1])/wAnalytic[1:,:-1])**2)/Nelements)
            schemes_global_rms_error.append(error)
        
        endTime = timeModule.time()
        print """        Finished solving all schemes and computed corresponding analytical solution.
        CPU time: {0}""".format(endTime - startTime), """ seconds"""
        print """        Calculated Percent Gloabel RMS-errors: 
        
        FTCS: {0},
        Crank-Nicolson: {1}, 
        Laasonen: {1}
        """.format(schemes_global_rms_error[0]*100, schemes_global_rms_error[1]*100, schemes_global_rms_error[2]*100)
        
        animate ='y'#= raw_input("        Do you want to animate resuts ? (y/n) ")
        
        if animate=='y':
            createAnimation(wSolutions, wAnalytic, schemes_name, r, time, jump = int(round(0.002/dt)))
        
                #plot solutions at different times
        nsnapshotlist=np.linspace(0, nmax, 8)

        lstyle = ['r-', 'b:', 'c.', 'y-.']
        legendList = []
        titlestring = 'test_MES: Profile at times t = '
        plt.figure()
        for n in nsnapshotlist:
            titlestring = titlestring + str(round(n*dt, 3))+ ', '
            for i, scheme_name in enumerate(schemes_name):
                if i == 0:
                    plt.plot(r, wAnalytic[n, :], lstyle[i])
                    
                    plt.plot(r, wSolutions[i, n, :], lstyle[i + 1])
                    if n==0:
                        legendList.append('wanaltic')
                        legendList.append(scheme_name)
                else:
                    plt.plot(r, wSolutions[i, n, :], lstyle[i + 1])
                    if n==0:
                        legendList.append(scheme_name)
                
        plt.title(titlestring)
        plt.legend(legendList)
        #plt.title('Flow in pipe')
        plt.xlabel('r-coordinate [-]')
        plt.ylabel('Velocity  [-]')

        
    def test_MES_convergence():
        """Code verification using the Method of Exact solution and observed order of accuracy:
        the RMS-errors are caluculated for space and time seperately 
        (when calculating space accuracy, the time is hold constant/unchanged and vice verca).
        The grid is then refined Ntds number of times with ratio 2 for space and 4 for time
        (recall that for this problem we use dt = D*h**2). The order-approximation is then calculated by comparing
        RMS-errors for the different gridrefinements using order-approx = log10(RMS_error(n)/RMS_error(n+1))/log10(refinementRatio)"""
        from math import sqrt
        from numpy import log2, log10
        startTime = timeModule.time()
        #discretize:
        N = 4 # No. of parts
        h = 1./N # length of r-step
        D = 0.4 # numerical diffusion number
        T = 0.1
        dt = D*h**2
        time = np.arange(0, T + dt , dt)
        print time
        Ntds = 5 #update h, and dt Ntds times

        r = np.linspace(0, 1, N +1)
        wstart = 1 - r**2 # initial profile is the pouseille profil
        
        schemes_theta = [0, 0.5, 1]
        schemes_name = ['FTCS', 'Crank', 'Laasonen']
        
        print """\n######## Started test: 'test_MES_convergence()': #########
         
        Solving startup of a flow in a pipe (Szymanski's problem) using theta schemes:
        {0}, {1}, {2}""".format(schemes_name[0], schemes_name[1], schemes_name[2])
        print """        Number of gridrefinements: {0}, dr-refinement ratio: {1} dt-refinement ratio: {2}""".format(Ntds, 2, 4)
        print """        Acceptence criteria: Consistency and Observed order of accuracy\n"""
        #empty lists containing RMS-errors, and order-approximation for all schemes
        schemes_rms_error_time = [[], [], []]
        schemes_rms_error_distance = [[], [], []]
        schemes_rms_error = [[], [], []]
        schemes_order = [[], [], []]
        schemes_order_time = [[], [], []]
        schemes_order_space = [[], [], []]
        schemes_max_error = [[], [], []]
        dtList = []
        hList = []
        
        #iterate all theta-schemes, solve them with solveNextTimestep tdma solver, and store solutions
        eps = 2.2204*10**-16 #tollerance to be used in calculations of analytical
        for n in range(Ntds):
            dtList.append(dt)
            hList.append(h)
            jumpTime = 4**n
            jumpSpace = 2**n
            itTimestart = timeModule.time()
            for scheme, theta in enumerate(schemes_theta):
                wold = wstart
                errorTime = 0
                temporalElements = 0
                error = 0
                errorElements = 0
                maxerror = 0
                
                for i, t in enumerate(time[1:]):             
                    wAnalytic = np.zeros_like(r) # the entire analytical solution is only needed at t = T
                    
                        
                    for k, r_variable in enumerate(r):
                        wAnalytic[k] = analyticSolution(r_variable, t, eps)
                            
                    wnew = solveNextTimestepStartup(theta, D, N, wold)
                    wold = wnew
                    
                    error += np.sum((wnew - wAnalytic)**2)
                    errorElements += len(wnew)
                    temperror = max(abs((wnew - wAnalytic)))
                    if temperror > maxerror:
                        maxerror = temperror
                    if float(i + 1)/jumpTime == temporalElements + 1.0:
                        #if temporalElements<3:
                        if scheme == 0:
                            print "sampling at r = {0}, t = {1}".format(r[N/2], t)

                        errorTime += (wnew[N/2] - analyticSolution(r[N/2], t, eps))**2 # add error for each timestep 
                        temporalElements += 1
                errorTime = sqrt(errorTime/(temporalElements)) #rms error time
                error = sqrt(error/(errorElements))
                errorDistance = sqrt(np.sum((wnew[0::jumpSpace]-wAnalytic[0::jumpSpace])**2)/(len(r[0::jumpSpace]))) #rms error space
                schemes_rms_error_time[scheme].append(errorTime)
                schemes_rms_error[scheme].append(error)
                schemes_rms_error_distance[scheme].append(errorDistance)
                schemes_max_error[scheme].append(maxerror)
                
                if n>0:
                    # Calculate order approximation
                    schemes_order_time[scheme].append(log10(schemes_rms_error_time[scheme][n-1]/schemes_rms_error_time[scheme][n])/log10(4))
                    schemes_order[scheme].append(log10(schemes_rms_error[scheme][n-1]/schemes_rms_error[scheme][n])/log10(2))
                    schemes_order_space[scheme].append(log10(schemes_rms_error_distance[scheme][n-1]/schemes_rms_error_distance[scheme][n])/log10(2))
            
            itTimeend = timeModule.time()
            itTime = itTimeend - itTimestart    
            print """        iteration {0} of {1} finsihed :dr/h = {2}, dt = {3}, tsample = {4}, CPU-time = {5}""".format(n + 1, Ntds, h, dt, t, itTime)
            #Update numerical values and discretization
            N *= 2 # 
            h = 1./N 
            dt = D*h**2
            time = np.arange(0, T + dt, dt)
            r = np.linspace(0, 1, N + 1)
            wstart = 1 - r**2

        endTime = timeModule.time()
        
        print """\n        Finished convergencetest CPU time: {0}""".format(endTime - startTime)
        #order_approx_time = log10(schemes_rms_error_time[:-1]/schemes_rms_error_time[1:])/log10(4)
        
        
        plotErrorAndOrder(schemes_name, schemes_rms_error_distance, schemes_rms_error,
                      schemes_order_space, schemes_order, Ntds)
        
        #from sympy_Newton_solver import optimize_error_cxct, optimize_error
        hxList = np.asarray(hList)
        htList = np.asarray(dtList)
        #legendList = []
        #lstyle = ['b', 'r', 'g', 'm']
        #fig , axarr = plt.subplots(2, 4, squeeze=False)
        #fig2 , axarr2 = plt.subplots(2, 3, squeeze=False)
        
        from sympy import symbols, lambdify, latex
        h_x, h_t = symbols('h_x h_t')
        errorfuncs = {}
        legends = []
        
        Ntds = len(schemes_rms_error[0])
        
        for k, scheme in enumerate(schemes_name):
            
            with open(scheme + '_error.txt', 'w') as filename:
                for n in range(len(hxList)):
                    
                    
                    filename.write(str(hxList[n]))
                    filename.write('    ')
                    filename.write(str(htList[n]))
                    filename.write('    ')
                    filename.write(str(schemes_rms_error[k][n]))
                    filename.write('    ')
                    filename.write(str(schemes_rms_error_distance[k][n]))
                    filename.write('    ')
                    filename.write(str(schemes_rms_error_time[k][n]))
                    filename.write('    ')
                    filename.write(str(schemes_max_error[k][n]))
                    filename.write('\n')
                    
            filename.close()
#             if solver=='macCormack':
#                 theoretical = 2
#         
#             if solver=='lax_friedrich':
#                 theoretical = 1
#         
#             if solver=='ftbs':
#                 theoretical = 1
#             print "\n"
#             print solver
#             gx0, gt0 = optimize_error_cxct(schemes_rms_error[k], hxList, htList, solver)
#             print "gx0, gt0: ", gx0, gt0 
#             gx, p, gt, q = optimize_error(schemes_rms_error[k], hxList, htList, gx0, gt0, solver)
#             print "gx, p, gt, q: ", gx, p, gt, q 
#             errorexpr = gx*h_x**p + gt*h_t**q
#             
#             print errorexpr
#             errorX = lambdify([h_x], gx*h_x**p, np)
#             errort = lambdify([h_t], gt*h_t**q, np)
#             errorfuncs[solver] = errorexpr
#             tempfunc = lambdify([h_x, h_t], errorexpr, np)
#             errorcalc = tempfunc(hxList, htList)
        

    
    def test_Stability():
        """This is a test of the Stability of the theta schemes:
        Version1 of FTCS has a varying stability limit D, somewhere between 0.25 and 0.5, depending on dt and dr.
        The ustability is accosiated with the the discretization of the BC dw/dr(r=0) = 0. Version2 of FTCS has 
        a different discretization of this BC, and has a stability-limit of D = 0.5 for any r. Crank-Nicolson and Laasonen
        are unconditionally stable.
        The test can create an animation and also plots solutions for different schemes and D-values."""
        
        print """\n######## Started test: 'test_Stability()': #########
        """
        N = 10 # No. of parts
        h = 1./N # length of r-step
        D = 0.45 # numerical diffusion number
        dt = D*h**2
        nmax = 60 # No. of time-steps
        T = nmax*dt
        time = np.linspace(0, T, nmax + 1)
        r = np.linspace(0, 1, N +1)
        wstart = 1 - r**2 # initial profile is the pouseille profil
        eps = 2.2204*10**-16 #tollerance to be used in calculations of analytical solution
        solvers = [solveNextTimestepStartup, solveNextTimestepStartupV2]
        
        animate = 'n'#= raw_input("        Do you want to animate resuts ? (y/n) ")
        
        if animate=='y':
            
            solver = solvers[0]
            schemes_name = ['FTCS', 'Crank-Nicolson', 'Laasonen']
            schemes_theta = [0, 0.5, 1]
            schemes_global_rms_error = []
            #initialize solution matrices
            wSolutions = np.zeros((len(schemes_theta), len(time), N + 1))
            wAnalytic = np.zeros((len(time), N + 1))
    
            
            #iterate all theta-schemes, solve them with solveNextTimestep tdma solver, and store solutions
            
            for scheme, theta in enumerate(schemes_theta):
                wold = wstart
                for i, t in enumerate(time[1:]):
                    wnew = solver(theta, D, N, wold)
                    wold = wnew
                    wSolutions[scheme, i + 1, :] = 1 - r**2 - wnew
                    
                    if scheme==0:
                        for k, r_variable in enumerate(r):
                            wAnalytic[i + 1, k] = 1 - r_variable**2 - analyticSolution(r_variable, t, eps)
                
                Nelements = (len(time)-1)*(N)
                error = sqrt(np.sum(((wSolutions[scheme, 1:, :-1] - wAnalytic[1:,:-1])/wAnalytic[1:,:-1])**2)/Nelements)
                schemes_global_rms_error.append(error)
    
            createAnimation(wSolutions, wAnalytic, schemes_name, r, time, jump = int(round(0.002/dt)))
        
        startTime = timeModule.time()
        schemes_names = ['FTCSV1', 'FTCSV2', 'Crank-Nicolson', 'Laasonen']
        schemes_theta = [0, 0, 0.5, 1]
        schemes_solver = [solvers[0], solvers[1], solvers[0], solvers[0]]
        NumericaldiffusionList =[0.4, 0.45, 0.5, 0.53, 1]
        
        # change some default values to make plots more readable 
        LNWDT=1; FNT=10
        plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT
        f2, axarr = plt.subplots(len(NumericaldiffusionList), len(schemes_names), sharex=True, sharey=True)
        lstyle = ['b', 'r', 'g', 'c']
        
        for row, D in enumerate(NumericaldiffusionList):
            
            dt = D*h**2
            nmax = 150 # No. of time-steps
            T = nmax*dt
            time = np.linspace(0, T, nmax + 1)
            wAnalytic = np.zeros((len(time), N + 1))
            for scheme, theta in enumerate(schemes_theta):
                wold = wstart
                drawing = 1
                for i, t in enumerate(time[1:]):

                    if scheme==0:
                        wa = np.zeros_like(r)
                        for k, r_variable in enumerate(r):
                            wa[k] = 1 - r_variable**2 - analyticSolution(r_variable, t, eps)
                        wAnalytic[i] = wa
                    wnew = schemes_solver[scheme](theta, D, N, wold)
                    wold = wnew
                    wnew = 1 - r**2 - wnew
                    error = np.sum((wnew-wAnalytic[i])**2) # used as criteria for stability
                    if error < 0.1 and (round(i/5)-drawing)==0:
                        #only draw each fifth solution and if the solution is stable
                        drawing = drawing + 1
                        axarr[row][scheme].plot(wnew, r, lstyle[scheme])
                        axarr[row][scheme].set_ylim(0, 1.1)
                        axarr[row][scheme].set_xlim(0, 1)
                if row == 0:
                    
                    axarr[row][scheme].set_title(schemes_names[scheme])
                if row == len(NumericaldiffusionList)-1:
                    axarr[row][scheme].set_xlabel('Velocity [-]')
                if scheme == 0:
                    axarr[row][scheme].set_ylabel('radius [-]')
                    axarr[row][scheme].text(-0.75, 0.5, '$D={0}$'.format(D))
        #plt.suptitle('test_Stability: Results from Stability test with different schemes and Difussion number D')
#        plt.savefig('../fig-ch5/startup_stability_transparent1.png', transparent=True)
        endTime = timeModule.time()
        print """        Finished Stability test
        CPU time: {0}""".format(endTime - startTime)
        #plt.show()
        
        
    def test_CPUtime():
        """This is a test of the speed of the solvers: Startup and StartupV3. """
        print """\n######## Started test: 'test_CPUtime()': #########
        """
        N = 200 # No. of parts
        h = 1./N # length of t-step
        D = 0.4 # numerical diffusion number
        dt = D*h**2
        nmax = 2000 # No. of time-steps
        theta = 0.5
        r = np.linspace(0, 1, N +1)
        wold = 1 - r**2
        
        startTimeStartup = timeModule.time()
        solveAllTimestepsStartup(theta, D, nmax, N, wold)
        endTimeStartup = timeModule.time()
        CPUTimeStartup = endTimeStartup - startTimeStartup
        
        startTimeStartupV3 = timeModule.time()
        solveAllTimestepsStartupV3(theta, D, nmax, N, wold)
        endTimeStartupV3 = timeModule.time()
        CPUTimeStartupV3 = endTimeStartupV3 - startTimeStartupV3
        
        print """        CPU-times:
        Startup: {0}
        StartupV2: {1}
        Ratio: {2}
        """.format(CPUTimeStartup, CPUTimeStartupV3, CPUTimeStartup/CPUTimeStartupV3)
        
        #plt.show()

    #test_MES()
    #test_MES_convergence()
    test_Stability()
    #test_CPUtime()
    plt.show()
                      
        

    
    

