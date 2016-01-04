# chapter5/src-ch6/startup.py;TRIdiagonalSolvers.py @ git@lrhgit/tkt4140/allfiles/digital_compendium/chapter5/src-ch6/TRIdiagonalSolvers.py;


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
        if t< 0:
            print "t is less then zero"
        temp = exp(-t*(n*pi)**2)*sin(n*pi*y)/n
        sumValue += temp
    u = 1 - y - (2/pi)*sumValue
    return u


def solveNextTimestepCouette(theta, D, N, uold):
    """ Method that sets up the tridiagonal theta-scheme for the transient couettflow.
        At time t=t0 the plate starts moving at y=0
        The method solves only for the next time-step, using standard tdma solver.
        The Governing equation is:
        du/dt = d^2(u)/dx^2 , u = u(y,t), 0 < y < 1
        Boundary conditions: u(0, t) = 1, u(1, t) = 0
        Initial condition: u(t, 0) = 0 t<0,  u(t, 0) = 1 t>0
        Args:
            theta(float): parameter between 0 and 1 seperating different schemes (FTCS, Laasonen, Crank-Nicolson...)
            D(float): Numerical diffusion number
            N(int): number of equations
            uold(array): solution from previous iteration

    
        Returns:
            unew(array): solution at time t^n+1
    """
    u0 = uold[0]
    uold = uold[1:] 
    N = N - 1 
    tmp1 = D*(1. - theta)
    
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)

    for j in range(1,N):
        a[j] = -D*theta
        b[j] = (1. + 2*D*theta)
        c[j] = -D*theta
        d[j] = tmp1*(uold[j - 1] + uold[j + 1]) + (1. - 2*tmp1)*uold[j]
        
    b[0] = (1 + 2*D*theta)
    c[0] = -D*theta
    d[0] = tmp1*(u0 + uold[1]) + (1. - 2*tmp1)*uold[0] + D*theta*u0
    unew = tdma(a, b, c, d) # solve with tdma
    unew = np.append(np.append(1., unew), 0.) # add BCs
    
    return unew

def solveNextTimestepCouetteCrank(theta, D, N, uold):
    """ Method that sets up the tridiagonal theta-scheme for the transient couettflow.
        At time t=t0 the plate starts moving at y=0
        The method solves only for the next time-step, using standard tdma solver.
        The Governing equation is:
        du/dt = d^2(u)/dx^2 , u = u(y,t), 0 < y < 1
        Boundary conditions: u(0, t) = 1, u(1, t) = 0
        Initial condition: u(t, 0) = 0 t<0,  u(t, 0) = 1 t>0
        Args:
            theta(float): parameter between 0 and 1 seperating different schemes (FTCS, Laasonen, Crank-Nicolson...)
            D(float): Numerical diffusion number
            N(int): number of equations
            uold(array): solution from previous iteration

    
        Returns:
            unew(array): solution at time t^n+1
    """
    u0 = uold[0]
    uold = uold[1:] 
    N = N - 1
    
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)

    for j in range(1,N):
        a[j] = 1.
        b[j] = -2*(1. + 1./D)
        c[j] = 1.
        d[j] = -(uold[j - 1] + uold[j + 1]) + 2*(1. - 1/D)*uold[j]
        
    b[0] = -2*(1. + 1/D)
    c[0] = 1
    d[0] = -(u0 + uold[1]) + 2*(1. - 1/D)*uold[0] - u0
    unew = tdma(a, b, c, d) # solve with tdma
    unew = np.append(np.append(1., unew), 0.) # add BCs
    
    return unew

            
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import time as timeModule
    from TRIdiagonalSolvers import tdma, tripiv
    from Visualization import createAnimation, plotErrorAndOrder
    from math import sqrt
    
    # change some default values to make plots more readable 
    LNWDT=3; FNT=11
    plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT
    
    def test_MES():
        """Code verification using the Method of Exact solution"""
        print """\n######## Started test 'test_MES()': #########"""
        N = 40 # No. of parts
        h = 1./N # length of r-step
        D = 0.2 # numerical diffusion number
        dt = D*h**2
        T = 0.05
        time = np.arange(0, T + dt, dt)
        y = np.linspace(0, 1, N +1)
        
        schemesTheta = [0, 0.5, 1]
        schemesName = ['FTCS', 'Crank-Nicolson', 'Laasonen']
        uStart = np.zeros_like(y) # initial profile is the pouseille profil
        uStart[0] = 1
        #initialize solution matrices
        uSolutions = np.zeros((len(schemesTheta),len(time), N + 1))
        uSolutions[:,0 , 0] = 1
        uAnalytic = np.zeros((len(time), N + 1))
        uAnalytic[0, 0] = 1
        
        #iterate all theta-schemes, solve them with solveNextTimestep tdma solver, and store solutions
        eps = 2.2204*10**-16 #tollerance to be used in calculations of analytical solution
        for scheme, theta in enumerate(schemesTheta):
            uOld = uStart
            for i, t in enumerate(time[1:]):
                
                if scheme==0:
                    
                    # Calculate analytical solution at discrete times and points
                    for k, y_variable in enumerate(y):
                        uAnalytic[i + 1, k] = analyticSolution(y_variable, t)
                
                # calculate numerical solutions at time t        
                uNew = solveNextTimestepCouetteCrank(theta, D, N, uOld)
                uOld = uNew
                uSolutions[scheme, i + 1, :] =  uNew
                
        animate = raw_input("        do yo want to animate? y/n")
        if animate == 'y':
            createAnimation(uSolutions, uAnalytic, schemesName, y, time, jump=2, symmetric=False)
        
                
        plt.figure()
        lstyle = ['r--', 'b:', 'c-.', 'y-.']
        lstyleAnalytic = 'k--'
        nsnapshotlist=np.linspace(0, len(time)-1, 8)
        legendList = []
        for n in nsnapshotlist:
            for k, scheme in enumerate(schemesName):
                if k == 0:
                    plt.plot(y, uAnalytic[n, :], lstyleAnalytic)
                    if n ==0:
                        legendList.append('Analytic')
                plt.plot(y, uSolutions[k, n, :], lstyle[k])
                if n==0:
                    legendList.append(scheme)
        plt.legend(legendList)
        
    def test_MES_convergence():
        """Code verification using the Method of Exact solution and observed order of accuracy:
        the RMS-errors are caluculated for space and time seperately 
        (when calculating space accuracy, the time is hold constant/unchanged and vice verca).
        The grid is then refined Ntds number of times with ratio 2 for space and 4 for time
        (recall that for this problem we use dt = D*h**2). The order-approximation is then calculated by comparing
        RMS-errors for the different gridrefinements using order-approx = log10(RMS_error(n)/RMS_error(n+1))/log10(refinementRatio)"""
        from math import sqrt
        from numpy import log2, log10, log
        from scipy import interpolate
        
        startTime = timeModule.time()
        
        """Code verification using the Method of Exact solution"""
        print """\n######## Started test: 'test_MES_convergence()': #########
        """
        
        N = 20 # No. of parts
        h = 1./N # length of r-step
        D = 0.4 # numerical diffusion number
        dt = D*h**2
        T = 0.005
        time = np.arange(0, T + dt, dt)
        y = np.linspace(0, 1, N +1)
        
        schemesTheta = [0, 0.5, 1.]
        schemesName = ['FTCS','Crank-Nicolson', 'Laasonen']

        uStart = np.zeros_like(y) # initial profile is the pouseille profil
        uStart[0] = 1
        
        errorList = [] #empty list of space Error
        TemporalErrorList = []
        orderList = []
        orderList2 = []
        hxList = []
        htList = []
        for p in range(len(schemesTheta)): 
            # create empty lists for all schemes for order and spaceError
            errorList.append([])
            TemporalErrorList.append([])
            orderList.append([])
            orderList2.append([])
        Ntds = 5
        #iterate all theta-schemes, Ntds number of times and caluculate spatial and temporal error:

        for n in range(Ntds):
            uAnalytic = np.zeros_like(y)
            for k, theta in enumerate(schemesTheta):
                uOld = uStart
                error = 0
                Elements = 0
                temporalElements = 0
                temperror = 0
                for i, t in enumerate(time[1:]):
                            
                    uNew = solveNextTimestepCouette(theta, D, N, uOld)
                    uOld = uNew
                    
                    for j, y_variable in enumerate(y):
                            uAnalytic[j] = analyticSolution(y_variable, t)
                    error += np.sum((uNew[1:-1]-uAnalytic[1:-1])**2)
                    Elements += len(uNew[1:-1])
                    temperror += (uNew[N/2]- uAnalytic[N/2])**2
                    temporalElements += 1
                error = sqrt(error/Elements)
                temperror = sqrt(temperror/temporalElements)
                
                errorList[k].append(error)
                TemporalErrorList[k].append(temperror)
                if n > 0:
                    orderList[k].append(np.log(errorList[k][n-1]/errorList[k][n])/np.log(2))
                    orderList2[k].append(np.log(TemporalErrorList[k][n-1]/TemporalErrorList[k][n])/np.log(2))
            #Update numerical values and discretization
            print "\n        finished iteration {0} of {1}, h = {2}, dt = {3}, tsample = {4}".format(n+1, Ntds, h, dt, t)

            hxList.append(h)
            htList.append(dt)
            N *= 2 # 
            h = 1./N 
            dt = D*h**2
            time = np.arange(0, T + dt, dt)
            y = np.linspace(0, 1, N + 1)
            uStart = np.zeros_like(y)
            uStart[0] = 1.
        
            legendList = []
            lstyle = ['b', 'r', 'g', 'm']
            fig , axarr = plt.subplots(2, 2, squeeze=False)
            for k, scheme_name in enumerate(schemesName):
                axarr[0][0].plot(errorList[k],lstyle[k])
                axarr[1][0].plot(orderList[k],lstyle[k])
                axarr[0][1].plot(TemporalErrorList[k],lstyle[k])
                axarr[1][1].plot(orderList2[k],lstyle[k])
                legendList.append(scheme_name)
            plt.suptitle('test_MES_convergence(): Results from convergence test using Method of Exact Solution')
        
            axarr[1][0].axhline(1.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
            axarr[1][0].axhline(2.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
            axarr[1][0].set_ylim(0, 5)

            axarr[0][0].set_ylabel('rms Error')
            axarr[0][0].set_title('space Error')
            axarr[1][0].set_ylabel('rms Error')

            axarr[1][0].set_ylabel('order')
            axarr[1][0].set_title('space order')
            axarr[0][0].legend(legendList, frameon=False)
        
        
        endTime = timeModule.time()
        print """\n        Finished Convergence test. CPU-time = {0}""".format(endTime - startTime)
        
    
        

    def test_Stability():
        """This is a test of the Stability of the theta schemes:
        Version1 of FTCS has a varying stability limit D, somewhere between 0.25 and 0.5, depending on dt and dr.
        The ustability is accosiated with the the discretization of the BC dw/dr(r=0) = 0. Version2 of FTCS has 
        a different discretization of this BC, and has a stability-limit of D = 0.5 for any r. Crank-Nicolson and Laasonen
        are unconditionally stable.
        The test can create an animation and also plots solutions for different schemes and D-values."""
        
        print """\n######## Started test: 'test_Stability()': #########
        """
        N = 50 # No. of parts
        h = 1./N # length of r-step
        D = 0.45 # numerical diffusion number
        dt = D*h**2
        nmax = 259 # No. of time-steps
        T = nmax*dt
        time = np.linspace(0, T, nmax + 1)
        y = np.linspace(0, 1, N +1)
        ustart = np.zeros_like(y)
        ustart[0] = 1
        eps = 2.2204*10**-16 #tollerance to be used in calculations of analytical solution
        
        startTime = timeModule.time()
        schemesNames = ['FTCS', 'Crank-Nicolson', 'Laasonen']
        schemesTheta = [0, 0.5, 1]
        NumericaldiffusionList =[0.4, 0.45, 0.5, 0.503, 0.51]
        
        # change some default values to make plots more readable 
        LNWDT=0.8; FNT=11
        plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT
        f2, axarr = plt.subplots(len(NumericaldiffusionList), len(schemesNames), sharex=True, sharey=True)
        lstyle = ['b', 'r', 'g', 'c']
        
        for row, D in enumerate(NumericaldiffusionList):
            uAnalytic = np.zeros((len(time), N + 1))
            for scheme, theta in enumerate(schemesTheta):
                uold = ustart
                drawing = 0
                jump = 2
                for i, t in enumerate(time[1:]):
                    if scheme==0:
                        #only calculate analytical solution once
                        wa = np.zeros_like(y)
                        for k, y_variable in enumerate(y):
                            wa[k] = analyticSolution(y_variable, t)
                        uAnalytic[i] = wa
                    unew = solveNextTimestepCouette(theta, D, N, uold)
                    uold = unew
                    error = np.sum((unew-uAnalytic[i])**2) # used as criteria for stability
                    if error < 1.1 and (round(i/jump)-drawing)==0:
                        #only draw each "jump" solution and if the solution is stable
                        drawing = drawing + 1
                        axarr[row][scheme].plot(y, unew, lstyle[scheme])
                        axarr[row][scheme].set_xlim(0, 1)
                        axarr[row][scheme].set_ylim(0, 1.1)
                        if i == 12:
                            # as time develops draw fewer plots
                            drawing = 4
                            jump = 4
                        if i == 32:
                            drawing = 5
                            jump = 8
                        if i == 96:
                            drawing = 6
                            jump = 16
                if row == 0:
                    axarr[row][scheme].set_title(schemesNames[scheme])
                if row == len(NumericaldiffusionList)-1:
                    axarr[row][scheme].set_xlabel('radius [-]')
                if scheme == 0:
                    axarr[row][scheme].set_ylabel('Velocity [-]')
                    axarr[row][scheme].text(-0.5, 0.5, 'D = {0}'.format(D))
        plt.suptitle('test_Stability: Results from Stability test with different schemes and numerical Difussion number D')
        
        endTime = timeModule.time()
        print """        Finished Stability test
        CPU time: {0}""".format(endTime - startTime)
        #plt.show()
        
        
        #plt.show()

    test_MES()
    #test_MES_convergence()
    #test_Stability()
    plt.show()
                        
        

    
    

