# src-ch5/startup.py

import numpy as np
import matplotlib.pylab as plt
from matplotlib import animation
import time as timeModule

LNWDT=3; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT

def animateSphere(r,Tmatrix, ntot):
        
    fig = plt.figure()
    ax = plt.axes(xlim=(0,5),ylim = (0,310))
    line1, = ax.plot([],[], lw=2)
    #time_text = ax.text(1,100,'',transform=ax.transAxes)
    plt.xlabel('r [cm]')
    plt.ylabel('T [$^o C$')
        
        
    def init():
        line1.set_data([],[])
        #time_text.set_text('')
        return line1,# time_text
        
    def animate(i):
            
        Rnow = Tmatrix[:,i*jump]
        line1.set_data(r,Rnow)
        #time_text = time_text.set_text(str(i*dtau))
        return line1,# time_text
        
    jump = 5        
    Npictures = ntot/jump
        
    anim = animation.FuncAnimation(fig, animate,init_func=init,frames=Npictures,interval=1,blit=False)
    plt.show()

def createAnimation(NumericalSolutions, analyticalSolution, scheme_names, r_vector, time, jump=1, symmetric=True):

    """Method that creates and saves the animation of uNumeric and uAnalytic. 
    """
    #print """\n        Creating animation"""
     
    #startTime = timeModule.time()
    fig = plt.gcf()
    if symmetric:
        xlim=(np.min(analyticalSolution), np.max(analyticalSolution))
        ylim=(-r_vector[-1],r_vector[-1])
    else:
        xlim=(np.min(analyticalSolution), np.max(analyticalSolution))
        ylim=(r_vector[0],r_vector[-1])
    ax = plt.axes(xlim=xlim, ylim=ylim)
    time_text = ax.text(0.5, 0.95, '', transform=ax.transAxes)
    dt = time[1]-time[0]
    
    lines=[]     # list for plot lines for solvers and analytical solutions
    legends=[]   # list for legends for solvers and analytical solutions
    lstyle = ['r-', 'b--', 'c.', 'y:']
    line, = ax.plot([], [], lstyle[0]) #add extra plot line for analytical solution
    lines.append(line)
    legends.append('Analytical')
    
    lstyleit = 1
    for scheme in scheme_names:
        line, = ax.plot([], [],lstyle[lstyleit])
        lines.append(line)
        legends.append(scheme)
        lstyleit += 1
    
    #plt.title(titlestring)
    plt.xlabel('Velocity [-]')
    plt.ylabel('r-coordinate [-]')
    #plt.title('test_MES: Animation from test_MES')
    plt.legend(legends, loc='best', frameon=False)
     
    # initialization function: plot the background of each frame
    def init():
        #line.set_data([], [])
        time_text.set_text('')
        for line in lines:
            line.set_data([], [])
        return lines, time_text
    
    # animation function.  This is called sequentially
    
    if symmetric:
        r_vector2 = -1*r_vector[::-1]
        r_vector = np.append(r_vector2, r_vector)
    
    def animate_alt(i):
        i = i*jump
        time = i*dt
        time_text.set_text('time = %.3f' % time)
        for k, line in enumerate(lines):
            if (k==0):
                v_vector = analyticalSolution[i,:]
                if symmetric:
                    v_vector2 = v_vector[::-1]
                
                    v_vector = np.append(v_vector, v_vector2)
                line.set_data(v_vector, r_vector)
            else:
                
                v_vector = NumericalSolutions[k-1,i,:]
                if symmetric:
                    v_vector2 = v_vector[::-1]
                    v_vector = np.append(v_vector, v_vector2)
                line.set_data(v_vector, r_vector)

        return lines, time_text

    if jump<1:
        jump = 1
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate_alt, init_func=init, frames=len(time)/jump, interval=10, blit=False)
#    Writer = animation.writers['ffmpeg']
#    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
#    anim.save('../mov-ch5/couette_FTCS.mov', writer=writer)
    #endTime = timeModule.time()
    #print """        Animation created, CPU time: {0}""".format(endTime - startTime)
    
    plt.show()
    
def plotErrorAndOrder(schemesName, spaceErrorList,temporalErrorList,
                      spaceOrderList, temporalOrderList, Ntds):
    legendList = []
    lstyle = ['b', 'r', 'g', 'm']
    fig , axarr = plt.subplots(2, 2, squeeze=False)
    for k, scheme_name in enumerate(schemesName):
        axarr[0][0].plot(np.log10(np.asarray(spaceErrorList[k])),lstyle[k])
        axarr[0][1].plot(np.log10(np.asarray(temporalErrorList[k])),lstyle[k])
        axarr[1][0].plot(spaceOrderList[k],lstyle[k])
        axarr[1][1].plot(temporalOrderList[k],lstyle[k])
        legendList.append(scheme_name)
    plt.suptitle('test_MES_convergence(): Results from convergence test using Method of Exact Solution')

    axarr[1][0].axhline(1.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
    axarr[1][0].axhline(2.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
    axarr[1][1].axhline(1.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
    axarr[1][1].axhline(2.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
    axarr[1][0].set_ylim(0, 5)
    axarr[1][1].set_ylim(0, 5)
    axarr[0][0].set_ylabel('rms Error')
    axarr[0][0].set_title('space Error')
    axarr[1][0].set_ylabel('rms Error')
    axarr[0][1].set_title('temporal Error')
    axarr[1][0].set_ylabel('order')
    axarr[1][0].set_title('space order')
    axarr[1][1].set_ylabel('order')
    axarr[1][1].set_title('temporal order')
    axarr[0][1].legend(legendList, frameon=False)
    

    



