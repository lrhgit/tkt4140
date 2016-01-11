# src/src-ch5/sphereAnimation.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation



def animateSphere(r,Tmatrix, Tmatrix_analytic, ntot, dtau):
        
    fig = plt.figure()
    ax = plt.axes(xlim=(0,6),ylim = (0,350))
    line1, = ax.plot([],[], lw=2, color='b')
    line2, = ax.plot([],[], lw=2, linestyle='--', color='r')
    line3, = ax.plot([],[], lw=2, linestyle=':', color='k')
    line4, = ax.plot([],[], lw=2, color='k')
    
    #time_text = ax.text(1,100,'',transform=ax.transAxes)
    plt.xlabel('r [cm]')
    plt.ylabel('T [$^o C$')
    plt.legend(['sphere-Analytic', 'sphere-FTCS', 'water'], frameon= False, loc='best')
    #plt.xticks([0, 1, 2, 3, 4, 5])
    
    xh = np.linspace(5, 6, 11)
    yh = np.ones_like(xh)*20
    
    xv = np.array([5, 5])
    yv = np.array([0, 350])
    def init():
        
        line1.set_data([],[])
        line2.set_data([],[])
        line3.set_data([],[])
        line4.set_data([],[])
        #time_text.set_text('')
        return line1
        
    def animate(i):
            
        Rnow = Tmatrix[i*jump, :]
        Rnow_analytic = Tmatrix_analytic[i*jump,:]
        
        line1.set_data(r,Rnow_analytic)
        line2.set_data(r,Rnow)
        line3.set_data(xh, yh)
        line4.set_data(xv, yv)
        
        plt.title('t = {0} s'.format(int(i*jump*dtau)))
        plt.xticks([0, 1, 2, 3, 4, 5])
        plt.yticks([20, 50, 100, 150, 200, 250, 300])
        #time_text = time_text.set_text(str(i*dtau))
        return line1, line2
        
    jump = 5        
    Npictures = ntot/jump
        
    anim = animation.FuncAnimation(fig, animate,init_func=init,frames=Npictures,interval=1,blit=False)
#    Writer = animation.writers['ffmpeg']
#    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
#    anim.save('../fig-ch5/sphere.mov', writer=writer)
    plt.show()



