'''
Created on 23. okt. 2014

@author: leifh
'''
from numpy import linspace, array, append, logspace, zeros_like
import numpy as np  
from matplotlib.pyplot import loglog, xlabel, ylabel, grid, savefig, show, rc

def CDsphere(Re):
    "Computes the drag coefficient of a sphere as a function of the Reynolds number Re."
    # Curve fitted after fig . A -56 in Evett & Liu :% " Fluid Mechanics & Hydraulics ",
    # Schaum ' s Solved Problems McGraw - Hill 1989.
    
    from numpy import log10, array, polyval    
    
    if Re <= 0.0:
        CD = 0.0
    elif Re > 8.0e6:
        CD = 0.2
    elif Re > 0.0 and Re <= 0.5:
        CD = 24.0/Re
    elif Re > 0.5 and Re <= 100.0:
        p = array([4.22, -14.05, 34.87, 0.658])
        CD = polyval(p, 1.0/Re) 
    elif Re > 100.0 and Re <= 1.0e4:
        p = array([-30.41, 43.72, -17.08, 2.41])
        CD = polyval(p, 1.0/log10(Re))
    elif Re > 1.0e4 and Re <= 3.35e5:
        p = array([-0.1584, 2.031, -8.472, 11.932])
        CD = polyval(p, log10(Re))
    elif Re > 3.35e5 and Re <= 5.0e5:
        x1 = log10(Re/4.5e5)
        CD = 91.08*x1**4 + 0.0764
    else:
        p = array([-0.06338, 1.1905, -7.332, 14.93])
        CD = polyval(p, log10(Re))
    return CD


if __name__ == '__main__':              #Check whether this file is executed (name==maine) or imported as module
    
    ReNrs = logspace(-2,7,num=500)
    CD    = zeros_like(ReNrs)

    counter = 0
    for Re in ReNrs:
        CD[counter]=CDsphere(Re)
        counter += 1

    # set fontsize prms 
    fnSz = 16; font = {'size'   : fnSz}; rc('font', **font)          
    
    # plot the function    
    loglog(ReNrs, CD)
    xlabel('$Re$',fontsize=fnSz)
    ylabel('$C_D$',fontsize=fnSz) 
    grid('on','both','both')
    savefig('example_sphere.png')
    show()
