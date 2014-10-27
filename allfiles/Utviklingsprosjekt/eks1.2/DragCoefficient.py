'''
Created on 23. okt. 2014

@author: leifh
'''
from numpy import linspace,array,append,logspace,zeros_like,where,vectorize,\
    logical_and
import numpy as np  
from matplotlib.pyplot import loglog,xlabel,ylabel,grid,savefig,show,rc,hold
from test.test_heapq import func_names

def cd_sphere(Re):
    "Computes the drag coefficient of a sphere as a function of the Reynolds number Re."
    # Curve fitted after fig . A -56 in Evett & Liu :% " Fluid Mechanics & Hydraulics ",
    # Schaum ' s Solved Problems McGraw - Hill 1989.
    
    from numpy import log10,array,polyval    
    
    if Re <= 0.0:
        CD = 0.0
    elif Re > 8.0e6:
        CD = 0.2
    elif Re > 0.0 and Re <= 0.5:
        CD = 24.0/Re
    elif Re > 0.5 and Re <= 100.0:
        p = array([4.22,-14.05,34.87,0.658])
        CD = polyval(p,1.0/Re) 
    elif Re > 100.0 and Re <= 1.0e4:
        p = array([-30.41,43.72,-17.08,2.41])
        CD = polyval(p,1.0/log10(Re))
    elif Re > 1.0e4 and Re <= 3.35e5:
        p = array([-0.1584,2.031,-8.472,11.932])
        CD = polyval(p,log10(Re))
    elif Re > 3.35e5 and Re <= 5.0e5:
        x1 = log10(Re/4.5e5)
        CD = 91.08*x1**4 + 0.0764
    else:
        p = array([-0.06338,1.1905,-7.332,14.93])
        CD = polyval(p,log10(Re))
    return CD

def cd_sphere_vector(Re):
    "Computes the drag coefficient of a sphere as a function of the Reynolds number Re."
    # Curve fitted after fig . A -56 in Evett & Liu :% " Fluid Mechanics & Hydraulics ",
    # Schaum ' s Solved Problems McGraw - Hill 1989.

    from numpy import log10,array,polyval
    CD = zeros_like(Re)
   
    CD = where(Re<0,0.0,0.0)     # condition 1
    
    CD = where((Re > 0.0) & (Re <=0.5),24/Re,CD) # condition 2

    p = array([4.22,-14.05,34.87,0.658])
    CD = where((Re > 0.5) & (Re <=100.0),polyval(p,1.0/Re),CD) #condition 3

    p = array([-30.41,43.72,-17.08,2.41])
    CD = where((Re >100.0)  & (Re <=1.0e4) ,polyval(p,1.0/log10(Re)),CD) #condition 4

    p = array([-0.1584,2.031,-8.472,11.932])
    CD = where((Re > 1.0e4)  &  (Re <=3.35e5),polyval(p,log10(Re)),CD) #condition 5

    CD = where((Re > 3.35e5) & (Re <=5.0e5),91.08*(log10(Re/4.5e5))**4 + 0.0764,CD) #condition 6

    p  = array([-0.06338,1.1905,-7.332,14.93])
    CD = where((Re > 5.05e5)  &  (Re <=8.0e6),polyval(p,log10(Re)),CD) #condition 7
    
    CD = where(Re>8.0e6,0.2,CD)  # condition 8

    return CD

def cv_sphere_vector(Re):
    "Computes the drag coefficient of a sphere as a function of the Reynolds number Re."
    # Curve fitted after fig . A -56 in Evett & Liu :% " Fluid Mechanics & Hydraulics ",
    # Schaum ' s Solved Problems McGraw - Hill 1989.

    from numpy import log10,array,polyval
    CD = zeros_like(Re)
    CD = where(Re<0,0.0,0.0)
    CD = where(Re>8.0e6,0.2,CD)

    CD = where((Re > 0.0) & (Re <=0.5),24/Re,CD)

    p = array([4.22,-14.05,34.87,0.658])
    CD = where((Re > 0.5) & (Re <=100.0),polyval(p,1.0/Re),CD)

    p = array([-30.41,43.72,-17.08,2.41])
    CD = where((Re >100.0)  & (Re <=1.0e4) ,polyval(p,1.0/log10(Re)),CD)

    p = array([-0.1584,2.031,-8.472,11.932])
    CD = where((Re > 1.0e4)  &  (Re <=3.35e5),polyval(p,log10(Re)),CD)

    CD = where((Re > 3.35e5) & (Re <=5.0e5),91.08*(log10(Re/4.5e5))**4 + 0.0764,CD)

    p  = array([-0.06338,1.1905,-7.332,14.93])
    CD = where((Re > 5.05e5)  &  (Re <=8.0e6),polyval(p,log10(Re)),CD)

    return CD

def cd_sphere_vector_boolean(Re):
    from numpy import log10,array,polyval
       
    condition1 = Re < 0
    condition2 = logical_and(0 < Re, Re <= 0.5)
    condition3 = logical_and(0.5 < Re, Re <= 100.0)
    condition4 = logical_and(100.0 < Re, Re <= 1.0e4)
    condition5 = logical_and(1.0e4 < Re, Re <= 3.35e5)
    condition6 = logical_and(3.35e5< Re, Re <= 5.0e5)
    condition7 = logical_and(5.0e5 < Re, Re <= 8.0e6)
    condition8 = Re > 8.0e6
    
    cd = zeros_like(Re)
    cd[condition1] = 0.0
    
    cd[condition2] = 24/Re[condition2]
    
    p = array([4.22,-14.05,34.87,0.658])
    cd[condition3] = polyval(p,1.0/Re[condition3]) 
    
    p = array([-30.41,43.72,-17.08,2.41])
    cd[condition4] = polyval(p,1.0/log10(Re[condition4]))
    
    p = array([-0.1584,2.031,-8.472,11.932])
    cd[condition5] = polyval(p,log10(Re[condition5]))
    
    cd[condition6] = 91.08*(log10(Re[condition6]/4.5e5))**4 + 0.0764
    
    p  = array([-0.06338,1.1905,-7.332,14.93])
    cd[condition7] = polyval(p,log10(Re[condition7]))
    
    cd[condition8] = 0.2
    
    return cd
    

if __name__ == '__main__':              #Check whether this file is executed (name==maine) or imported as module
    
    ReNrs = logspace(-2,7,num=50000)
    CD    = zeros_like(ReNrs)

    import time
    
    t0 = time.clock()
    counter = 0
    for Re in ReNrs:
        CD[counter]=cd_sphere(Re)
        counter += 1
    t1 = time.clock()
    
    dt1 = t1-t0
    
    
    cd_sphere_auto_vector=vectorize(cd_sphere) # make a vectorized version of the function automatically
    t0 = time.clock()           
    CD2 = cd_sphere_auto_vector(ReNrs) # compute CD for all ReNrs
    t1 = time.clock()
    dt2 =t1-t0
    
    t0 = time.clock()               
    CD3 = cd_sphere_vector(ReNrs) # compute CD with our vectorized function
    t1 = time.clock()
    dt3 =t1-t0

    t0 = time.clock()               
    CD4 = cd_sphere_vector_boolean(ReNrs) # compute CD with our vectorized function
    t1 = time.clock()
    dt4=t1-t0
              
    print 'elapsed time', dt1, dt2, dt3, dt4
    
    #funcs =[cd_sphere,cd_sphere_vector,cd_sphere_auto_vector]     
    funcs =[cd_sphere_vector,cd_sphere_vector_boolean]     
    
    timings = {}
    for func in funcs:
        name = func.func_name   
        t0 = time.clock()
        r = func(ReNrs) 
        t1 = time.clock()
        timings[name] = t1 - t0
       
    print 'elapbsed time', dt1, dt2, dt3
            
    # set fontsize prms 
    fnSz = 16; font = {'size'   : fnSz}; rc('font',**font)          
    
    # plot the function    
    loglog(ReNrs,CD)
    hold('on')
    loglog(ReNrs,CD2)
    loglog(ReNrs,CD3)
    loglog(ReNrs,CD4)
    
    xlabel('$Re$',fontsize=fnSz)
    ylabel('$C_D$',fontsize=fnSz) 
    grid('on','both','both')
 #   savefig('example_sphere.png')
    show()