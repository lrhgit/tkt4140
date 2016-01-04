# ../Kap6/advection_schemes.py

import numpy as np
import os, glob
import matplotlib.pyplot as plt
from scipy import interpolate
from numpy import where
from math import sin

LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT




def init_step(x):
    """Assigning a value of 1.0 for values less than 0.1"""
    f = np.zeros_like(x)
    f[np.where(x <= 0.1)] = 1.0
    return f

def init_sine(x, T):
    """A smooth sin^2 function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = 0
    x_right = T/2
    xm = (x_right-x_left)/2.0
    f = where((x>x_left) & (x<x_right), np.sin(np.pi*(x-x_left)/(x_right-x_left))**2,f) 
    return f

def ftbs(u): # forward time backward space
    u[1:-1] = (1-c)*u[1:-1] + c*u[:-2]
    return u[1:-1]

# Lax-Wendroff
def lax_wendroff(u): 
    u[1:-1] = c/2.0*(1+c)*u[:-2] + (1-c**2)*u[1:-1] - c/2.0*(1-c)*u[2:]
    return u[1:-1]

# Lax-Friedrich Flux formulation
def lax_friedrich_Flux(u):
    u[1:-1] = (u[:-2] +u[2:])/2.0 -  dt*(F(u[2:])-F(u[:-2]))/(2.0*dx)
    return u[1:-1] 

# Lax-Friedrich Advection
def lax_friedrich(u):
    u[1:-1] = (u[:-2] +u[2:])/2.0 -  c*(u[2:] - u[:-2])/2.0
    return u[1:-1] 

# macCormack for advection quation
def macCormack(u):
    up = u.copy()
    up[:-1] = u[:-1] - c*(u[1:]-u[:-1])
    u[1:] = .5*(u[1:]+up[1:] -  c*(up[1:]-up[:-1]))
    return u[1:-1] 

a = 1.0 # wave speed
tmin, tmax = 0.0, 1.0 # start and stop time of simulation
xmin, xmax = 0.0, 1.5 # start and end of spatial domain
Nx = 100 # number of spatial points
CFL_list = np.linspace(0.1,1,10) #0.9 # courant number, need c<=1 for stability

x = np.linspace(xmin, xmax, Nx+1) # discretization of space
dx = float((xmax-xmin)/Nx) # spatial step size
f = init_step
sampletimes = np.linspace(0,1,11)
sampletimes[-1] = 1
solvers = [ftbs, macCormack, lax_friedrich]
if not os.path.isdir('solutions'):
    os.mkdir('solutions')
    print "######## Solutions directory does not exist. Creating one. ########"
else:
    
    path = 'solutions/step/'
    print "######## Solutions directory for step exists. Deleting it########"
    for solver in os.listdir(path):
        for cflsample in os.listdir(path + solver):
            for filename in os.listdir(path + solver + "/" + cflsample):
                if ".txt" in filename:
                    os.remove(path + solver + "/" + cflsample+ "/" + filename)
            os.removedirs(path + solver + "/" + cflsample)
    print "####### Deleted step" 
    path = 'solutions/sine/'
    print "######## Solutions directory for sine exists. Deleting it########"
    for solver in os.listdir(path):
        for Periodsample in os.listdir(path + solver):
            for cflsample in os.listdir(path + solver + "/" + Periodsample):
                for filename in os.listdir(path + solver + "/" + Periodsample + "/" + cflsample):
                    if ".txt" in filename:
                        os.remove(path + solver + "/" + Periodsample + "/" + cflsample+ "/" + filename)
        
    
                os.removedirs(path + solver + "/" + Periodsample + "/" + cflsample)
    print "##### Deleted solutions directory. Creating a new one: "
    os.mkdir('solutions')

tol = 0.0000002 
os.mkdir('solutions/step/') 
for k, solver in enumerate(solvers): # Solve for all solvers in list
    
    os.mkdir('solutions/step/' + solver.func_name)
    #path = os.path.join(solver.func_name)
    
    
    for c in CFL_list:
        u = f(x)
        # Discretize
        dt = c/a*dx # stable time step calculated from stability requirement
        if str(c) == '0.6' or str(c) == '0.8':
            tol = 0.0201
        else:
            tol = 0.0201  
        time = np.arange(tmin, tmax + dt, dt) # discretization of time
        os.mkdir('solutions/step/' + solver.func_name + '/' + str(c)) 
        cpath = 'solutions/step/' + solver.func_name + '/' + str(c) + '/'
        #with open(path + '0')
        sample = 0
        tpath = cpath + str(sampletimes[sample]) + '.txt'
        #os.mkdir(cpath + str(sampletimes[sample]))
 
        with open(tpath,'w') as filename:
            for indice, x_variable in enumerate(x):
                filename.write(str(x_variable))
                filename.write(' ')
                filename.write(str(u[indice]))
                filename.write('\n')
        filename.close()
        sample +=1
        
        for i, t in enumerate(time[1:]):
      
                  
            u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
              
            u[1:-1] = solver(u[:]) # calculate numerical solution of interior
            u[-1] = u_bc(x[-1] - a*dt) # interpolate along a characteristic to find the boundary value
            #u[0] = f(0)
            
            if abs(t - sampletimes[sample]) < tol:
                tpath = cpath + str(sampletimes[sample]) + '.txt'
                with open(tpath,'w') as filename:
                    for indice, x_variable in enumerate(x):
                        filename.write(str(x_variable))
                        filename.write(' ')
                        filename.write(str(u[indice]))
                        filename.write('\n')
                    filename.close()
                sample +=1
            if abs(sample - len(sampletimes)) < tol:
                sample = 0

f = init_sine
os.mkdir('solutions/sine/') 
Periods = np.linspace(0.2, 1, 5)
for k, solver in enumerate(solvers): # Solve for all solvers in list
    
    os.mkdir('solutions/sine/' + solver.func_name)
    #path = os.path.join(solver.func_name)
    
    for T in Periods:
        
        os.mkdir('solutions/sine/' + solver.func_name + '/' + str(T))
        for c in CFL_list:
            u = f(x, T)
            # Discretize
            dt = c/a*dx # stable time step calculated from stability requirement
            if str(c) == '0.6' or str(c) == '0.8':
                tol = 0.0201
            else:
                tol = 0.0201  
            time = np.arange(tmin, tmax + dt, dt) # discretization of time
            os.mkdir('solutions/sine/' + solver.func_name + '/' + str(T) + '/' + str(c)) 
            cpath = 'solutions/sine/'  + solver.func_name + '/' + str(T) + '/' + str(c) + '/'
            #with open(path + '0')
            sample = 0
            tpath = cpath + str(sampletimes[sample]) + '.txt'
            #os.mkdir(cpath + str(sampletimes[sample]))
     
            with open(tpath,'w') as filename:
                for indice, x_variable in enumerate(x):
                    filename.write(str(x_variable))
                    filename.write(' ')
                    filename.write(str(u[indice]))
                    filename.write('\n')
            filename.close()
            sample +=1
            
        
            for i, t in enumerate(time[1:]):
          
                      
                u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
                  
                u[1:-1] = solver(u[:]) # calculate numerical solution of interior
                u[-1] = u_bc(x[-1] - a*dt) # interpolate along a characteristic to find the boundary value
                #u[0] = f(0)
                
                if abs(t - sampletimes[sample]) < tol:
                    tpath = cpath + str(sampletimes[sample]) + '.txt'
                    with open(tpath,'w') as filename:
                        for indice, x_variable in enumerate(x):
                            filename.write(str(x_variable))
                            filename.write(' ')
                            filename.write(str(u[indice]))
                            filename.write('\n')
                        filename.close()
                    sample +=1
                if abs(sample - len(sampletimes)) < tol:
                    sample = 0

