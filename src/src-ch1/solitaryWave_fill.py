# src-ch1/solitaryWave.py
#import matplotlib; matplotlib.use('Qt4Agg')
import matplotlib.pylab as plt
#plt.get_current_fig_manager().window.raise_()
import numpy as np

#### set default plot values: ####
LNWDT=3; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT

""" This script solves the problem with the solitary wave:

        y'' = a*3*y*(1-y*3/2)
        
        y(0) = 1, y'(0) = 0
        
    or as a system of first order differential equations (y0 = y, y1 = y'):
        
        y0' = y'
        y1' = a*3*y0*(1-y0*3/2)
        
        y0(0) = 1, y1(0) = 0
        
"""
a = 2./3
h = 0.2 # steplength dx
x_0, x_end = 0, 0.6

x = np.arange(x_0, x_end + h, h) # allocate x values

#### solution vectors: ####
Y0_euler = np.zeros_like(x) # array to store y values
Y1_euler = np.zeros_like(x) # array to store y' values

Y0_heun = np.zeros_like(x)
Y1_heun = np.zeros_like(x)

#### initial conditions: ####
Y0_euler[0] = 1 # y(0) = 1
Y1_euler[0] = 0 # y'(0) = 0

Y0_heun[0] = 1 
Y1_heun[0] = 0 


#### solve with euler's method ####

for n in range(len(x) - 1):
    y0_n = Y0_euler[n] # y at this timestep
    y1_n = Y1_euler[n] # y' at this timestep
    
    "Fill in lines below"
    f0 = 
    f1 = 
    "Fill in lines above"
    
    Y0_euler[n + 1] = y0_n + h*f0
    Y1_euler[n + 1] = y1_n + h*f1

#### solve with heun's method: ####

for n in range(len(x) - 1):
    y0_n = Y0_heun[n] # y0 at this timestep (y_n)
    y1_n = Y1_heun[n] # y1 at this timestep (y'_n)
    
    "Fill in lines below"
    f0 = 
    f1 = 
    
    y0_p = 
    y1_p = 
    
    f0_p = 
    f1_p = 
    "Fill in lines above"
    
    Y0_heun[n + 1] = y0_n + 0.5*h*(f0 + f0_p)
    Y1_heun[n + 1] = y1_n + 0.5*h*(f1 + f1_p)
    

Y0_taylor = 1 - x**2/2 + x**4/6
Y1_taylor = -x + (2./3)*x**3

Y0_analytic = 1./(np.cosh(x/np.sqrt(2))**2)


#### Print and plot solutions: ####

print "a) euler's method: y({0})={1}, y'({2})={3}".format(x_end, round(Y0_euler[-1], 4), x_end, round(Y1_euler[-1], 4))
print "b) heun's method: y({0})={1}, y'({2})={3}".format(x_end, round(Y0_heun[-1], 4), x_end, round(Y1_heun[-1], 4))
print "c) Taylor series: y({0})={1}, y'({2})={3}".format(x_end, round(Y0_taylor[-1], 4), x_end, round(Y1_taylor[-1], 4))
print "d) Analytical solution: y({0})={1}".format(x_end, round(Y0_analytic[-1], 4))

plt.figure()
plt.plot(x, Y0_euler, 'r-o')
plt.plot(x, Y0_heun, 'b-^')
plt.plot(x, Y0_taylor, 'g-*')
plt.plot(x, Y0_analytic, 'k--')

eulerLegend = 'euler, y({0})={1}'.format(x_end, round(Y0_euler[-1], 4))
heunLegend = 'heun, y({0})={1}'.format(x_end, round(Y0_heun[-1], 4))
taylorLegend = 'taylor, y({0})={1}'.format(x_end, round(Y0_taylor[-1], 4))
analyticLegend = 'analytic, y({0})={1}'.format(x_end, round(Y0_analytic[-1], 4))

plt.legend([eulerLegend, heunLegend, taylorLegend, analyticLegend], loc='best', frameon=False)
plt.show()





