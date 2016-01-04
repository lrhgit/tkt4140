
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def optimize_error_cxct(error, hx, ht, scheme):
    
    if scheme=='FTCS':
        p = 2
        q = 1
        x0 = np.array([1, - 3 ])
        
    elif scheme=='Crank':
        p = 1.981
        q = 2
        x0 = np.array([0.0692, -0.4 ]) 
    elif scheme=='Laasonen':
        p = 2
        q = 1
        x0 = np.array([1, - 3 ])
    
    D = ht[0]/hx[0]**2
    def func(x, gx, gt):
        dt = D*x**2
        return gx*x**p + gt*dt**q
    
    
    
    xdata = hx
    ydata = error
    cx, ct = curve_fit(func, xdata, ydata, x0)[0]

    return cx, ct

def optimize_error(error, hx, ht, gx0, gt0, scheme):
    
    
    if scheme=='FTCS':
        p0 = 2
        q0 = 1
    elif scheme=='Crank':
        p0 = 1.98480884701
        q0 = 2
    elif scheme=='Laasonen':
        p0 = 2
        q0 = 1
    
    D = ht[0]/hx[0]**2
    print "D: ", D
    print "hx[0], ht[0],  D*hx[0]**2: ", hx[0], ht[0],  D*hx[0]**2
    def func(x, gx, p, gt, q):
        dt = D*x**2
        return gx*x**p + gt*dt**q
    

    x0 = np.array([gx0, p0, gt0, q0])
    print scheme
    print "x0", x0
    xdata = hx
    ydata = error
    gx, p, gt, q = curve_fit(func, xdata, ydata, x0)[0]
    gx, p, gt, q = round(gx, 2), round(p, 3), round(gt, 2), round(q, 3)
    return gx, p, gt, q

def optimize_error_gxp(error, hx, scheme):
    
    
    if scheme=='FTCS':
        p0 = 2
        q0 = 1
    elif scheme=='Crank':
        p0 = 2
        gx0 = 0.5
    elif scheme=='Laasonen':
        p0 = 2
        q0 = 1
    
    
    def func(x, gx, p):
        
        return gx*x**p
    

    x0 = np.array([gx0, p0])
    print scheme
    print "x0", x0
    xdata = hx
    ydata = error
    gx, p = curve_fit(func, xdata, ydata, x0)[0]
    
    return gx, p


def optimize_error_gtq(error, hx, ht, scheme, gx0, p0 ):
    
    
    
    D = ht[0]/hx[0]**2
    def func(x, gx, p, gt):
        dx = np.sqrt(x/D)
        return gx*dx**p + gt*x**2
    

    x0 = np.array([gx0, p0, -gx0])
    print scheme
    print "x0", x0
    xdata = ht
    ydata = error
    gx, p, gt = curve_fit(func, xdata, ydata, x0)[0]
    gx, p, gt =  round(gx, 2), round(p, 3), round(gt, 2),
    return gx, p, gt

FTCS = {}
Crank ={}
Laasonen = {}
solvers = ['FTCS', 'Crank', 'Laasonen']
for k, solver in enumerate(solvers):
    hx = []
    ht = []
    globalError = []
    spaceError = []
    temporalError = []
    maxError = []
    
    with open(solver + '_errors.txt', 'r') as filename:
        for line in filename:
            #print line
            splitline = line.split()
            hx.append(float(splitline[0]))
            ht.append(float(splitline[1]))
            globalError.append(float(splitline[2]))
            spaceError.append(float(splitline[3]))
            temporalError.append(float(splitline[4]))
            maxError.append(float(splitline[5]))
    rx = hx[0]/hx[1]
    rt = ht[0]/ht[1]

    eval(solver)['globalError'] =  globalError  
    eval(solver)['globalOrder'] = np.log(np.asarray(globalError[:-1])/np.asarray(globalError[1:]))/np.log(rx)
    eval(solver)['spaceError'] =  spaceError
    eval(solver)['spaceOrder'] = np.log(np.asarray(spaceError[:-1])/np.asarray(spaceError[1:]))/np.log(rx)
    eval(solver)['temporalError'] =  temporalError
    eval(solver)['temporalOrder'] = np.log(np.asarray(temporalError[:-1])/np.asarray(temporalError[1:]))/np.log(rt)
    eval(solver)['maxError'] =  maxError
    eval(solver)['maxOrder'] = np.log(np.asarray(maxError[:-1])/np.asarray(maxError[1:]))/np.log(rx)

    
print hx
print ht
print rx
print rt
print FTCS['globalError']
print FTCS['globalOrder']
print Crank['globalError']
print Crank['globalOrder']
print Laasonen['globalError']
print Laasonen['globalOrder']



#plt.plot(curvefiterror, 'b--')
#x = Newton_solver_gxgt(error, hx, ht, [cx, ct], p_value=2, q_value=2)
#Newton_solver_All(macCormack['globalError'], hx, ht, [cx, 2, ct, 2])

print "\n"


#x = Newton_solver_gxgt(macCormack['globalError'], hx, ht, [1, 1], p_value=2, q_value=2)

    #gx, p, gt, q = optimize_error(Crank['globalError'], hxList, htList, gx0, gt0, 'Crank')
    #print 'gx, p, gt, q', gx, p, gt, q
legendList = []
lstyle = ['b', 'r', 'g', 'm']
fig , axarr = plt.subplots(2, 4, squeeze=False)

legends = []

Ntds= len(eval(solver)['spaceError'])

for k, solver in enumerate(solvers):

    
    axarr[0][0].plot(eval(solver)['spaceError'],lstyle[k])
    axarr[0][1].plot(eval(solver)['temporalError'],lstyle[k])
    axarr[0][2].plot(eval(solver)['maxError'],lstyle[k])
    axarr[0][3].plot(eval(solver)['globalError'],lstyle[k])
    axarr[1][0].plot(eval(solver)['spaceOrder'],lstyle[k])
    axarr[1][1].plot(eval(solver)['temporalOrder'],lstyle[k])
    axarr[1][2].plot(eval(solver)['maxOrder'],lstyle[k])
    axarr[1][3].plot(eval(solver)['globalOrder'],lstyle[k])
    legendList.append(solver)






plt.suptitle('Results from error analysis on advection equation')

axarr[1][0].axhline(1.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
axarr[1][0].axhline(2.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
axarr[1][1].axhline(1.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
axarr[1][1].axhline(2.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')

axarr[1][2].axhline(1.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
axarr[1][2].axhline(2.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
axarr[1][3].axhline(1.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
axarr[1][3].axhline(2.0, xmin=0, xmax=Ntds-2, linestyle=':', color='k')

axarr[1][0].set_ylim(0, 3)
axarr[1][1].set_ylim(0, 3)
axarr[1][2].set_ylim(0, 3)
axarr[1][3].set_ylim(0, 3)
axarr[0][0].set_ylabel('Error')
axarr[0][0].set_title('space Error')
axarr[1][0].set_ylabel('Error')
axarr[0][1].set_title('temporal Error')
axarr[1][0].set_ylabel('order')
axarr[1][0].set_title('space order')
axarr[1][1].set_ylabel('order')
axarr[1][1].set_title('temporal order')

axarr[0][0].set_ylabel('Error')
axarr[0][2].set_title('max Error')

axarr[0][3].set_title('global Error')
axarr[1][2].set_title('max order')
axarr[1][1].set_ylabel('order')
axarr[1][3].set_title('global order')
axarr[0][1].legend(legendList, frameon=False)



hxList = np.asarray(hx)
htList = np.asarray(ht)

from sympy import symbols, lambdify, latex
h_x, h_t = symbols('h_x h_t')
errorfuncs = {}

gx0, gt0 = optimize_error_cxct(Laasonen['globalError'], hxList, htList, 'Laasonen')
gx, p, gt, q = optimize_error(Laasonen['globalError'], hxList, htList, gx0, gt0, 'Laasonen')

errorexpr = gx*h_x**p + gt*h_t**q
errorX =  gx*h_x**p
errorT =  gt*h_t**q
errorfuncs['Laasonen']= [errorexpr, errorX, errorT]
gx0, gt0 = optimize_error_cxct(FTCS['globalError'], hxList, htList, 'FTCS')

gx, p, gt, q = optimize_error(FTCS['globalError'], hxList, htList, gx0, gt0, 'FTCS')

errorexpr = gx*h_x**p + gt*h_t**q
errorX =  gx*h_x**p
errorT =  gt*h_t**q
errorfuncs['FTCS']= [errorexpr, errorX, errorT]

gx0, p0 = optimize_error_gxp(Crank['globalError'], hxList, 'Crank')
print "gx0, p0: ", gx0, p0
gx, p, gt = optimize_error_gtq(Crank['globalError'], hxList, htList, 'Crank', gx0, p0 )
print "gx, p, gt: ", gx, p, gt

errorexpr = gx*h_x**p + gt*h_t**2
errorX =  gx*h_x**p
errorT =  gt*h_t**2
errorfuncs['Crank']= [errorexpr, errorX, errorT]



fig2 , axarr2 = plt.subplots(2, 3, squeeze=False)
continu = True
if continu:
    for k, solver in enumerate(solvers):

        exprList = errorfuncs[solver]
        
        errorfunc = lambdify([h_x,h_t], exprList[0], np)
        errorfuncX = lambdify([h_x], exprList[1], np)
        errorfuncT = lambdify([h_t], exprList[2], np)
        
        axarr2[0][k].plot(np.log10(eval(solver)['globalError']),lstyle[k])
        axarr2[0][k].plot(np.log10(errorfunc(hxList, htList)),lstyle[k] + '--')
        axarr2[0][k].set_title(solver)
        axarr2[0][k].legend([solver, '$'+ latex(exprList[0]) + '$'], loc = 'best', frameon=False)
        #axarr2[0][k].set_ylim(-0.01, max(eval(solver)['globalError'])*1.1)
        
        axarr2[1][k].plot(eval(solver)['globalError'],lstyle[k])
        axarr2[1][k].plot(errorfuncX(hxList),lstyle[k] + '--')
        axarr2[1][k].plot(errorfuncT(htList),lstyle[k] + ':')
        #axarr2[1][k].set_title(solver)
        axarr2[1][k].legend([solver, '$'+ latex(exprList[1]) + '$', '$'+ latex(exprList[2]) + '$'],loc = 'best',  frameon=False)
        axarr2[1][k].set_ylim(-max(eval(solver)['globalError'][-5:])*1.1, max(eval(solver)['globalError'][-5:])*1.1)

axarr2[0][0].set_ylabel('log10(E/e)')
axarr2[1][0].set_ylabel('E/e')
        
plt.show()


