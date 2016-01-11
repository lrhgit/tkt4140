
import numpy as np
import matplotlib.pyplot as plt
from sympy_Newton_solver import Newton_solver_All, Newton_solver_gxgt, Newton_solver_pq

macCormack = {}
lax_friedrich ={}
ftbs = {}
solvers = ['macCormack', 'lax_friedrich', 'ftbs']
for k, solver in enumerate(solvers):
    hx = []
    ht = []
    spaceError = []
    temporalError = []
    maxError = []
    globalError = []
    
    with open(solver + '_errors.txt', 'r') as filename:
        for line in filename:
            #print line
            splitline = line.split()
            stripline = line.strip()
            try:
                if splitline[0] == 'hx:':
                    hx.append(float(splitline[1][1:-1]))
                    for k in range(2, 2 + len(splitline[2:])):
                        hx.append(float(splitline[k][:-1]))
                elif splitline[0] == 'ht:':
                    ht.append(float(splitline[1][1:-1]))
                    for k in range(2, 2 + len(splitline[2:])):
                        ht.append(float(splitline[k][:-1]))
                                  
                elif splitline[0] == 'spaceError:':
                    spaceError.append(float(splitline[1][1:-1]))
                    for k in range(2, 2 + len(splitline[2:])):
                        spaceError.append(float(splitline[k][:-1]))
    
                elif splitline[0] == 'temporalError:':
                    temporalError.append(float(splitline[1][1:-1]))
                    for k in range(2, 2 + len(splitline[2:])):
                        temporalError.append(float(splitline[k][:-1]))
    
                elif splitline[0] == 'maxError:':

                    maxError.append(float(splitline[1][1:-1]))
                    for k in range(2, 2 + len(splitline[2:])):
                        maxError.append(float(splitline[k][:-1]))
    
                elif splitline[0] == 'globalError:':

                    globalError.append(float(splitline[1][1:-1]))
                    for k in range(2, 2 + len(splitline[2:])):
                        globalError.append(float(splitline[k][:-1]))
                    
            except:
                print "empty line"
    rx = hx[0]/hx[1]
    rt = ht[0]/ht[1]
    eval(solver)['spaceError'] =  spaceError
    eval(solver)['spaceOrder'] = np.log(np.asarray(spaceError[:-1])/np.asarray(spaceError[1:]))/np.log(rx)
    eval(solver)['temporalError'] =  temporalError
    eval(solver)['temporalOrder'] = np.log(np.asarray(temporalError[:-1])/np.asarray(temporalError[1:]))/np.log(rt)
    eval(solver)['maxError'] =  maxError
    eval(solver)['maxOrder'] = np.log(np.asarray(maxError[:-1])/np.asarray(maxError[1:]))/np.log(rx)
    eval(solver)['globalError'] =  globalError  
    eval(solver)['globalOrder'] = np.log(np.asarray(globalError[:-1])/np.asarray(globalError[1:]))/np.log(rx)
    
print "\n"
print macCormack['globalError']
print "\n"
print lax_friedrich['globalError']
print "\n"
print ftbs['globalError']
print "\n"
print hx
print ht
from sympy_Newton_solver import optimize_error_cxct, optimize_error

#plt.plot(curvefiterror, 'b--')
#x = Newton_solver_gxgt(error, hx, ht, [cx, ct], p_value=2, q_value=2)
#Newton_solver_All(macCormack['globalError'], hx, ht, [cx, 2, ct, 2])

print "\n"


#x = Newton_solver_gxgt(macCormack['globalError'], hx, ht, [1, 1], p_value=2, q_value=2)
hxList = np.asarray(hx)
htList = np.asarray(ht)
legendList = []
lstyle = ['b', 'r', 'g', 'm']
fig , axarr = plt.subplots(2, 4, squeeze=False)
fig2 , axarr2 = plt.subplots(2, 3, squeeze=False)

from sympy import symbols, lambdify, latex
h_x, h_t = symbols('h_x h_t')
errorfuncs = {}
legends = []

Ntds= len(eval(solver)['spaceError'])

for k, solver in enumerate(solvers):
    if solver=='macCormack':
        theoretical = 2

    if solver=='lax_friedrich':
        theoretical = 1

    if solver=='ftbs':
        theoretical = 1
    gx0, gt0 = optimize_error_cxct(eval(solver)['globalError'][-5:], hxList[-5:], htList[-5:], solver)
    gx, p, gt, q = optimize_error(eval(solver)['globalError'][-5:], hxList[-5:], htList[-5:], gx0, gt0, solver)
    
    errorexpr = gx*h_x**p + gt*h_t**q
    errorX = lambdify([h_x], gx*h_x**p, np)
    errort = lambdify([h_t], gt*h_t**q, np)
    errorfuncs[solver] = errorexpr
    tempfunc = lambdify([h_x, h_t], errorexpr, np)
    errorcalc = tempfunc(hxList, htList)

    axarr2[0][k].plot(np.log10(np.asarray(eval(solver)['globalError'])),lstyle[k])
    axarr2[0][k].plot(np.log10(np.asarray(errorcalc)),lstyle[k] + '--')
    axarr2[0][k].set_title(solver)
    axarr2[0][k].legend([solver, '$'+ latex(errorexpr) + '$'], loc = 'best', frameon=False)
    axarr2[0][k].set_ylim(min(np.log10(np.asarray(eval(solver)['globalError']))), max(np.log10(np.asarray(eval(solver)['globalError']))))
    
    axarr2[1][k].plot(np.log10(np.asarray(eval(solver)['globalError'][-5:])),lstyle[k])
    axarr2[1][k].plot(np.log10(abs(np.asarray(errorX(hxList[-5:])))),lstyle[k] + '--')
    axarr2[1][k].plot(np.log10(abs(np.asarray(errort(htList[-5:])))),lstyle[k] + ':')
    #axarr2[1][k].set_title(solver)
    axarr2[1][k].legend([solver, '$'+ latex(gx*h_x**p) + '$', '$'+ latex(gt*h_t**q) + '$'], frameon=False)
    #axarr2[1][k].set_ylim(-max(eval(solver)['globalError'][-5:])*1.1, max(eval(solver)['globalError'][-5:])*1.1)

#     axarr2[2][k].plot(eval(solver)['globalOrder'],lstyle[k])
#     axarr2[2][k].plot(np.log(tempfunc(hxList[:-1], htList[:-1])/tempfunc(hxList[1:], htList[1:]))/np.log(rx),lstyle[k] + '--')
#     axarr2[2][k].axhline(theoretical, xmin=0, xmax=Ntds-2, linestyle=':', color='k')
#     
#     #axarr2[1][k].set_title(solver)
#     axarr2[2][k].legend([solver, '$'+ latex(errorexpr) + '$', 'theoretical'], frameon=False)
#     axarr2[2][k].set_ylim(0, 3)
    
    axarr[0][0].plot(eval(solver)['spaceError'],lstyle[k])
    axarr[0][1].plot(eval(solver)['temporalError'],lstyle[k])
    axarr[0][2].plot(eval(solver)['maxError'],lstyle[k])
    axarr[0][3].plot(eval(solver)['globalError'],lstyle[k])
    axarr[1][0].plot(eval(solver)['spaceOrder'],lstyle[k])
    axarr[1][1].plot(eval(solver)['temporalOrder'],lstyle[k])
    axarr[1][2].plot(eval(solver)['maxOrder'],lstyle[k])
    axarr[1][3].plot(eval(solver)['globalOrder'],lstyle[k])
    print "\n"
    print solver
    print eval(solver)['globalOrder']
    legendList.append(solver)


axarr2[0][0].set_ylabel('rms Error')
axarr2[1][0].set_ylabel('rms Error')





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
plt.show()

print 
        



