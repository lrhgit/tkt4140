# chapter6/src-ch6/advection_optimize_error.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sympy import symbols, lambdify, latex

def optimize_error_cxct(errorList, hxList, 
                        htList, p=1.0, q=1.0):
    """ Function that optimimze the values Cx and Ct in the expression
        E = epsilon = Cx hx^p + Ct ht^q, assuming p and q are known, 
        using scipy's optimization tool curve_fit
    
        Args:
            errorList(array): array of calculated numerical discretization errors E
            hxList(array): array of spatial step lengths corresponding to errorList
            htList(array): array of temporal step lengths corresponding to errorList
            p (Optional[float]): spatial order. Assumed to be equal to theoretical value
            q (Optional[float]): temporal order. Assumed to be equal to theoretical value

        Returns:
            Cx0(float): the optimized values of Cx
            Ct0(float): the optimized values of Ct
    """
    
    def func(h, Cx, Ct):
        """ function to be matched with ydata:
            The function has to be on the form func(x, parameter1, parameter2,...,parametern)
            where where x is the independent variable(s), and parameter1-n are the parameters to be optimized.
        """
        return Cx*h[0]**p + Ct*h[1]**q
    
    
    x0 = np.array([1,2]) # initial guessed values for Cx and Ct
    xdata = np.array([hxList, htList])# independent variables
    ydata = errorList # data to be matched with expression in func
    Cx0, Ct0 = curve_fit(func, xdata, ydata, x0)[0] # call scipy optimization tool curvefit

    return Cx0, Ct0

def optimize_error(errorList, hxList, htList, 
                   Cx0, p0, Ct0, q0):
    """ Function that optimimze the values Cx, p Ct and q in the expression
        E = epsilon = Cx hx^p + Ct ht^q, assuming p and q are known, 
        using scipy's optimization tool curve_fit
    
        Args:
            errorList(array): array of calculated numerical discretization errors E
            hxList(array): array of spatial step lengths corresponding to errorList
            htList(array): array of temporal step lengths corresponding to errorList
            Cx0 (float): initial guessed value of Cx
            p (float): initial guessed value of p
            Ct0 (float): initial guessed value of Ct
            q (float): initial guessed value of q

        Returns:
            Cx(float): the optimized values of Cx
            p (float): the optimized values of p
            Ct(float): the optimized values of Ct
            q (float): the optimized values of q
    """
    
    
    def func(h, gx, p, gt, q):
        """ function to be matched with ydata:
            The function has to be on the form func(x, parameter1, parameter2,...,parametern)
            where where x is the independent variable(s), and parameter1-n are the parameters to be optimized.
        """
        return gx*h[0]**p + gt*h[1]**q
    
    
    x0 = np.array([Cx0, p0, Ct0, q0]) # initial guessed values for Cx, p, Ct and q
    xdata = np.array([hxList, htList]) # independent variables
    ydata = errorList # data to be matched with expression in func
    
    gx, p, gt, q = curve_fit(func, xdata, ydata, x0)[0] # call scipy optimization tool curvefit
    gx, p, gt, q = round(gx,2), round(p, 2), round(gt,2), round(q, 2)
    
    return gx, p, gt, q

# Program starts here:

# empty lists, to be filled in with values from 'advection_scheme_errors.txt'
hxList = [] 
htList = []

E_ftbs = []
E_lax_friedrich = []
E_lax_wendroff = []

lineNumber = 1

with open('advection_scheme_errors.txt', 'r') as FILENAME:
    """ Open advection_scheme_errors.txt for reading.
        structure of file:
        hx    ht    E_ftbs    E_lax_friedrich    E_lax_wendroff
        
        with the first line containing these headers, and the next lines containing
        the corresponding values.
    """
    # iterate all lines in FILENAME:
    for line in FILENAME:
        if lineNumber ==1:
            # skip first line which contain headers
            lineNumber += 1
        else:
            lineList = line.split() # sort each line in a list: lineList = [hx, ht, E_ftbs, E_lax_friedrich, E_lax_wendroff]
            
            # add values from this line to the lists
            hxList.append(float(lineList[0]))
            htList.append(float(lineList[1]))
            
            E_ftbs.append(float(lineList[2]))
            E_lax_friedrich.append(float(lineList[3]))
            E_lax_wendroff.append(float(lineList[4]))
            
            lineNumber += 1
      
FILENAME.close()

# convert lists to numpy arrays:
hxList = np.asarray(hxList) 
htList = np.asarray(htList) 

E_ftbs = np.asarray(E_ftbs) 
E_lax_friedrich = np.asarray(E_lax_friedrich) 
E_lax_wendroff = np.asarray(E_lax_wendroff) 


ErrorList = [E_ftbs, E_lax_friedrich, E_lax_wendroff]
schemes = ['ftbs', 'lax_friedrich', 'lax_wendroff']

p_theoretical = [1, 2, 2] # theoretical spatial orders
q_theoretical = [1, 1, 2] # theoretical temporal orders

h_x, h_t = symbols('h_x h_t')

XtickList = [i for i in range(1, len(hxList)+1)]
Xticknames = [r'$(h_x , h_t)_{0}$'.format(str(i)) for i in range(1, len(hxList)+1)]
lstyle = ['b', 'r', 'g', 'm']
legendList = []

# Optimize Cx, p, Ct, q for every scheme
for k, E in enumerate(ErrorList):
    
    # optimize using only last 5 values of E and h for the scheme, as the first values may be outside asymptotic range:
    
    Cx0, Ct0 = optimize_error_cxct(E[-5:], hxList[-5:], htList[-5:], 
                                   p=p_theoretical[k], q=q_theoretical[k]) # Find appropriate initial guesses for Cx and Ct 
    
    Cx, p, Ct, q = optimize_error(E[-5:], hxList[-5:], htList[-5:], 
                                  Cx0, p_theoretical[k], Ct0, q_theoretical[k]) # Optimize for all parameters Cx, p, Ct, q

    # create sympy expressions of e, ex and et:
    errorExpr = Cx*h_x**p + Ct*h_t**q
    
    print errorExpr
    errorExprHx = Cx*h_x**p 
    errorExprHt = Ct*h_t**q
    
    # convert e, ex and et to python functions:
    errorFunc = lambdify([h_x, h_t], errorExpr)
    errorFuncHx = lambdify([h_x], errorExprHx)
    errorFuncHt = lambdify([h_t], errorExprHt)
    
    # plotting:
    fig , ax = plt.subplots(2, 1, squeeze=False)
    
    ax[0][0].plot(XtickList, np.log10(E),lstyle[k])
    ax[0][0].plot(XtickList, np.log10(errorFunc(hxList, htList)),lstyle[k] + '--')
    
    ax[1][0].plot(XtickList[-5:], E[-5:],lstyle[k])
    ax[1][0].plot(XtickList[-5:], errorFunc(hxList, htList)[-5:],lstyle[k] + '--')
    ax[1][0].plot(XtickList[-5:], errorFuncHx(hxList[-5:]), lstyle[k] + '-.')
    ax[1][0].plot(XtickList[-5:], errorFuncHt(htList[-5:]),lstyle[k] + ':')
    
    # legends and labels:
    
    e_x_latex = str(Cx) + '\,h_x^'+ '{' + str(p) + '}\,'
    e_t_latex = str(Ct) + '\,h_t^'+ '{' + str(q) + '}\,'
    ax[0][0].set_title(schemes[k])
    if schemes[k]== 'macCormack':
        add = ' '
    else:
        add = ' + '
    ax[0][0].legend(['$\epsilon = E$', '$\epsilon = '  + e_x_latex + add +  e_t_latex +'$'], loc = 'best', frameon=False)
    #ax[0][0].set_ylim(-min(np.log10(E)), max(np.log10(E)))
    ax[0][0].set_xticks(XtickList)
    ax[0][0].set_xlabel('$(h_x , h_t)_n$')
    #ax[0][0].set_xticklabels(Xticknames)
    
    

    #axarr2[1][k].set_title(solver)
    ax[1][0].legend(['$\epsilon = E$', '$\epsilon = '  + e_x_latex + add +  e_t_latex +'$', '$\epsilon =' + e_x_latex+ '$', '$\epsilon =' + e_t_latex+ '$'], frameon=False)
    ax[1][0].set_xticks(XtickList[-5:])
    ax[1][0].set_xlabel('$(h_x , h_t)_n$')

    ax[0][0].set_ylabel(r'$log_{10}(\epsilon)$')
    ax[1][0].set_ylabel('$\epsilon$')
    fig.tight_layout()
#    plt.savefig('../figs/Optimize_errors_'+str(schemes[k])+'.png', transparent=True, frameon=False) # transparent=True
    # Convert to Bokeh
#     import bokeh.mpl, bokeh.plotting as bpl
#     p = bokeh.mpl.to_bokeh(notebook=False, xkcd=False)
#     #p = bokeh.mpl.to_bokeh()
#     bpl.output_file('tmp.html', mode='cdn')
#     bpl.save(p)
    #bpl.show(p)


plt.show()
 
        



