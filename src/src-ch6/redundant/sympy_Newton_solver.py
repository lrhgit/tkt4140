from sympy import (
    symbols,   # define mathematical symbols for symbolic computing
    diff,      # differentiate expressions
    integrate, # integrate expressions
    Rational,  # define rational numbers
    lambdify,  # turn symbolic expressions into Python functions
    cos,
    sin,
    exp,
    Matrix 
    )
import numpy as np
def Newton_solver_All(error, hx, ht, x0):
    
    gx, p, gt, q  = symbols('gx p gt q')
    hx1, hx2, hx3, hx4, ht1, ht2, ht3, ht4 = symbols('hx1 hx2 hx3 hx4 ht1 ht2 ht3 ht4')
    epsilon1 = gx*hx1**p + gt*ht1**q - error[-4]
    epsilon2 = gx*hx2**p + gt*ht2**q - error[-3]
    epsilon3 = gx*hx3**p + gt*ht3**q - error[-2]
    epsilon4 = gx*hx4**p + gt*ht4**q - error[-1]

    epsilon = [epsilon1, epsilon2, epsilon3, epsilon4]
    x = [gx, p, gt, q]
    knowns = [hx1, hx2, hx3, hx4, ht1, ht2, ht3, ht4]
    
    def f(i,j):
        return diff(epsilon[i], x[j])
        
    Jacobi = Matrix(4, 4, f)

    JacobiInv = Jacobi.inv()
    epsilonfunc = lambdify([x, knowns], epsilon)
    JacobiInvfunc = lambdify([x, knowns], JacobiInv)
    x = x0
    knowns = [hx[-4], hx[-3], hx[-2], hx[-1], ht[-4], ht[-3], ht[-2], ht[-1]]

    
    print "knowns: ", knowns
    
    try:
        Jinv = np.asarray(JacobiInvfunc(x, knowns))
    except:
        print "print unable to invert jacobi matrix"
        Jinv = None
    F = np.asarray(epsilonfunc(x, knowns))
    
    for n in range(8):
        if Jinv != None:
            xnew = x - np.dot(Jinv, F)
            F = np.asarray(epsilonfunc(xnew, knowns))
            try:
                Jinv = np.asarray(JacobiInvfunc(x, knowns))
            except:
                print " unable to invert jacobi matrix"
                Jinv = None
            
            x = xnew
            print "n, x: ", n, x

def Newton_solver_gxgt(error, hx, ht, x0, p_value=1, q_value=1):
    
    gx, gt, p, q, hx1, hx2, ht1, ht2  = symbols('gx gt p q qhx1 hx2 ht1 ht2')
    epsilon1 = gx*hx1**p + gt*ht1**q - error[-2]
    epsilon2 = gx*hx2**p + gt*ht2**q - error[-1]


    epsilon = [epsilon1, epsilon2]
    x = [gx, gt]
    knowns = [p, q, hx1, hx2, ht1, ht2]
    
    def f(i,j):
        return diff(epsilon[i], x[j])
        
    Jacobi = Matrix(2, 2, f)

    JacobiInv = Jacobi.inv()
    epsilonfunc = lambdify([x, knowns], epsilon)
    JacobiInvfunc = lambdify([x, knowns], JacobiInv)
    x = x0
    knowns = [p_value, q_value, hx[-2], hx[-1], ht[-2], ht[-1]]
    F = np.asarray(epsilonfunc(x, knowns))
    
    for n in range(8):
        Jinv = np.asarray(JacobiInvfunc(x, knowns))
        F = np.asarray(epsilonfunc(x, knowns))
        x = x - np.dot(Jinv, F)

        print "n, x: ", n, x
    
    return x

def Newton_solver_pq(error, hx, ht, gx, gt, x0):
    
    p, q  = symbols('p q')
    
    epsilon1 = gx*hx[-2]**p + gt*ht[-2]**q - error[-2]
    epsilon2 = gx*hx[-1]**p + gt*ht[-1]**q - error[-1]


    epsilon = [epsilon1, epsilon2]
    x = [p, q]
    
    def f(i,j):
        return diff(epsilon[i], x[j])
        
    Jacobi = Matrix(2, 2, f)

    JacobiInv = Jacobi.inv()
    epsilonfunc = lambdify(x, epsilon)
    JacobiInvfunc = lambdify(x, JacobiInv)
    x = x0
    
    
    

    F = np.asarray(epsilonfunc(x))
    
    for n in range(8):
        Jinv = np.asarray(JacobiInvfunc(x))
        F = np.asarray(epsilonfunc(x))
        x = x - np.dot(Jinv, F)

    return x

def optimize_error_cxct(error, hx, ht, scheme):
    from scipy.optimize import curve_fit
    
    if scheme=='macCormack':
        p = 2
        q = 2
    elif scheme=='lax_friedrich':
        p = 2
        q = 1
    elif scheme=='ftbs':
        p = 1
        q = 1
    
    ht_hx = ht[0]/hx[0]
    def func(x, gx, gt):
        dt = ht_hx*x
        return gx*x**p + gt*dt**q
    
    
    x0 = np.array([1,2])
    xdata = hx
    ydata = error
    cx, ct = curve_fit(func, xdata, ydata, x0)[0]

    return cx, ct

def optimize_error(error, hx, ht, gx0, gt0, scheme):
    from scipy.optimize import curve_fit
    
    if scheme=='macCormack':
        p0 = 2
        q0 = 2
    elif scheme=='lax_friedrich':
        p0 = 2
        q0 = 1
    elif scheme=='ftbs':
        p0 = 1
        q0 = 1
    
    ht_hx = ht[0]/hx[0]
    def func(x, gx, p, gt, q):
        dt = ht_hx*x
        return gx*x**p + gt*dt**q
    
    
    x0 = np.array([gx0, p0, gt0, q0])
    xdata = hx
    ydata = error
    gx, p, gt, q = curve_fit(func, xdata, ydata, x0)[0]
    gx, p, gt, q = round(gx,1), round(p, 2), round(gt,1), round(q, 2)
    return gx, p, gt, q


