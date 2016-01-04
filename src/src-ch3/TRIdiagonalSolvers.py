# chapter3/src-ch3/TRIdiagonalSolvers.py

import numpy as np
from matplotlib.pyplot import *

# change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT
font = {'size' : 16}; rc('font', **font)


def tdma(a, b, c, d):
    """Solution of a linear system of algebraic equations with a
        tri-diagonal matrix of coefficients using the Thomas-algorithm.

    Args:
        a(array): an array containing lower diagonal (a[0] is not used)
        b(array): an array containing main diagonal 
        c(array): an array containing lower diagonal (c[-1] is not used)
        d(array): right hand side of the system
    Returns:
        x(array): solution array of the system
    
    """
    
    n = len(b)
    x = np.zeros(n)
    
    # elimination:
    
    for k in range(1,n):
        q = a[k]/b[k-1]
        b[k] = b[k] - c[k-1]*q
        d[k] = d[k] - d[k-1]*q
    
    # backsubstitution:
    
    q = d[n-1]/b[n-1]
    x[n-1] = q
    
    for k in range(n-2,-1,-1):
        q = (d[k]-c[k]*q)/b[k]
        x[k] = q
    
    
    return x

def tripiv(a, b, c, d):
    """Solution of a linear system of algebraic equations with a
        tri-diagonal matrix of coefficients using the Thomas-algorithm with pivoting.

    Args:
        a(array): an array containing lower diagonal (a[0] is not used)
        b(array): an array containing main diagonal 
        c(array): an array containing lower diagonal (c[-1] is not used)
        d(array): right hand side of the system
    Returns:
        x(array): solution array of the system
    
    """
    
    n = len(b)
    x = np.zeros(n)
    fail = 0
    
    # reordering
    
    a[0] = b[0]
    b[0] = c[0]
    c[0] = 0
    
    # elimination:
    
    l = 0
    
    for k in range(0,n):
        q = a[k]
        i = k
        
        if l < n-1:
            l = l + 1
    
        for j in range(k+1,l+1):
            q1 = a[j]
            if (np.abs(q1) > np.abs(q)):
                q = q1
                i = j
        if q == 0:
            fail = -1
        
        if i != k:
            q = d[k]
            d[k] = d[i]
            d[i] = q
            q = a[k]
            a[k] = a[i]
            a[i] = q
            q = b[k]
            b[k] = b[i]
            b[i] = q
            q = c[k]
            c[k] =c[i]
            c[i] = q
        for i in range(k+1,l+1):
            q = a[i]/a[k]
            d[i] = d[i]-q*d[k]
            a[i] = b[i]-q*b[k]
            b[i] = c[i]-q*c[k]
            c[i] = 0
            

    
    # backsubstitution
    
    x[n-1] = d[n-1]/a[n-1]
    x[n-2] = (d[n-2]-b[n-2]*x[n-1])/a[n-2]
    
    for i in range(n-3,-1,-1):
        
        q = d[i] - b[i]*x[i+1]
        x[i] = (q - c[i]*x[i+2])/a[i]
    
    return x

if __name__ == '__main__':
    

    
    xanalytic = np.array([1., 2., 3.])
 
    def test_tdma():
        """Unit test of tdma and tripivequation of 3x3 system with known solution of x1 = 1, x2 = 2, x3 = 3"""
        a = np.array([0,1,1])
        b = np.array([1,1,1])
        c = np.array([2,1,0])
        d = np.array([5,6,5])
        x = tdma(a,b,c,d)
        print "x= ", x
        a = np.array([0,1,1])
        b = np.array([1,1,1])
        c = np.array([2,1,0])
        d = np.array([5,6,5])
        x2 = tripiv(a, b, c, d)
        print "x2= ", x2
        
        return x, x2
    
    

    x,x2 = test_tdma()

    print "xtdma - xanalytic ",x-xanalytic
    print "xtripiv - xanalytic ",x2-xanalytic
    
    
    
    