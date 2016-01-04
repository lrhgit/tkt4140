# chapter7/src-ch7/test_solvers.py
import scipy 
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import numpy as np
from math import sin, sinh, pi

def laplace_Dirichlet_analytic(Ttop, x, y):
    """ Function that calculates the analytic solution to the laplace heat equation
        
            d^2(T)/dx^2 + d^2(T)/dy^2 = 0
            
        with Dirichlet boundary conditions: 
            T(x, 0) = 0
            T(0, y) = 0 y<1.5
            T(1, y) = 0 y<1.5
            T(x, 1.5) = Ttop   
        on a square grid with xmax = 1 and ymax= 1
        
        Args:
            Ttop(float): Temperature at top, T(x, 1.5)
            x(float): x-coordinate
            y(float): y-coordinate
        Returns:
            T(float): Temperature at T(x, y)
            
    """
    
    sum = 0
    
    for n in range(1,50, 1):
        lambdan = pi*n
        An = 2*((-1)**(n+1)+1)/(lambdan*sinh(1.5*lambdan))
        
        term = An*sinh(lambdan*y )*sin(lambdan*x)
        sum += term

    return Ttop*sum
    

def test_MES():
    """ Compare numerical solutions with analytic solution"""
    from laplace_heat1_gen2 import setup_LaplaceDiriclhet_yx
    
    Ttop=100

    xmax=1.0
    ymax=1.5
     
    # Set simulation parameters
    #need hx=(1/nx)=hy=(1.5/ny)
    
    Nx = 20
    
    h=xmax/Nx
    Ny = int(ymax/h)
    
    nx = Nx-1
    ny = Ny-1
    A, b  = setup_LaplaceDiriclhet_yx(Ttop, nx, ny)
    
    Temp = scipy.sparse.linalg.spsolve(A, b) # solve with sparse solvers
    X = np.linspace(h, xmax - h, nx)
    Y = np.linspace(h, ymax - h, ny)

    TempAnalytic = np.zeros_like(Temp)
    i = 0
    indice = []
    for x in X:
        for y in Y:
            TempAnalytic[i] = laplace_Dirichlet_analytic(Ttop, x, y)
            i += 1
            indice.append([x, y])
    
        
    assert np.linalg.norm(TempAnalytic-Temp, np.inf)<1

def test_XY_vs_YX_numbering():
    """ Verify that the different numbering (in laplace Dirichlet example); with respect to x, and with respect to y
        give the same result
    """
    from laplace_heat1_gen2 import setup_LaplaceDiriclhet_yx, setup_LaplaceDiriclhet_xy
    from Visualization import convert_xy_yx
    # Set temperature at the top
    Ttop=100
    
    xmax=1.0
    ymax=1.5
     
    # Set simulation parameters
    #need hx=(1/nx)=hy=(1.5/ny)
    
    Nx = 4
    h=xmax/Nx
    Ny = int(ymax/h)
    
    nx = Nx-1
    ny = Ny-1
        
    A_yx, b_yx  = setup_LaplaceDiriclhet_yx(Ttop, nx, ny)
    A_xy, b_xy = setup_LaplaceDiriclhet_xy(Ttop, nx, ny)
    
    TempXY = scipy.sparse.linalg.spsolve(A_xy, b_xy) # solve with sparse solvers
    TempYX = scipy.sparse.linalg.spsolve(A_yx, b_yx) 
    
    TempXY_converted = convert_xy_yx(TempXY, nx, ny) # convert Temperature array from xy format to yx format
    
    tol = 1e-12
    assert np.linalg.norm(TempXY_converted-TempYX, np.inf)

    



