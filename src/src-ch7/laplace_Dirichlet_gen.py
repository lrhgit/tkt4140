# src-ch7/laplace_Diriclhet_gen.py
import numpy as np
import scipy 
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pylab as plt
import time
from math import sinh, sin, pi

#import matplotlib.pyplot as plt

# Change some default values to make plots more readable on the screen
LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT


def setup_LaplaceDiriclhet_yx(Ttop, nx, ny):
    """ Function that returns A matrix and b vector  of the laplace Dirichlet heat problem A*T=b using central differences
        and assuming dx=dy, based on numbering with respect to y-dir, e.g:
        
            100 100 100 100             -4  1  0  1  0  0          0
            0   T3  T6  0                1 -4  1  0  1  0          0
            0   T2  T5  0                0  1 -4  0  0  1         -100
            0   T1  T4  0     --> A =    1  0  0 -4  1  0   , b =  0
            0   0   0   0                0  1  0  1 -4  1          0
                                         0  0  1  0  1 -4         -100
            
            T = [T1, T2, T3, T4, T5, T6]^T
            
        Args:
            nx(int): number of elements in each row in the grid, nx=2 in the example above
            ny(int): number of elements in each column in the grid, ny=3 in the example above
        Returns:
            A(matrix): Sparse matrix A, in the equation A*T = b
            b(array): RHS, of the equation A*T = b
    """
    n = (nx)*(ny) #number of unknowns
    d = np.ones(n) # diagonals
    b = np.zeros(n) #RHS
    
    d0 = d.copy()*-4
    d1 = d[0:-1].copy()
    dny = d[0:-ny].copy()
    
    d1[ny-1::ny] = 0 # every ny element on first diagonal is zero
    
    b[ny-1::ny] = - Ttop
    
    A = scipy.sparse.diags([d0, d1, d1, dny, dny], [0, 1, -1, ny, -ny], format='csc')
    return A, b


def setup_LaplaceDiriclhet_xy(Ttop, nx, ny):
    """ Function that returns A matrix and b vector  of the laplace Dirichlet heat problem A*T=b using central differences
        and assuming dx=dy, based on numbering with respect to x-dir, e.g:
        
            100 100 100 100             -4  1  1  0  0  0          0
            0   T5  T6  0                1 -4  0  1  0  0          0
            0   T3  T4  0                1  0 -4  1  1  0          0
            0   T1  T2  0     --> A =    0  1  1 -4  0  1    ,b =  0
            0   0   0   0                0  0  1  0 -4  1         -100
                                         0  0  0  1  1 -4         -100
            
            T = [T1, T2, T3, T4, T5, T6]^T
            
        Args:
            nx(int): number of elements in each row in the grid, nx=2 in the example above
            ny(int): number of elements in each column in the grid, ny=3 in the example above
        Returns:
            A(matrix): Sparse matrix A, in the equation A*T = b
            b(array): RHS, of the equation A*t = b
    """
    
    n = (nx)*(ny) #number of unknowns
    d = np.ones(n) # diagonals
    b = np.zeros(n) #RHS
    
    d0 = d.copy()*-4
    d1 = d.copy()[0:-1]
    dnx = d.copy()[0:-nx]
    
    d1[nx-1::nx] = 0 # every nx element on first diagonal is zero
    
    b[-nx:] = -Ttop

    A = scipy.sparse.diags([d0, d1, d1, dnx, dnx], [0, 1, -1, nx, -nx], format='csc')
    
    return A, b

    
if __name__ == '__main__':
    
    from Visualization import plot_Surface_yx, visualize_setup_yx, visualize_setup_xy, convert_xy_yx
    
    # Main program
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
    print "\n"
    visualize_setup_yx(TempYX, nx, ny, Ttop) # visualize the numbering of the XY setup on a printed grid, with corresponding results
    print "\n"
    visualize_setup_xy(TempXY, nx, ny, Ttop) # visualize the numbering of the YX setup on a printed grid, with corresponding results
    
    max_difference = max(np.asarray(TempXY_converted - TempYX))
    print "\n"
    print "max(abs(TempXY-TempYX)) < 10^-10: ", max_difference<1e-10
    print "np.linalg.norm(TempXY-TempYX, np.inf): ", np.linalg.norm(TempXY_converted-TempYX, np.inf)
    print "\n"
    plot_Surface_yx(TempYX, Ttop, xmax, ymax, Nx, Ny, nx, ny)
    plt.show()
    



