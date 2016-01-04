import numpy as np
import scipy 
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pylab as plt
import time
from math import sinh

#import matplotlib.pyplot as plt

# Change some default values to make plots more readable on the screen
LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT

def setup_LaplaceNeumann_xy(Ttop, Tright, nx, ny):
    """ Function that returns A matrix and b vector  of the laplace Neumann heat problem A*T=b using central differences
        and assuming dx=dy, based on numbering with respect to x-dir, e.g:
        
                1   1   1             -4  2  2  0  0  0          0
            T6  T5  T6  0              1 -4  0  2  0  0          0
            T4  T3  T4  0              1  0 -4  2  1  0          0
            T2  T1  T2  0     --> A =  0  1  1 -4  0  1    ,b =  0
                T3  T4                 0  0  1  0 -4  2         -1
                                       0  0  0  1  1 -4         -1
            
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
    d1_lower = d.copy()[0:-1]
    d1_upper = d1_lower.copy()
    
    dnx_lower = d.copy()[0:-nx]
    dnx_upper = dnx_lower.copy()
    
    d1_lower[nx-1::nx] = 0 # every nx element on first diagonal is zero; starting from the nx-th element
    d1_upper[nx-1::nx] = 0
    d1_upper[::nx] = 2 # every nx element on first upper diagonal is two; stating from the first element. 
                 # this correspond to all equations on border (x=0, y)
    
    dnx_upper[0:nx] = 2 # the first nx elements in the nx-th upper diagonal is two; 
                        # This correspond to all equations on border (x, y=0) 
    
    b[-nx:] = -Ttop
    b[nx-1::nx] += -Tright

    A = scipy.sparse.diags([d0, d1_upper, d1_lower, dnx_upper, dnx_lower], [0, 1, -1, nx, -nx], format='csc')
    
    return A, b

if __name__ == '__main__':
    from Visualization import plot_SurfaceNeumann_xy
    # Main program
    # Set temperature at the top
    Ttop=1
    Tright = 0.0
    xmax=1.0
    ymax=1.
     
    # Set simulation parameters
    #need hx=(1/nx)=hy=(1.5/ny)
    
    Nx = 10
    h=xmax/Nx
    Ny = int(ymax/h)
    
    A, b = setup_LaplaceNeumann_xy(Ttop, Tright, Nx, Ny)
    
    Temp = scipy.sparse.linalg.spsolve(A, b)
    
    plot_SurfaceNeumann_xy(Temp, Ttop, Tright, xmax, ymax, Nx, Ny)

#    figfile='LaPlace_vNeumann.png'
#    plt.savefig(figfile, format='png',transparent=True)
    plt.show()


