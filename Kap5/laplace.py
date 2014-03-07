import numpy
#from numpy import arange    

import matplotlib.pyplot as plt
import matplotlib.mlab as ml
from scipy import weave
import scipy as sc
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg


class Grid:
    """A simple grid class that stores the details and solution of the
    computational grid."""
    def __init__(self, nx=10, ny=10, xmin=0.0, xmax=1.0,
                 ymin=0.0, ymax=1.0):
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.dx = float(xmax-xmin)/(nx-1)
        self.dy = float(ymax-ymin)/(ny-1)
        self.u = numpy.zeros((nx, ny), 'd')
        self.nx = nx
        self.ny = ny
        # used to compute the change in solution in some of the methods.
        self.old_u = self.u.copy()

    def setBCFunc(self, func):
        """Sets the BC given a function of two variables."""
        xmin, ymin = self.xmin, self.ymin
        xmax, ymax = self.xmax, self.ymax
        x = numpy.arange(xmin, xmax + self.dx*0.5, self.dx)
        y = numpy.arange(ymin, ymax + self.dy*0.5, self.dy)
        self.u[0 ,:] = func(xmin,y)
        self.u[-1,:] = func(xmax,y)
        self.u[:, 0] = func(x,ymin)
        self.u[:,-1] = func(x,ymax)

    def setBC(self, l, r, b, t):        
        """Sets the boundary condition given the left, right, bottom
        and top values (or arrays)"""        
        self.u[0, :] = l
        self.u[-1, :] = r
        self.u[:, 0] = b
        self.u[:,-1] = t
        self.old_u = self.u.copy()

    def computeError(self):
        """Computes absolute error using an L2 norm for the solution.
        This requires that self.u and self.old_u must be appropriately
        setup."""
        v = (self.u - self.old_u).flat
        return numpy.sqrt(numpy.dot(v,v))


class LaplaceSolver:
    """A simple Laplacian solver that can use different schemes to
    solve the problem."""
    def __init__(self, grid, stepper='numeric'):
        self.grid = grid
        self.setTimeStepper(stepper)

    def slowTimeStep(self, dt=0.0):
        """Takes a time step using straight forward Python loops."""
        g = self.grid
        nx, ny = g.u.shape
        dx2, dy2 = g.dx**2, g.dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        u = g.u

        err = 0.0
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                tmp = u[i,j]
                u[i,j] = ((u[i-1, j] + u[i+1, j])*dy2 +
                          (u[i, j-1] + u[i, j+1])*dx2)*dnr_inv
                diff = u[i,j] - tmp
                err += diff*diff
#        g.u = u
        return numpy.sqrt(err)
    
    def numericTimeStep(self, dt=0.0):
        """Takes a time step using a NumPy expression."""
        g = self.grid
        dx2, dy2 = g.dx**2, g.dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        u = g.u
        g.old_u = u.copy() # needed to compute the error.

        # The actual iteration
        u[1:-1, 1:-1] = ((u[0:-2, 1:-1] + u[2:, 1:-1])*dy2 +
                         (u[1:-1,0:-2] + u[1:-1, 2:])*dx2)*dnr_inv

        return g.computeError()
    
    def blitzTimeStep(self, dt=0.0):
        """Takes a time step using a NumPy expression that has been
        blitzed using weave."""
        g = self.grid
        dx2, dy2 = g.dx**2, g.dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        u = g.u
        g.old_u = u.copy()

        # The actual iteration
        expr = "u[1:-1, 1:-1] = ((u[0:-2, 1:-1] + u[2:, 1:-1])*dy2 + "\
               "(u[1:-1,0:-2] + u[1:-1, 2:])*dx2)*dnr_inv"
        weave.blitz(expr, check_size=0)

        return g.computeError() 

    def setTimeStepper(self, stepper='numeric'):
        """Sets the time step scheme to be used while solving given a
        string which should be one of ['slow', 'numeric', 'blitz',
        'inline', 'fastinline', 'fortran']."""
        if stepper == 'slow':
            self.timeStep = self.slowTimeStep
        elif stepper == 'blitz':
             self.timeStep = self.blitzTimeStep
        else:
            self.timeStep = self.numericTimeStep

    def solve(self, n_iter=0, eps=1.0e-16):
        err = self.timeStep()
        count = 1
        returnValue = None

        while err > eps and count < n_iter:
            err = self.timeStep()
            count = count + 1
            
        returnValue = (count, err)
        return returnValue

def BC(x, y):    
    """Used to set the boundary condition for the grid of points.
    Change this as you feel fit."""    
    return (x**2 - y**2)


## Set up simulation and solve 
mygrid=Grid(nx=80,ny=80)
#mygrid.setBCFunc(BC) 
mygrid.setBC(1.0,0.0,0.0,0.0)
stepper='slow'
stepper='numeric'
stepper='blitz'
steppers=[]
steppers.append('slow')
steppers.append('numeric')
steppers.append('blitz')

#for step in steppers:
#    mysol=LaplaceSolver(mygrid,step)
#    print 'completed: ',step

mysolver=LaplaceSolver(mygrid,stepper)
myeps = 1.0e-4
myiter= 400

solveres=mysolver.solve(myiter,myeps)
print 'Completed {} iterations with error = {}'.format(solveres[0],solveres[1])

#Visualitzation of results
x = numpy.arange(mygrid.xmin,mygrid.xmax,mygrid.dx)
y = numpy.arange(mygrid.ymin,mygrid.ymax,mygrid.dy)
#X, Y = numpy.meshgrid(x, y)

plt.figure() 
#plt.imshow(mygrid.u)

#plt.figure()

umax = numpy.max(mygrid.u); umin = numpy.min(mygrid.u)
nclevels = 100              # number of contours
dc = (umax - umin)/nclevels # contour level increment

#Compute solution with a direct solver
n=mygrid.nx*mygrid.ny
nd=5
diagonals=numpy.zeros((5,n))
diagonals[0,:] = -1                       #all elts in first row is set to 1
diagonals[1,:] = -1
diagonals[2,:] =  4
diagonals[3,:] = -1
diagonals[4,:] = -1
As = sc.sparse.spdiags(diagonals, [-mygrid.nx,-1,0,1,mygrid.nx], n, n,format='csc') #sparse matrix instance, numerating in x-direction first

d = numpy.zeros(n)
d[0:mygrid.nx-1] = 1.0
ud = sc.sparse.linalg.spsolve(As,d) #
ud = numpy.reshape(ud,(mygrid.nx,mygrid.ny))
(udmin,udmax) = (numpy.min(ud),numpy.max(ud))
udc=(udmax-udmin)/nclevels

print 'umax =', udmax, umax
print 'umin =', udmin, umin

#CS = plt.contourf(x,y,mygrid.u[0:-1,0:-1],levels=numpy.arange(umin,umax,dc))
plt.figure(1)
#CS = plt.contourf(x,y,mygrid.u[1:,1:],levels=numpy.arange(umin,umax,dc))


#plt.figure(2)
CS2 = plt.contourf(x,y,ud[1:,1:],levels=numpy.arange(udmin,udmax,udc))

#print ud.shape
#print 
plt.show()
plt.close()


du = mygrid.u - ud

print 'maxdu=',numpy.max(du)
print 'mindu=',numpy.min(du)

