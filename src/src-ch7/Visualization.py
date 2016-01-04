# chapter7/src-ch7/Visualization.py
import numpy as np
import matplotlib.pylab as plt


# Change some default values to make plots more readable on the screen
LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT

def plot_Surface_yx(Temp, Ttop, xmax, ymax, Nx, Ny, nx, ny):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    # surfaceplot:
    x = np.linspace(0, xmax, Nx + 1)
    y = np.linspace(0, ymax, Ny + 1)
    
    X, Y = np.meshgrid(x, y)
    
    T = np.zeros_like(X)
    
    T[-1,:] = Ttop
    
    for j in range(1, ny+1):
        for i in range(1, nx + 1):
            T[j, i] = Temp[j + (i-1)*ny - 1]
    

    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, T, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    ax.set_zlim(0, Ttop+10)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('T [$^o$C]')

    
    
    nx=4
    xticks=np.linspace(0.0,xmax,nx+1)
    ax.set_xticks(xticks)
    
    ny=6
    yticks=np.linspace(0.0,ymax,ny+1)
    ax.set_yticks(yticks)
    
    nTicks=5
    dT=int(Ttop/nTicks)
    Tticklist=range(0,Ttop+1,dT)
    ax.set_zticks(Tticklist)
    

#    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

def plot_SurfaceNeumann_xy(Temp, Ttop, Tright, xmax, ymax, Nx, Ny, nxTicks=4, nyTicks=4):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    # surfaceplot:
    x = np.linspace(0, xmax, Nx + 1)
    y = np.linspace(0, ymax, Ny + 1)
    
    X, Y = np.meshgrid(x, y)
    
    T = np.zeros_like(X)
    
    T[-1,:] = Ttop
    T[:,-1] = Tright
    k = 1
    for j in range(Ny):
        T[j,:-1] = Temp[Nx*(k-1):Nx*k]
        k+=1

    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, T, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    ax.set_zlim(0, Ttop)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('T [$^o$C]')
    
    xticks=np.linspace(0.0,xmax,nxTicks+1)
    ax.set_xticks(xticks)
    
    yticks=np.linspace(0.0,ymax,nyTicks+1)
    ax.set_yticks(yticks)
    

    
def visualize_setup_xy(Temp, nx, ny, Ttop):
    
    # Generate boards; board contain the setup of the problem, with numbering of T, board2 containg the corresponding values
    Ttop = str(Ttop)
    board = []
    board2 = []
    
    board.append([])
    board2.append([])
    
    for column in range(nx + 2):
        board[0].append(Ttop+'  ')
        board2[0].append(Ttop+'  ')
        
    k = nx*ny 
    for row in range(ny):
        board.append(['0    '])
        board2.append(['0    '])
        k += -nx
        h = 1
        
        for column in range(nx):
            board[row+1].append('T{0}   '.format(k+h)[0:5])
            board2[row+1].append('{0}   '.format(str(Temp[k+h-1]))[0:5])
            h+=1
            
        board[row+1].append('0    ')
        board2[row+1].append('0    ')
        
    board.append([])
    board2.append([])
    
    for column in range(nx + 2):
        board[-1].append('0    ')
        board2[-1].append('0    ')
    
    for i, row in enumerate(board):
        #print board[row]
        board[i].append('    ')
        for column in board2[i]:
            board[i].append(column)
        
    
    print "        Setup XY:"+'     '*(nx+3) +'results:'
    print_board(board)
    
    


def visualize_setup_yx(Temp, nx, ny, Ttop):
    
    # Generate boards; board contain the setup of the problem, with numbering of T, board2 containg the corresponding values
    
    Ttop = str(Ttop)
    board = []
    board2 = []
    
    board.append([])
    board2.append([])
    
    for column in range(nx + 2):
        board[0].append(Ttop+'  ')
        board2[0].append(Ttop+'  ')
    k = ny 
    
    for row in range(ny):
        board.append(['0    '])
        board2.append(['0    '])
        
        h = 0
        
        for column in range(nx):
            board[row+1].append('T{0}   '.format(k+h*ny)[0:5])
            board2[row+1].append('{0}   '.format(str(Temp[k+h*ny-1]))[0:5])
            h+=1
            
        k += -1
        board[row+1].append('0    ')
        board2[row+1].append('0    ')
        
    board.append([])
    board2.append([])
    
    for column in range(nx + 2):
        board[-1].append('0    ')
        board2[-1].append('0    ')
    
    for i, row in enumerate(board):
        
        board[i].append('    ')
        for column in board2[i]:
            board[i].append(column)
    
    print "        Setup YX:"+'     '*(nx+3) +'results:'
    
    print_board(board)
    
def print_board(board):
    
  for row in board:
      
    print " ".join(row)  
    
def convert_xy_yx(Tempxy, nx, ny):
    
    Temp = [0]*len(Tempxy)
    k = 0
    k2 = 0
    for j in range(nx):
        k3 = k2
        for i in range(ny):
            Temp[k] = Tempxy[k3]
            k3 +=nx
            k +=1
        k2+=1
    
    return Temp
    
         

    
