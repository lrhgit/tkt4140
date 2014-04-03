from numpy import*
from Tkinter import *
import time

dt = 1
dx = 1
D = 1
alpha = D*dt/dx**2
nodes = 500
source_node = nodes/2
u = zeros(nodes)
u[source_node] = 1
x = arange(0, nodes, 1)
a = ones(nodes)*(-alpha)
b = ones(nodes)*(1+2.0*alpha)
c = ones(nodes)*(-alpha)
a[0]=0; a[nodes-1] = 0
b[0]=1; b[nodes-1]=1
c[0]=0; c[nodes-1]=0
cprime = zeros(nodes)
dprime = zeros(nodes)

##--Tkinter stuff--
t_start = time.time()#for referencing
height = 300.0
width = 400.0
center = height//2
x_factor = width/(nodes-1)
y_factor = 200
##------------------
    
def paint(canvas, parent):
    #see NMM equations (9.25) and (9.26)
    cprime[0] = c[0]/b[0]
    dprime[0] = u[0]/b[0]
    for i in range(1, nodes, 1):
        cprime[i] = c[i]/(b[i]-a[i]*cprime[i-1])
        dprime[i] = (u[i]-a[i]*dprime[i-1])/(b[i]-a[i]*cprime[i])
    #now, the backward pass:
    u[nodes-1] = dprime[nodes-1]
    for i in range (nodes-2, 0, -1):
        u[i] = dprime[i] - cprime[i] * u[i+1]        
    ##--Tkinter stuff--
    xy = []
    for i in range(0, nodes):
        xy.append((int)(i*x_factor))
        xy.append((int)(u[i]*y_factor)+center)
    #time.sleep(0.001)    
    canv.coords("curve", *xy)
    parent.after_idle(paint,parent,canvas)
    ##------------------    
#--Tkinter stuff:--
root = Tk()
root.title("Animated solution to 1D diffusion equation")
root.bind('q','exit')
canv = Canvas(width=width, height=height, bg='white')
canv.pack()
canv.create_line(tag = "curve", *zeros(2*width), fill='blue')    
root.after(200,paint,root,canv)
root.mainloop()    
##------------------
    
    
    