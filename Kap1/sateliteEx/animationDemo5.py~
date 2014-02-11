import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random

t=0
r=3.0
n=0
A=[]

for x in range(10):
    for y in range(10):
        A.append([random.uniform(0,1),random.uniform(0,1)])
A = np.array(A).transpose()

fig = plt.figure()
line, = plt.plot(A[0],A[1], "x", color="blue")

def update():
    for i in range(100):
        A[0], A[1] = r*A[0]*(1-A[0]), r*A[1]*(1-A[1])
        yield A

def draw(data):
    line.set_xdata(data[0])
    line.set_ydata(data[1])
    return line,

ani = animation.FuncAnimation(fig, draw, update, interval=1000, blit=False)

plt.show()
