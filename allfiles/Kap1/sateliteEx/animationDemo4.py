import numpy as np
import matplotlib.pyplot as plt
plt.ion()

import random
t=0
r=3.0
n=0
A=[]
for x in range(10):
    for y in range(10):
        A.append([random.uniform(0,1),random.uniform(0,1)])

for m in range(len(A)):
    plt.scatter(A[m][0],A[m][1], s=250)
    plt.draw()
plt.pause(1)

while n<=100:
    for m in range(len(A)):
        A[m][0]=r*A[m][0]*(1-A[m][0])
        A[m][1]=r*A[m][1]*(1-A[m][1])
    for m in range(len(A)):
        plt.scatter(A[m][0],A[m][1], s=250)
    plt.draw()
    plt.pause(1)
    plt.clf()

