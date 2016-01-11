# src/src-ch5/sphere.py


#import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from sphereAnimation import animateSphere
from math import pi, sin, exp

#change some default values to make plots more readable 
LNWDT=3; FNT=11
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT

def analytic(r, t):
    fac1 = 0.5*pi/b
    sign = 1
    criteria = 1
    epsi = 1e-5
    s = 0
    m = -1
    f1, f2 = b*1e-5, b*(1-1e-5)
    if (r >= f1) and (r <= f2):
        while (criteria > epsi):
            m = m + 2
            bm = m*fac1
            bmr = r*bm
            temp = exp(-alfa*bm**2*t)/m**2
            term = sign*sin(bmr)*temp
            s = s + term
            sign = -sign
            criteria = abs(temp/s)
        Tvalue = Tv + 8*b*(Tk - Tv)*s/(pi**2*r)
    
    if r < f1:
        while (criteria > epsi):
            m = m + 2
            bm = m*fac1
            term = sign*exp(-alfa*bm**2*t)/m
            s = s + term
            sign = -sign
            criteria = abs(term/s)
        Tvalue = Tv + 4*(Tk - Tv)*s/pi
    
    if r > f2:
        while (criteria > epsi):
            m = m + 2;  
            bm = m*fac1;  
            term = exp(-alfa*bm**2*t)/m**2
            s = s + term
            criteria = abs(term/s)
            
        Tvalue = Tv + 8*(Tk - Tv)*s/pi**2
    
    return Tvalue


"""script calculationg the temperature variation in a sphere with
    radius b, cooling in water.
    Initial temperature = Tk
    Temperature in water = Tv = constant
    Calculations are done with FTCS scheme, with r - steplength dr
    and numerical Fourier number D given. 
    Using 2. order backward difference in BC for r = b = 5 cm
    rb = 1: Using separate equation
        for r = 0. Stable for D < 1/3
    rb = 2: Using 2. order forward differance
        for r = 0. stable for D < 1/3
    The analytical solution is calculated by method kanalyt. The condition
    H*b/k = 1 has been used, where h is the heat transfer number and k the conductivity.
"""

b = 5. # radius of sphere [cm]
N = 50 # number of parts
r_array = np.linspace(0, b, N+1)
alfa = 0.04 # diffusivity [cm^2/s
dr = b/N
D = 0.4 # numerical fourrier number
dtau = D*dr**2/alfa # time step [s]
tauend = 600 # stop time [s]
ntot = int(tauend/dtau) + 1
K = 0.1 # conductivity [W/cm/C]
H = 0.02 # heat transfer number [W/cm^2/C]
dra = dr*H/K
ta = np.zeros(N + 1) # initialize analytic solution
Tk = 300 # start temperature in the sphere
Tv = 20 # Temperature in the water

Told = Tk*np.ones(N + 1) # initial values
Tnew = Tk*np.ones(N + 1)
Tmatrix = np.zeros((ntot, N + 1)) # initialize solution matrix
Tmatrix[0,:] = Told
Tmatrix_analytic = np.zeros((ntot, N + 1))
Tmatrix_analytic[0,:] = Told
temp_analytic = np.zeros_like(r_array)
rb = 2


for k in range(1,ntot):
    tau = k*dtau
    for j in range(1,N):

        Tnew[j] = (1 - 2*D)*Told[j] + D*((1 - 1./j)*Told[j - 1] + (1 + 1./j)*Told[j + 1]) #FTCS
        
    if rb == 2:
        Tnew[0] = (4*Tnew[1]-Tnew[2])/3.
    elif rb == 1:
        Tnew[0] = (1 - 6*D)*Told[0] + 6*D*Told[1]

    Tnew[N] = (4*Tnew[N - 1] - Tnew[N - 2] + 2*dra*Tv)/(3. + 2*dra)
    
    Told = Tnew.copy()

    for i , r in enumerate(r_array):
        temp_analytic[i]=analytic(r, tau)
        
    Tmatrix_analytic[k:] = temp_analytic
    Tmatrix[k:]= Tnew
    

TsnapshotList = [0, 10, 100, 300, 600]
legendList = []

for Tsnapshot in TsnapshotList:
    Temp = Tmatrix[Tsnapshot,:]
    time = Tsnapshot*dtau
    plot(Temp)
    legendList.append(str(int(time)) + ' s')

for Tsnapshot in TsnapshotList:
    Temp = Tmatrix_analytic[Tsnapshot,:]
    time = Tsnapshot*dtau
    plot( Temp, 'k--')
    #legendList.append(str(int(time)) + ' s')   
 
legend(legendList,loc='best',frameon=False)
title('Temperature in Sphere')
#axis([0, 5, 150, 310])
xlabel('r [cm]')
ylabel('T [$^o C$]')
#savefig('sphere.pdf')

        
showanimation = True

if showanimation:
    animateSphere(r_array, Tmatrix, Tmatrix_analytic, ntot, dtau)
show()


