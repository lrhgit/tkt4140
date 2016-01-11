from sympy import symbols, Rational, simplify
import numpy as np
import matplotlib.pylab as plt
# a, uj, ujp, ujm, c, dt, dx = symbols('a uj ujp ujm c dt dx')
# 
# Fp = (uj + Rational(1,2)*(1-c)*(ujp-uj))
# Fm = (ujm + Rational(1,2)*(1-c)*(uj-ujm))
# Lax1 = uj - c*Rational(1,2)*(ujp - ujm) + c**2*Rational(1,2)*(ujp - 2*uj + ujm)
# Lax2 = c*Rational(1,2)*(1+c)*ujm + (1-c**2)*uj - c*Rational(1,2)*(1-c)*ujp
# F = c*(Fp-Fm)
# Lax3 = uj - F
# 
# print "F: ", F
# print  simplify(F)
# print simplify(Lax1)
# print simplify(Lax2)
# print simplify(Lax3)
# print simplify(Lax1-Lax2)
# print simplify(Lax1-Lax3)
# print simplify(Lax2-Lax3)

def calck_smoothness_array(u):
    r = np.ones_like(u[1:-1])
    uj = u[1:-1]
    ujp = u[2:]
    ujm = u[0:-2]
    numerator = uj-ujm
    denominator = ujp-uj
    r = np.where(abs(denominator) >1e-10, numerator/denominator, r) # both gradients are big -> normal smoothness criteria
    #r = np.where(r<0, 0, r)
    
    r = np.append(1, np.append(r, 1))
  
    return r

def calck_phi(r):
    method = limiter
    phi = np.ones_like(r)
    if method=='minmod':
        phi = np.where((r>0) & (abs(r)>=1), 1, phi)
        phi = np.where((r>0) & (abs(r)<1), r, phi)
        phi = np.where((r<=0), 0, phi)
    elif method=='superbee':
        phi = np.where((r<=0), 0, phi)
        phi = np.where((r>0) & (r<=0.5), 2*r, phi)
        phi = np.where((r>0.5) & (r<=1), 1, phi)
        phi = np.where((r>1)  & (r<=2), r, phi)
        phi = np.where((r>2), 2, phi)

    elif method=='Fredrik':
        phi = np.where((r<=0), 0, phi)
        phi = np.where((r>0) & (r<=0.5), 2*r, phi)
        phi = np.where((r>0.5), 1, phi)
    elif method=='van_leer':
        phi = (r + abs(r))/(1 + abs(r))
    elif method=='lax_wendroff':
        phi[:] = 1
    elif method=='upwind':
        
        phi[:] = 0 

    
    return phi

def superbee(thethaList):
    phi = np.ones_like(thethaList)

    phi = np.where((thethaList<=0), 0, phi)
    phi = np.where((thethaList>0) & (thethaList<=0.5), 2*thethaList, phi)
    phi = np.where((thethaList>0.5) & (thethaList<=1), 1, phi)
    phi = np.where((thethaList>1)  & (thethaList<=2), thethaList, phi)
    phi = np.where((thethaList>=2), 2, phi)
    
    return phi
def minmod(thethaList):
    phi = np.ones_like(thethaList)
    phi = np.where((thethaList>0) & (abs(thethaList)>=1), 1, phi)
    phi = np.where((thethaList>0) & (abs(thethaList)<1), thethaList, phi)
    phi = np.where((thethaList<=0), 0, phi)
    return phi

def VanLeer(thethaList):
    
    phi = (thethaList + abs(thethaList))/(1 + abs(thethaList))
    return phi

def warming(thethaList):
    return thethaList

def TDV(thethaList):
    warm = np.ones_like(thethaList)*2
    
    warm = np.where(thethaList<1, 2*thethaList, warm)
    warm = np.where(thethaList>=1, 2, warm)
    
    return warm

def Lax(thethaList):
    
    return np.ones_like(thethaList)

def init_step(x):
    """Assigning a value of 1.0 for values less than 0.1"""
    f = np.zeros_like(x)
    f[np.where(x <= 0.1)] = 1.0
    return f

def init_sine2(x):
    """A smooth sin^2 function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = 0.45
    x_right = 0.55
    f = np.where((x>x_left) & (x<x_right), np.sin(np.pi*(x-x_left)/(x_right-x_left))**4,f) 
    return f

def init_box(x):
    """A smooth sin^2 function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = dx
    x_right = 0.2
    f = np.where((x>x_left) & (x<x_right), 1,f) 
    return f

def init_box_sine(x):
    """A smooth sin^2 function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = dx
    x_right = 0.2
    f = np.where((x>x_left) & (x<x_right), 1,f)
    x_left2 = 0.4
    x_right2 = 0.65
    f = np.where((x>x_left2) & (x<x_right2), np.sin(np.pi*(x-x_left2)/(x_right2-x_left2))**2,f) 
    return f
def init_box_sine_sine2(x):
    """A smooth sin^2 function between x_left and x_right"""
    f = np.zeros_like(x)
    x_left = dx
    x_right = 0.2
    f = np.where((x>x_left) & (x<x_right), 1,f)
    x_left1 = 0.3
    x_right1 = 0.4
    f = np.where((x>x_left1) & (x<x_right1), np.sin(np.pi*(x-x_left1)/(x_right1-x_left1))**4,f) 
    x_left2 = 0.5
    x_right2 = 0.75
    f = np.where((x>x_left2) & (x<x_right2), np.sin(np.pi*(x-x_left2)/(x_right2-x_left2))**2,f) 
    return f

LNWDT=2.5; FNT=12
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT


limiter = 'superbee'
x = np.linspace(0,1,101)
dx = x[1]-x[0]

u = init_box_sine(x)

r = calck_smoothness_array(u)
phi_sup = superbee(r)
phi_van = VanLeer(r)
phi_min = minmod(r)
phi_sup_class = calck_phi(r)
limiter = 'minmod'
phi_min_class = calck_phi(r)
limiter = 'van_leer'
phi_van_class = calck_phi(r)
#phi_sup = superbee(r)

fig , ax = plt.subplots(3, 1, squeeze=False)
ax[0][0].plot(x, u, 'b')
ax[0][0].plot(x, u, 'bo')
ax[1][0].plot(x, r, 'b')
ax[1][0].plot(x, r, 'bo')

ax[2][0].plot(x, phi_sup, 'b')
#ax[2][0].plot(x, phi_sup_class, 'y--')
ax[2][0].plot(x, phi_van, 'g')
#ax[2][0].plot(x, phi_van_class, 'c--')
ax[2][0].plot(x, phi_min, 'r')
#ax[2][0].plot(x, phi_min_class, 'm--')
ax[2][0].plot(x, phi_sup, 'bo')
ax[2][0].plot(x, phi_van, 'go')
ax[2][0].plot(x, phi_min, 'ro')
legends= ['superbee','van-Leer', 'minmod',]#legends= ['sup','sup2', 'van', 'van2', 'min', 'min2']
ax[0][0].set_ylabel('u')
ax[0][0].set_xlabel('x')
ax[1][0].set_ylabel('r')
ax[1][0].set_xlabel('u')
ax[2][0].set_ylabel('$\phi$')
ax[2][0].set_xlabel('r')
ax[2][0].legend(legends)
fig.tight_layout()
plt.savefig('../../figs/limiter_example.png') #, transparent=True
plt.show()



# print VanLeer(2.8)
# 
# r = np.linspace(0, 3, 301)
# sup = superbee(r)
# mm = minmod(r)
# lax = Lax(r)
# van = VanLeer(r)
# warm = warming(r)
# warm2 = TDV(r)
# 
# 
# plt.figure()
# plt.plot(r, warm, 'k')
# plt.plot(r,lax,'k')
# plt.fill_between(r, 0, warm2, color='grey')
# plt.annotate('Lax-Wendroff, $\phi=1$', xy=(2.5, 1), xycoords='data',
#                 xytext=(-70, 25), textcoords='offset points',
#                 arrowprops=dict(arrowstyle="->",
#                                 connectionstyle="arc3,rad=0.1"),
#                 )
# 
# plt.annotate('Warming and Beam, $\phi=r$', xy=(2.5, 2.5), xycoords='data',
#                 xytext=(-150, 30), textcoords='offset points',
#                 arrowprops=dict(arrowstyle="->",
#                                 connectionstyle="arc3,rad=0.1"),
#                 )
# 
# 
# 
# plt.xlabel('$r$')
# plt.ylabel('$\phi$')
# plt.xticks([0,1,2,3])
# plt.yticks([1,2,3])
# plt.savefig('../figs/TVD.png')
# plt.figure()
# 
# plt.fill_between(r, mm, sup, color='grey')
# plt.xlabel('$r$')
# plt.ylabel('$\phi$')
# plt.xticks([0,1,2,3])
# plt.yticks([1,2,3])
# plt.savefig('../figs/2nd_TVD.png')
# plt.figure()
# 
# plt.fill_between(r, mm, sup, color='grey')
# plt.xlabel('$r$')
# plt.ylabel('$\phi$')
# plt.xticks([0,1,2,3])
# plt.yticks([1,2,3])
# #plt.savefig('../figs/2nd_TVD.png')
# # plt.figure()
# plt.plot(r,sup, 'b')
# plt.plot(r,mm, 'r')
# plt.plot(r,van, 'y')
# plt.annotate('Superbee', xy=(1.5, 1.5), xycoords='data',
#                 xytext=(-60, 30), textcoords='offset points',
#                 arrowprops=dict(arrowstyle="->",
#                                 connectionstyle="arc3,rad=0.1"),
#                 )
# plt.annotate('Van-Leer', xy=(2.8, 1.4736), xycoords='data',
#                 xytext=(-60, 30), textcoords='offset points',
#                 arrowprops=dict(arrowstyle="->",
#                                 connectionstyle="arc3,rad=-0.1"),
#                 )
# plt.annotate('Min-Mod', xy=(2, 1.), xycoords='data',
#                 xytext=(-20, -40), textcoords='offset points',
#                 arrowprops=dict(arrowstyle="->",
#                                 connectionstyle="arc3,rad=-0.1"),
#                 )
# plt.xticks([0,1,2,3])
# plt.yticks([1,2])
# #plt.savefig('../figs/2nd_TVD_limiters2.png')
# 
# plt.show()