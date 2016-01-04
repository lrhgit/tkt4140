# ../Kap6/advection_schemes.py

import numpy as np
import matplotlib.pylab as plt
# function defining the initial condition
class FluxLimiters:
    
    def __init__(self, a, dx, dt, c, limiter):
        
        self.c = c
        self.dx = dx
        self.dt = dt  
        self.a = a
        self.limiter = limiter
        
    def calck_smoothness(self, u):
        thetha = np.ones_like(u)
        for k in range(len(u[1:-1])):
            j = k + 1
            uj = u[j]
            ujm = u[j-1]
            ujp = u[j+1]
            if abs(uj-ujm)>1e-10 and abs(ujp - uj)>1e-10:
                thetha[j] = (uj-ujm)/(ujp - uj)
            elif abs(uj-ujm)<1e-10 and abs(ujp - uj)<1e-10:
                thetha[j] = 1
            elif abs(uj-ujm)>1e-10 and abs(ujp - uj)<1e-10:
                top = np.sign(uj-ujm)
                bottom = np.sign(ujp - uj)
                if top == bottom:
                    thetha[j] = 1
                else:
                    thetha[j] = -1
            elif abs(uj-ujm)<1e-10 and abs(ujp - uj)>1e-10:
                thetha[j] = 0
                
        return thetha
    
    def calck_smoothness_array(self, u):
        thetha = np.zeros_like(u)
        uj = u[1:-1]
        ujp = u[2:]
        ujm = u[0:-2]
        top = uj-ujm
        bottom = ujp-uj
        top_sign = np.sign(top)
        bottom_sign = np.sign(bottom)
        thetha = np.where(abs(top)>1e-10 & abs(bottom) >1e-10, top/bottom, thetha) # both gradients are big -> normal smoothness criteria
        thetha = np.where(abs(top)>1e-10 & abs(bottom) <1e-10, top_sign*bottom_sign*1, thetha) # division by zero ->big gradient difference
        thetha = np.where(abs(top)<1e-10 & abs(bottom) >1e-10, 0, thetha) # zero gradient on top-> nonzero bottom -> zero
                
        return thetha
    
    def calck_phi(self, thethaList):
        method = self.limiter
        phi = np.ones_like(thethaList)
        if method=='minmod':
            phi = np.where((thethaList>0) & (abs(thethaList)>=1), 1, phi)
            phi = np.where((thethaList>0) & (abs(thethaList)<1), thethaList, phi)
            phi = np.where((thethaList<=0), 0, phi)
        elif method=='superbee':
            phi = np.where((thethaList<=0), 0, phi)
            phi = np.where((thethaList>0) & (thethaList<=0.5), 2*thethaList, phi)
            phi = np.where((thethaList>0.5) & (thethaList<=1), 1, phi)
            phi = np.where((thethaList>1)  & (thethaList<=2), thethaList, phi)
            phi = np.where((thethaList>2), 2, phi)
        elif method=='van_leer':
            phi = (thethaList + abs(thethaList))/(1 + abs(thethaList))
        elif method=='lax_wendroff':
            phi[:] = 1
        elif method=='upwind':
            
            phi[:] = 0 

        
        return phi
    
    def Flux(self, ua, ub, phi):
        #print "phi used in flux:", phi
        a = self.a
        c = self.c
        Flux = a*(ua + phi*0.5*(1-c)*(ub - ua))
        return Flux
    
    def solve(self, u):
        #global phi
        dt = self.dt
        dx = self.dx 
        thetha = self.calck_smoothness(u)
        phi = self.calck_phi(thetha)
        
        if self.limiter =='supjjerbee':
            fig , ax = plt.subplots(3, 1, squeeze=False)
            ax[0][0].plot(u)
            ax[1][0].plot(thetha)
            ax[2][0].plot(phi)
            plt.show()
            raw_input('hey')
        
        ua_pluss = u[1:-1]
        ub_pluss = u[2:]
        phi_pluss = phi[1:-1]
        
        ua_minus = u[:-2]
        ub_minus = u[1:-1]
        phi_minus = phi[:-2]
        
        u[1:-1] = u[1:-1] -  (dt/dx)*(self.Flux(ua_pluss, ub_pluss, phi_pluss) - self.Flux(ua_minus, ub_minus, phi_minus))
        return u[1:-1]
    
    def name(self):
        
        limiter = self.limiter
        
        if limiter =='lax_wendroff':
            name = 'Lax-W'
            
        elif limiter =='upwind':
            name = 'upwind'
            
        elif limiter =='minmod':
            name = 'Lax-W-MM'
            
        elif limiter =='van_leer':
            name = 'Lax-W-VanL'
            
        elif limiter =='superbee':
            name = 'Lax-W-sbee'
        else:
            name = limiter

        return name

class FluxLimiters2:
    
    def __init__(self, a, dx, dt, c, limiter):
        
        self.c = c
        self.dx = dx
        self.dt = dt  
        self.a = a
        self.limiter = limiter
        
    def calck_smoothness(self, u):
        thetha = np.zeros_like(u)
        for k in range(len(u[1:-1])):
            j = k + 1
            uj = u[j]
            ujm = u[j-1]
            ujp = u[j+1]
            if abs(uj-ujm)>1e-10 and abs(ujp - uj)>1e-10:
                thetha[j] = (uj-ujm)/(ujp - uj)
            elif abs(uj-ujm)<1e-10 and abs(ujp - uj)<1e-10:
                thetha[j] = 1
            elif abs(uj-ujm)>1e-10 and abs(ujp - uj)<1e-10:
                top = np.sign(uj-ujm)
                bottom = np.sign(ujp - uj)
                if top == bottom:
                    thetha[j] = 100
                else:
                    thetha[j] = -100
            elif abs(uj-ujm)<1e-10 and abs(ujp - uj)>1e-10:
                thetha[j] = 0
                
        return thetha
    
    def calck_smoothness_array(self, u):
        thetha = np.zeros_like(u)
        uj = u[1:-1]
        ujp = u[2:]
        ujm = u[0:-2]
        top = uj-ujm
        bottom = ujp-uj
        top_sign = np.sign(top)
        bottom_sign = np.sign(bottom)
        thetha = np.where(abs(top)>1e-10 & abs(bottom) >1e-10, top/bottom, thetha) # both gradients are big -> normal smoothness criteria
        thetha = np.where(abs(top)>1e-10 & abs(bottom) <1e-10, top_sign*bottom_sign*1, thetha) # division by zero ->big gradient difference
        thetha = np.where(abs(top)<1e-10 & abs(bottom) >1e-10, 0, thetha) # zero gradient on top-> nonzero bottom -> zero
                
        return thetha
    
    def calck_phi(self, thethaList):
        method = self.limiter
        phi = np.ones_like(thethaList)
        if method=='minmod':
            phi = np.where((thethaList>0) & (abs(thethaList)>1), 1, phi)
            phi = np.where((thethaList>0) & (abs(thethaList)<1), thethaList, phi)
            phi = np.where((thethaList<0), 0, phi)
        elif method=='superbee':
            phi = np.where((thethaList<0), 0, phi)
            phi = np.where((thethaList>0) & (thethaList<0.5), 2*thethaList, phi)
            phi = np.where((thethaList>0.5) & (thethaList<1), 1, phi)
            phi = np.where((thethaList>1)  & (thethaList<2), thethaList, phi)
            phi = np.where((thethaList>2), 2, phi)
        elif method=='van_leer':
            phi = (thethaList + abs(thethaList))/(1 + abs(thethaList))
        elif method=='lax_wendroff':
            phi[:] = 1
        elif method=='upwind':
            
            phi[:] = 0 

        
        return phi
    
    def Flux(self, ua, ub, phi):
        #print "phi used in flux:", phi
        F_low = self.Flux_upwind(ua, ub)
        F_high = self.Flux_lax_wendroff(ua, ub)
        
        F = F_low + phi*(F_high - F_low)
        return F
    
    def Flux_lax_wendroff(self, ua, ub):
        a = self.a
        c = self.c
        
        return a*(ua + 0.5*(1-c)*(ub - ua))
    
    def Flux_upwind(self, ua, ub):
        a = self.a
        
        return a*ua
    
    def solve(self, u):
        #global phi
        dt = self.dt
        dx = self.dx 
        thetha = self.calck_smoothness(u)
        phi = self.calck_phi(thetha)
        
        ua_pluss = u[1:-1]
        ub_pluss = u[2:]
        phi_pluss = phi[1:-1]
        
        ua_minus = u[:-2]
        ub_minus = u[1:-1]
        phi_minus = phi[:-2]
        
        u[1:-1] = u[1:-1] -  (dt/dx)*(self.Flux(ua_pluss, ub_pluss, phi_pluss) - self.Flux(ua_minus, ub_minus, phi_minus))
        return u[1:-1]
    
    def name(self):
        
        return self.limiter + '2'
