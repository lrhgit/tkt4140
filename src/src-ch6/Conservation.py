# ../Kap6/advection_schemes.py

import numpy as np
import matplotlib.pylab as plt
# function defining the initial condition
class FluxLimiters:
    
    def __init__(self, dx, dt, x, RHS, limiter, Flux='burger'):
        
        
        self.dx = dx
        self.dt = dt
        self.x = x
        self.RHS = RHS
        self.limiter = limiter
        if Flux =='burger':
            self.F = self.F_burger
        
    def calck_smoothness(self, u):
        r = np.ones_like(u)
        for k in range(len(u[1:-1])):
            j = k + 1
            uj = u[j]
            ujm = u[j-1]
            ujp = u[j+1]
            if abs(uj-ujm)>1e-10 and abs(ujp - uj)>1e-10:
                r[j] = (uj-ujm)/(ujp - uj)
            elif abs(uj-ujm)<1e-10 and abs(ujp - uj)<1e-10:
                r[j] = 1
            else:
                r[j]=0
            
            r = np.where(r<0, 0, r)
#             elif abs(uj-ujm)>1e-10 and abs(ujp - uj)<1e-10:
#                 top = np.sign(uj-ujm)
#                 bottom = np.sign(ujp - uj)
#                 if top == bottom:
#                     r[j] = 0
#                 else:
#                     r[j] = 0
#             elif abs(uj-ujm)<1e-10 and abs(ujp - uj)>1e-10:
#                 r[j] = 0
                
        return r
    
    def calck_smoothness_array(self, u):
        r = np.ones_like(u[1:-1])
        uj = u[1:-1]
        ujp = u[2:]
        ujm = u[0:-2]
        numerator = uj-ujm
        denominator = ujp-uj
        r = np.where(abs(denominator) >1e-10, numerator/denominator, r) # both gradients are big -> normal smoothness criteria
        r = np.where(r<0, 0, r)
        
        r = np.append(1, np.append(r, 1))
      
        return r
        
    
    def calck_phi(self, r):
        method = self.limiter
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
    
    def wavespeed(self, uj, ujp, tol=1e-10):
        """method that approximates the numerical wavespeed at j+0.5:
                        { F(u{j+1})-F(u{j})/(u{j+1}-u{j}) if u{j+1} != u{1}
            a{j+0.5} =  {
                        {F'(u{j}                          if u{j+1} == u{1}
            
            the method assures conservation
            Args:
                uj(array): an array containg the solution of u{j}
                ujp(array): an array containg the solution of u{j+1}
            Returns:
                A(array): the numerical wavespeed at j+0.5
        """
        F = self.F
        A = np.ones_like(uj)
        for n in range(len(uj)):
            if abs(ujp[n]-uj[n])<tol:
                A[n]=uj[n]
            else:
                A[n]=(F(ujp[n])-F(uj[n]))/(ujp[n]-uj[n])
    
        return A
    
    def F_LW(self, uj, ujp, a, phi):
        """method that calculates the general Flux term F_LW for the Lax-Wendroff Flux
            F(j+0.5,n+0.5) = F(uj) + 0.5*a*(1-a*dt/dx)*(u{j+1}-u{j})
            where a is an aproximation of the numerical wavespeed at j+0.5
        
            Args:
                u(array): an array containg the solution of u
            Returns:
                Flux(array): The Lax-Wendroff Flux .
        """
        dt = self.dt
        dx = self.dx
        F = self.F
        
        Flux = F(uj) + phi*0.5*a*(1-a*dt/dx)*(ujp-uj)
        
        return Flux
    
    def F_burger(self, u):
        
        return 0.5*u**2
    
    def solve(self, u, t):
        #global phi
        dt = self.dt
        dx = self.dx
        x = self.x
        RHS = self.RHS
        
        r = self.calck_smoothness_array(u)
        phi = self.calck_phi(r)
        
        phi_p = phi[1:-1]
        phi_m = phi[:-2]
        
        ujm = u[:-2].copy() #u(j-1)
        uj = u[1:-1].copy() #u(j)
        ujp = u[2:].copy() #u(j+1)
        ajp = self.wavespeed(uj, ujp)
        ajm = self.wavespeed(ujm, uj)
        
        vjp = ajp*dt/dx
        vjm = ajm*dt/dx
    
        Qjm = RHS(t-dt, x[:-2])
        Qj = RHS(t-dt, x[1:-1])
        Qjp = RHS(t-dt, x[2:])
        Qjpnp = RHS(t, x[1:-1])
        
        Source = 0.5*dt*(Qj+Qjpnp)-0.25*dt*(vjp*(Qjp+Qj)-vjm*(Qj+Qjm))
        
        u[1:-1] = uj -(dt/dx)*(self.F_LW(uj, ujp, ajp, phi_p) - self.F_LW(ujm, uj, ajm, phi_m)) + Source
        return u[1:-1]
    
    def name(self):
        
        limiter = self.limiter
        
        if limiter =='lax_wendroff':
            name = 'Lax-W'
            
        elif limiter =='upwind':
            name = 'upwind'
            
        elif limiter =='minmod':
            name = 'minmod'
            
        elif limiter =='van_leer':
            name = 'van-Leer'
            
        elif limiter =='superbee':
            name = 'Superbee'
        else:
            name = limiter

        return name

class Classical:
    
    def __init__(self, dx, dt, x, solver, RHS, Flux='burger'):
        
        
        self.dx = dx
        self.dt = dt
        self.x = x
        self.solver_name = solver
        self.solve = eval('self.'+ solver)
        self.RHS = RHS
        if Flux =='burger':
            self.F = self.F_burger
        
        
    
    
    def F_burger(self, u):
        
        
        return 0.5*u**2
    
    def macCormack(self, u, t):
        """method that solves u(n+1), for the scalar conservation equation with source term:
            du/dt + dF/dx = RHS,
            where F = 0.5u^2 for the burger equation
            with use of the MacCormack scheme
            
            Args:
                u(array): an array containg the previous solution of u, u(n). (RHS)
                t(float): an array 
            Returns:
                u[1:-1](array): the solution of the interior nodes for the next timestep, u(n+1).
        """
        dt = self.dt
        dx = self.dx
        x = self.x
        RHS = self.RHS
        F = self.F
        
        up = u.copy()
        up[:-1] = u[:-1] - (dt/dx)*(F(u[1:]) - F(u[:-1])) + dt*RHS(t-0.5*dt, x[:-1])
        u[1:] = .5*(u[1:] + up[1:] -  (dt/dx)*(F(up[1:]) - F(up[:-1])) + dt*RHS(t-0.5*dt, x[1:])) 
        
        return u[1:-1]
    
    def Lax_W_Two_Step(self, u, t):
        """method that solves u(n+1), for the scalar conservation equation with source term:
            du/dt + dF/dx = RHS,
            where F = 0.5u^2 for the burger equation
            with use of the Two-step Lax-Wendroff scheme
            
            Args:
                u(array): an array containg the previous solution of u, u(n).
                t(float): time at t(n+1) 
            Returns:
                u[1:-1](array): the solution of the interior nodes for the next timestep, u(n+1).
        """
        dt = self.dt
        dx = self.dx
        x = self.x
        RHS = self.RHS
        F = self.F
        
        ujm = u[:-2].copy() #u(j-1)
        uj = u[1:-1].copy() #u(j)
        ujp = u[2:].copy() #u(j+1)
        up_m = 0.5*(ujm + uj) - 0.5*(dt/dx)*(F(uj)-F(ujm)) + 0.5*dt*RHS(t-0.5*dt, x[1:-1] - 0.5*dx) #u(n+0.5dt,j-0.5dx)
        up_p = 0.5*(uj + ujp) - 0.5*(dt/dx)*(F(ujp)-F(uj)) + 0.5*dt*RHS(t-0.5*dt, x[1:-1] + 0.5*dx)#u(n+0.5dt,j+0.5dx)
        
        u[1:-1] = uj -(dt/dx)*(F(up_p) - F(up_m)) + dt*RHS(t-0.5*dt, x[1:-1])
        return u[1:-1]
    

    def lax_friedrich_Flux(self, u, t):
        """method that solves u(n+1), for the scalar conservation equation with source term:
            du/dt + dF/dx = RHS,
            where F = 0.5u^2 for the burger equation
            with use of the lax-friedrich scheme
            
            Args:
                u(array): an array containg the previous solution of u, u(n). (RHS)
                t(float): an array 
            Returns:
                u[1:-1](array): the solution of the interior nodes for the next timestep, u(n+1).
        """
        dt = self.dt
        dx = self.dx
        x = self.x
        RHS = self.RHS
        F = self.F
        
        u[1:-1] = (u[:-2] +u[2:])/2.0 -  dt*(F(u[2:])-F(u[:-2]))/(2.0*dx) + dt*(RHS(t, x[:-2]) + RHS(t, x[2:]))/2.0
        return u[1:-1]
    
    def ftbs(self, u, t):
        """method that solves u(n+1), for the scalar conservation equation with source term:
            du/dt + dF/dx = RHS,
            where F = 0.5u^2 for the burger equation
            with use of the forward in time backward in space (upwind) scheme
            
            Args:
                u(array): an array containg the previous solution of u, u(n). (RHS)
                t(float): an array 
            Returns:
                u[1:-1](array): the solution of the interior nodes for the next timestep, u(n+1).
        """
        dt = self.dt
        dx = self.dx
        x = self.x
        RHS = self.RHS
        F = self.F
        
        u[1:-1] = u[1:-1] -  (dt/dx)*(F(u[1:-1])-F(u[:-2])) + dt*RHS(t-0.5*dt, x[1:-1])
        
        return u[1:-1]
        
    def name(self):
        
        return self.solver_name
