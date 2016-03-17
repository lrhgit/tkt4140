import numpy as np
cimport numpy as np

def gauss(np.ndarray[np.float64_t, ndim=2] U, double reltol, double h, double omega):
    cdef Py_ssize_t i, j, it
    cdef double rel_res, du_max, df, dT, R
    
    it=0
    rel_res=1.0
        
    itmax=20
    
    while ((rel_res>reltol) and (it<=itmax)):
        for j in range(1,U.shape[0]-2):
                for i in range(1,U.shape[1]-2):
                    R = (U[j,i-1]+U[j-1,i]+U[j,i+1]+U[j+1,i]-4.0*U[j,i]) + h**2*(U[j,i]**2+1.0)
                    df=4-2*h**2*U[j,i]
                    dT =  omega*R/df
                    U[j,i]+=dT
                    du_max=np.max([np.abs(dT),du_max])
                   
        rel_res=du_max/np.max(np.abs(U))
        it+=1
    if (it==itmax): print 'Terminated after max iterations'    
    return U, rel_res, it

