import numpy as np
cimport numpy as np

def gauss(np.ndarray[np.float64_t, ndim=2] U, double reltol, double h, double omega):
    cdef Py_ssize_t i, j, it
    cdef double rel_res, dU_max, df, dU, R
    
    
    
    it=0
    rel_res=1.0
        
    itmax=100
    
    while ((rel_res>reltol) and (it<=itmax)):
        dU_max=0.0
        for j in range(1,U.shape[0]-2):
                for i in range(1,U.shape[1]-2):
                    R = (U[j,i-1]+U[j-1,i]+U[j,i+1]+U[j+1,i]-4.0*U[j,i]) + h**2*(U[j,i]**2+1.0)
                    df=4.0-2*h**2*U[j,i]
                    dU =  omega*R/df
                    U[j,i]+=dU
                    dU_max=np.max([np.abs(dU),dU_max])
                   
        rel_res=dU_max/np.max(np.abs(U[:,:]))
#         print 'rel_res', rel_res
        it+=1
    if (it>=itmax): print 'Terminated after max iterations'    
    return U, rel_res, it

