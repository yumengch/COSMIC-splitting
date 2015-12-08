def streamfunction(psi,dx,dy,dt,initialProfile):
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)
    #define the velocity field
    u[:-1,:] = -(psi[1:,:]-psi[:-1,:])/dy
    v[:,:-1] =  (psi[:,1:]-psi[:,:-1])/dx
    
    #calculate the Courant number
    cx = u*dt/dx
    cy = v*dt/dy
    return u,v,cx,cy