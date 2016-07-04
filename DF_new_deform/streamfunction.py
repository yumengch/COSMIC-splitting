def streamfunction(psi,dx,dy,dt,J):
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)
    #define the velocity field
    u[:-1,:] = -(psi[1:,:]-psi[:-1,:])/dy
    v[:,:-1] =  (psi[:,1:]-psi[:,:-1])/dx
    v[:,-1],v[-1,:] = v[:,0],v[0,:]
    u[-1,:],u[:,-1] = u[0,:],u[:,0]
    #calculate the Courant number
    cx = J*u*dt/dx
    cy = J*v*dt/dy
    return u,v,cx,cy