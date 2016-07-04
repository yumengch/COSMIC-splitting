def conservative(phiOld,cx,cy,nt,epsilon,dx,dy,J):
    #Calculate the number of grid points
    nx = len(phiOld[0,:]) - 1
    ny = len(phiOld[:,0]) - 1
    #the phi values at the boundary of each cell
    phi_mid = np.zeros_like(phiOld)
    #advective-form one-dimensional phi update
    XC = np.zeros_like(phiOld)
    phi_CX = np.zeros_like(phiOld)
    #conservative-form one-dimensional flux difference of each cell
    YC_CX = np.zeros_like(phiOld)
    mass = np.zeros((ny,nx))
    #the updated phi values
    phi = np.zeros_like(phiOld)

    for t in xrange(int(nt)):
        #Calculate the flux on each dimension
        for j in xrange(ny):
            phi_mid[j,:], mass[j,:]= PPM(phiOld[j,:]/J[j,:],J[j,:]*cx[j,:],nx,epsilon,dx)
            OUT = flux(nx,dx,J[j,:]*cx[j,:],phi_mid[j,:-1],mass[j,:])
            XC[j,:-1] = conf(nx,J[j,:]*cx[j,:],OUT)
            XC[j,-1] = XC[j,0]
        XC[-1,:] = XC[0,:]

        phi_CX= phiOld+ XC

        #Calculate the flux on each dimension using and phi_AX
        for i in xrange(nx):
            phi_mid[:,i],mass[:,i] = PPM(phi_CX[:,i],J[:,i]*cy[:,i],ny,epsilon,dy)
            OUT = flux(ny,dy,J[:,i]*cy[:,i],phi_mid[:-1,i],mass[:,i])
            YC_CX[:-1,i] = conf(ny,J[:,i]*cy[:,i],OUT)
            YC_CX[-1,i] = YC_CX[0,i]
        YC_CX[:,-1] = YC_CX[:,0]

        phi = (YC_CX+phi_CX)*J
        phiOld = phi.copy() #update the time step
        print 'the mass of ', t,'time steps is',phiOld[:-1,:-1].sum()#,phi_AX[:-2,:-2].sum(),YC_AX[:-2,:-2].sum(),XC[:-2,:-2].sum()
    return phi

def advective(phiOld,cx,cy,nt,epsilon,dx,dy):

    #Calculate the number of grid points
    nx = len(phiOld[0,:]) - 1
    ny = len(phiOld[:,0]) - 1
    #Define the arrays used in COSMIC
    #advective-form one-dimensional flux difference of each cell
    XA = np.zeros_like(phiOld)
    YA = np.zeros_like(phiOld)
    #advective-form one-dimensional phi update
    phi_AY = np.zeros_like(phiOld)
    #the updated phi values
    phi = np.zeros_like(phiOld)
    #calculate the integer and remnant part of the Courant number
    phi_mid= np.zeros_like(phiOld)
    mass = np.zeros((ny,nx))
    for t in xrange(int(nt)):
        #Calculate the flux on each dimension
        for i in xrange(nx):
            phi_mid[:,i],mass[:,i] = PPM(phiOld[:,i],cy[:,i],ny,epsilon,dy)
            OUT = flux(ny,dy,cy[:,i],phi_mid[:-1,i],mass[:,i])
            YA[:-1,i] = adf(ny,cy[:,i],OUT)
        YA[-1,:] = YA[0,:]
        YA[:,-1] = YA[:,0]
        phi_AY = phiOld + YA
        for j in xrange(ny):
            phi_mid[j,:], mass[j,:]= PPM(phi_AY[j,:],cx[j,:],nx,epsilon,dx)
            OUT = flux(nx,dx,cx[j,:],phi_mid[j,:-1],mass[j,:])
            XA[j,:-1] = adf(nx,cx[j,:],OUT)
        XA[-1,:] = XA[0,:]
        XA[:,-1] = XA[:,0]
        phi = phi_AY+XA
    return phi