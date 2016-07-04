def solid(x, y, xmin, ymin, nx, ny, Lx, Ly, nt, dt, distorted, t, change = False):
    '''Initial condition of solid body rotation test case'''
    def f(x,Lx):
        "the function that distort the grid"
        return np.where(x<=0.5*Lx,
                        (-1/np.sqrt(3))*x+Lx/4./np.sqrt(3)+0.5*Lx,
                        np.where(x>0.5*Lx,
                        1/np.sqrt(3)*(x-0.5*Lx)-Lx/4./np.sqrt(3)+0.5*Lx,0))

    def computational(X,Y,ymax):
        "The computational grids"
        fx = f(x,Lx)
        return np.where(Y>0.5*Ly, fx+(Y-0.5*Ly)*(Ly-fx)/(0.5*Ly),Y*fx/(0.5*Ly))
    #the necessary parameters
    dy,dx = Ly/ny,Lx/nx
    X,Y1 = np.meshgrid(x,y)
    phi = np.zeros((ny+1, nx+1))
    ymax = ymin+Ly
    
    #the computational domain
    if distorted:
        Y = computational(X,Y1,ymax)
        J = np.where(Y1>0.5*Ly, 0.5*Ly/(Ly-f(X,Lx)),0.5*Ly/f(X,Lx))
    else:
        Y = Y1
        J= 1

    #initialization of the tracer distribution
    A= np.pi*5/3/1000.
    r = 2500
    theta0 = np.pi/2.
    x0 = 0.5*nx*dx+r*np.cos(theta0)
    y0 = 0.5*ny*dy+ r*np.sin(theta0)
    phi = np.exp(- 0.5*(((X-x0)/500)**2 + ((Y-y0)/500)**2))

    #initialization of the streamfunction
    if distorted:
        Y = computational(X-0.5*dx,Y1-0.5*dy,ymax)
    r = ((X-0.5*dx)-(0.5*nx)*dx)**2+((Y-0.5*dy)-0.5*Ly)**2
    r1 = 4000.
    r2 = 4500.
    psi = A*r
    Y = computational(X,Y1,ymax)
    return phi,X,Y,J,psi

def terrain(x,z,xmin,zmin,nx,nz,Lx,Lz,t,nt,dt):
    def h(x):
        "The mountain height as a function of x"
        h0 = 6e3
        a = 25e3
        lam = 8e3
        return np.where(np.abs(x) <= a,
                        h0*np.cos(np.pi*x/lam)**2 * np.cos(np.pi*x*0.5/a)**2,
                        0)

    def sigmaHeights(x,Z,zmax):
        "The height of sigma coordinate Z"
        hx = h(x)
        return hx + Z*(zmax - hx)/zmax

    # parameters of the flow
    u0 = 10.
    z1 = 7e3
    z2 = 8e3
    dz,dx = Lz/nz,Lx/nx
    X,Z1 = np.meshgrid(x,z)
    phi = np.zeros((nz+1, nx+1))
    zmax = zmin+Lz
    # 1d arrays (x and z locations for different variables
    Z = sigmaHeights(X,Z1,zmax)
    psi= np.zeros((nz+1,nx+1))
    Ax = 25e3
    Az = 3e3
    x0 = -50e3
    z0 = 12e3
    r = np.sqrt(((X-x0)/Ax)**2 + ((Z-z0)/Az)**2)
    phi = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 
    J = zmax/(zmax - h(X-0.5*dx))
    J[0,:] = J[-1,:]
    J[:,0] = J[:,-1]
    Z = sigmaHeights(X-0.5*dx,Z1-0.5*dz,zmax)
    # Stream function at doubly staggered grid points
    for i in xrange(nx+1):
        for j in xrange(nz+1):   
            # Z[j,i] = Z1[j,i]-0.5*dz
            if (Z[j,i] <= z1):
                psi[j,i] = 0
            elif (Z[j,i] <= z2) and (Z[j,i]>z1):
                psi[j,i] = -0.5*u0*(Z[j,i]-z1-(z2-z1)/np.pi*np.sin(np.pi*((Z[j,i])-z1)/(z2-z1)))
                # psi[j,i] = -0.5*u0*(((z2-z1)*np.sin(np.pi*(Z[j,i]-z1)/(z2-z1))/np.pi)+z1-Z[j,i])
            else:
                psi[j,i] = -0.5*u0*(2*(Z[j,i]) - z1 - z2)
                # psi[j,i] = u0*Z[j,i]
    # print Z
    psi[:,0] = psi[:,-1]
    Z = sigmaHeights(X,Z1,zmax)
    return phi,X,Y,J,psi
    
def deformational(x,y,xmin,ymin,nx,ny,Lx,Ly,nt,dt, distorted, t, change):
    '''Initial condition of deformational flow test case'''
    def f(x,Lx):
        "the function that distort the grid"
        return np.where(x<=0.25*Lx,
                        -(1./np.sqrt(3))*(x)+Lx/8./np.sqrt(3),
                        np.where((x<0.5*Lx)&(x>=0.25*Lx),
                        1./np.sqrt(3)*(x-0.25*Lx)-Lx/8./np.sqrt(3),
                        np.where((x<0.75*Lx)&(x>=0.5*Lx),
                        -(1./np.sqrt(3))*(x-0.5*Lx)+Lx/8./np.sqrt(3),1./np.sqrt(3)*(x-0.75*Lx)-Lx/8./np.sqrt(3))))
     
    def computational(X,Y,ymax):
        "The computational grids"
        fx = f(X,Lx)
        return np.where(Y>=0, fx+Y*(1-fx/ymax),fx+Y*(1-fx/ymin))

    #the necessary parameters
    psi = np.zeros((ny+1,nx+1))
    phi = np.zeros((ny+1,nx+1))
    X,Y1 = np.meshgrid(x,y)
    dx,dy = Lx/nx, Ly/ny
    ymax = ymin+Ly
    u0 = 10.

    #the initialization of the test case for distorted grids:
    if distorted:
        Y = computational(X,Y1,ymax)
        J = np.where(Y1>=0, ymax/(ymax-f(x,Lx)),ymin/(ymin-f(x,Lx)))
        x0 = 5.*Lx/12.
        y0 = 0.
        x1 = 7.*Lx/12.
        Y = computational(X,Y1,ymax)
        phi = 0.95*np.exp(- 5.*((X-x0)**2 + (Y-y0)**2)) + 0.95*np.exp(- 5.*((X-x1)**2+ (Y-y0)**2))
        Y = computational(X-0.5*dx,Y1-0.5*dy,ymax)       
        psi = (u0/nt/dt)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*(((X-0.5*dx)/Lx) - (float(t)/nt))))**2)*((np.cos(np.pi*(Y-0.5*dy)/Lx))**2)      
        Y = computational(x,Y1,ymax)
    
    #if the flow is time-dependent
        if change:
            Y = computational(X-0.5*dx,Y1-0.5*dy,ymax)
            psi = (u0/nt/dt)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*(((X-0.5*dx)/Lx) - (float(t)/nt))))**2)*((np.cos(Y-0.5*dy))**2)*np.cos(np.pi*t/nt)-Lx*(Y-0.5*dy)/(nt*dt)
            return psi
        else:
            return phi,X,Y,J,psi
    
    #for orthogonal grids
    else:
        J = np.ones((ny+1,nx+1))
        Y = Y1
        x0 = 5.*Lx/12.
        y0 = 0.
        x1 = 7.*Lx/12.
        phi = 0.95*np.exp(- 5.*((X-x0)**2 + (Y-y0)**2)) + 0.95*np.exp(- 5.*((X-x1)**2+ (Y-y0)**2))

        psi = (u0/nt/dt)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*((X-0.5*dx)/Lx - float(t)/nt)))**2)*((np.cos(np.pi*(Y-0.5*dy)/Ly))**2)*np.cos(np.pi*t/nt)-Lx*(Y-0.5*dy)/(nt*dt)
        if change:
            psi = (u0/nt/dt)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*((X-0.5*dx)/Lx - float(t)/nt)))**2)*((np.cos(np.pi*(Y-0.5*dy)/Ly))**2)*np.cos(np.pi*t/nt)-Lx*(Y-0.5*dy)/(nt*dt)
            return psi
        else:
            return phi,X,Y,J,psi