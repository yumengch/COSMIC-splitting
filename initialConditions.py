    
    #---------------------------------------------------------------------------------
# Author: Yumeng Chen
# Initial conditions to test COSMIC splitting with solid body rotation test case,
# horizontal advection over orography test case and deformational flow test case.
# cx, cy is the Courant number in x- and y- direction
# mesh: 'orthogo' means orthogonal grids (defaulted in solid body rotation and 
#       deformational flow test case),
#       'quad', 'V' is used in solid body rotation test case for quadratic mesh 
#        and V shape mesh
#       'high' and 'low' means h0 = 6km and 3km as orography respectively 
#       'W' means W shape mesh in deformational flow test case
# J is the Jacobian for coordinate transformation
#---------------------------------------------------------------------------------
def solid(x, y, xmin, ymin, nx, ny, Lx, Ly ,t, nt, dt, mesh, change):
#---------------------------------------------------------------------------------
# Author: Yumeng Chen
# function: solid body rotation test case with angular velocity of 2A 
# and Gaussian distribution
#---------------------------------------------------------------------------------
    def f_quad(x,Lx):
        #-------------------------
        # quadratic mesh function
        #-------------------------
        return 2*(x-0.5*Lx)**2/Lx/np.sqrt(3) + 0.5*Ly-0.5*0.25*Lx/np.sqrt(3)
    def f_V(x,Lx):
        #-------------------------
        # V shape mesh function
        #-------------------------
        return np.where(x<=0.5*Lx,
                        (-1/np.sqrt(3))*x+Lx/4./np.sqrt(3)+0.5*Lx,
                        np.where(x>0.5*Lx,
                        1/np.sqrt(3)*(x-0.5*Lx)-Lx/4./np.sqrt(3)+0.5*Lx,0))

    def comput_SB(X, Y, ymax, fx):
        #-------------------------------------
        # non-orthogonal computational domain 
        # for solid body rotation test case
        #-------------------------------------
        return np.where(Y>0.5*Ly, fx+(Y-0.5*Ly)*(Ly-fx)/(0.5*Ly),Y*fx/(0.5*Ly))
    
    #-----------------------
    # Basic grid information
    #-----------------------
    dy,dx = Ly/ny,Lx/nx
    X,Y1 = np.meshgrid(x,y)
    ymax = ymin+Ly

    #----------------------------------------------
    # initial tracer position and angular velocity
    #----------------------------------------------
    A= np.pi*5/3/1000.
    r = 2500.
    theta0 = np.pi/2.
    x0 = 0.5*nx*dx+r*np.cos(theta0)
    y0 = 0.5*ny*dy+ r*np.sin(theta0)

    #-----------------------------
    # initial tracer distribution
    #-----------------------------
    if mesh == 'quad':
        Y = comput_SB(X,Y1,ymax,f_quad(X,Lx))
        J= np.where(Y1>0.5*Ly, 0.5*Ly/(Ly-f_quad(X,Lx)),0.5*Ly/f_quad(X,Lx))
    elif mesh == 'V':
        Y = comput_SB(X,Y1,ymax,f_V(X,Lx))
        J= np.where(Y1>0.5*Ly, 0.5*Ly/(Ly-f_V(X,Lx)),0.5*Ly/f_V(X,Lx))
    else:
        Y = Y1
        J = np.ones_like(X)
    # tracer distribution:
    phi = np.exp(- 0.5*(((X-x0)/500)**2 + ((Y-y0)/500)**2))

    # dx = X[:, 1:] - X[:, :-1]
    # print dx#, Y[:,0]
    #-----------------------------
    # Analytical solution 
    #-----------------------------
    phiExact = np.array([phi,phi,phi,phi,phi,phi])
    it = [int(nt/6.),int(nt/3.),int(nt/2.),int(2*nt/3.),int(5*nt/6), nt]
    for i in xrange(len(it)):
        x0 = 0.5*nx*dx+r*np.cos(theta0+2*A*it[i]*dt)
        y0 = 0.5*ny*dy+ r*np.sin(theta0+2*A*it[i]*dt)
        phiExact[i] = np.exp(- 0.5*(((X-x0) / (500))**2 + ((Y-y0) / (500))**2))

    #-----------------------------
    # streamfunction field
    #-----------------------------
    if mesh == 'quad':
        Y = comput_SB(X-0.5*dx,Y1-0.5*dy,ymax,f_quad(X,Lx))
    elif mesh == 'V':
        Y = comput_SB(X-0.5*dx,Y1-0.5*dy,ymax,f_V(X,Lx))
    else:
        Y = Y1
    # streamfunction definition    
    r = np.sqrt((X-0.5*Lx)**2+(Y-0.5*Ly)**2)
    psi = A*r**2



    if mesh == 'quad':
        Y = comput_SB(X,Y1,ymax,f_quad(X,Lx))
        # J= np.where(Y1>0.5*Ly, 0.5*Ly/(Ly-f_quad(X,Lx)),0.5*Ly/f_quad(X,Lx))
    elif mesh == 'V':
        Y = comput_SB(X,Y1,ymax,f_V(X,Lx))
        # J= np.where(Y1>0.5*Ly, 0.5*Ly/(Ly-f_V(X,Lx)),0.5*Ly/f_V(X,Lx))
    else:
        Y = Y1
        # J = np.ones_like(X)
    
    dx = np.zeros([ny+1, nx+1])
    dy = np.zeros([ny+1, nx+1])
    dx[:,:-1] = X[:, 1:] - X[:, :-1]
    dx[:,-1] = dx[:,0]
    dy[:-1, :] = Y[1:, :] - Y[:-1, :]
    dy[-1,:] = dy[0,:]

    #-----------------------------
    # Courant number
    #-----------------------------
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)

    # the velocity field
    u[:-1,:] = -(psi[1:,:]-psi[:-1,:])/dy[:-1,:]
    v[:,:-1] =  (psi[:,1:]-psi[:,:-1])/dx[:,:-1]
    u[-1,:] = u[0,:]
    v[:,-1] = v[:,0]
    # u[:,:] = 3.0
    # v[:,:] = 0.
    # #the Courant number
    # cx = u*dt/dx
    # cy = v*dt/dy

    return phi, phiExact, u, v, X, Y, J

def orography(x, z, xmin, zmin, nx, nz, Lx, Lz ,t, nt, dt, mesh, change):
#---------------------------------------------------------------------------------
# Author: Yumeng Chen
# function: horizontal advection over orography
# source: Schar, Christoph, et al. "A new terrain-following vertical coordinate 
#           formulation for atmospheric prediction models." Monthly Weather Review 
#                130.10 (2002): 2459-2480.
#---------------------------------------------------------------------------------
    def h(x, h0):
        #-------------------------
        # function for basic terrain
        # following coordinate (BTF)
        #-------------------------
        a = 25e3
        lam = 8e3
        return np.where(np.abs(x) <= a,
                        h0*np.cos(np.pi*x/lam)**2 * np.cos(np.pi*x*0.5/a)**2,
                        0)

    def sigmaHeights(x,Z,zmax):
        #------------------------------
        # computational domain for BTF
        #-----------------------------
        hx = h(x, h0)
        return hx + Z*(zmax - hx)/zmax

    #-----------------------
    # Basic grid information
    #-----------------------
    dz,dx = Lz/nz,Lx/nx
    X,Z1 = np.meshgrid(x,z+0.5*dz)
    zmax = zmin+Lz
    psi = np.zeros_like(X)

    #----------------------------------
    # initial flow and tracer parameters
    #----------------------------------
    u0 = 10.
    Ax = 25e3
    Az = 3e3
    x0 = -50e3


    if mesh == 'high':
        h0 = 6e3
        z0 = 12e3
        z1 = 7e3
        z2 = 8e3
    else:
        h0 = 3e3
        z0 = 9e3
        z1 = 4e3
        z2 = 5e3

    #-----------------------------
    # initial tracer distribution
    #-----------------------------
    Z = sigmaHeights(X,Z1,zmax)    
    r = np.sqrt(((X-x0)/Ax)**2 + ((Z-z0)/Az)**2)
    phi = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 

    #-----------------------------
    # Analytical solution 
    #-----------------------------
    phiExact = np.array([phi,phi])
    it = [int(nt/2.), nt]
    for i in xrange(len(it)):
        x0 = -50e3+it[i]*dt*u0
        r = np.sqrt(((X-x0)/Ax)**2 + ((Z-z0)/Az)**2)
        phiExact[i] = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 

    #-----------------------------
    # Jacobian
    #-----------------------------
    J = zmax/(zmax - h(X-0.5*dx,h0))
    J[0,:] = J[-1,:]
    J[:,0] = J[:,-1]

    #-----------------------------
    # streamfunction field
    #-----------------------------
    Z = sigmaHeights(X-0.5*dx,Z1-0.5*dz,zmax)
    # Stream function at doubly staggered grid points
    for i in xrange(nx+1):
        for j in xrange(nz+1):   
            # Z[j,i] = Z1[j,i]-0.5*dz
            if (Z[j,i] <= z1):
                psi[j,i] = 0.
            elif (Z[j,i] <= z2) and (Z[j,i]>z1):
                psi[j,i] = -0.5*u0*(Z[j,i]-z1-(z2-z1)/np.pi*np.sin(np.pi*((Z[j,i])-z1)/(z2-z1)))
            else:
                psi[j,i] = -0.5*u0*(2*(Z[j,i]) - z1 - z2)

    #-----------------------------
    # Courant number
    #-----------------------------
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)
    # the velocity field
    u[:-1,:] = -(psi[1:,:]-psi[:-1,:])/dz
    v[:,:-1] =  (psi[:,1:]-psi[:,:-1])/dx
    u[-1,:] = u[0,:]
    v[:,-1] = v[:,-2]
    
    #the Courant number
    cx = u*dt/dx
    cy = v*dt/dz

    return phi, phiExact, u, v, cx, cy, J

def deform(x, y, xmin, ymin, nx, ny, Lx, Ly ,t, nt, dt, mesh, change):
#----------------------------------------------------------------------------------------
# Author: Yumeng Chen
# adapted from: Lauritzen, P. H., W. C. Skamarock, M. J. Prather, and M. A. Taylor, 2012:
#  A standard test case suite for two-dimensional linear transport on the sphere.
#    Geosci. Model Dev., 5, 887–901, doi:10.5194/gmd-5-887-2012.
#----------------------------------------------------------------------------------------
    def f_W(x,Lx):
        #-------------------------
        # W shape mesh function
        #-------------------------
        return np.where(x<=0.25*Lx,
                        -(1./np.sqrt(3))*(x)+Lx/8./np.sqrt(3),
                        np.where((x<0.5*Lx)&(x>=0.25*Lx),
                        1./np.sqrt(3)*(x-0.25*Lx)-Lx/8./np.sqrt(3),
                        np.where((x<0.75*Lx)&(x>=0.5*Lx),
                        -(1./np.sqrt(3))*(x-0.5*Lx)+Lx/8./np.sqrt(3),1./np.sqrt(3)*(x-0.75*Lx)-Lx/8./np.sqrt(3))))
        
    def comput_deform(X,Y,ymax):
        #-------------------------
        # computational domain for
        # W shape mesh
        #-------------------------
        fx = f_W(X,Lx)
        return np.where(Y>=0, fx+Y*(1-fx/ymax),fx+Y*(1-fx/ymin))

    #-----------------------
    # Basic grid information
    #-----------------------
    X,Y1 = np.meshgrid(x,y)
    dx,dy = Lx/nx, Ly/ny
    ymax = ymin+Ly

    #--------------------------------------
    # initial stream and tracer parameters
    #--------------------------------------
    u0 = 10.
    x0 = 5.*Lx/12.
    y0 = 0
    x1 = 7.*Lx/12.

    #---------------------------------
    # time varying velocity in 
    #---------------------------------
    if change:
        if mesh == 'W':
            Y = comput_deform(X-0.5*dx,Y1-0.5*dy,ymax)
            J = np.where(Y1>=0, ymax/(ymax-f_W(x,Lx)),ymin/(ymin-f_W(x,Lx)))
        else:
            Y = Y1
            J = np.ones_like(X)
        psi = (u0/nt/dt)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*((X/Lx) - (float(t)/nt))))**2)*((np.cos(Y))**2)*np.cos(np.pi*t/nt)-Lx*Y/(nt*dt)
        #-----------------------------
        # Courant number
        #-----------------------------
        u = np.zeros_like(psi)
        v = np.zeros_like(psi)
        # the velocity field
        u[:-1,:] = -(psi[1:,:]-psi[:-1,:])/dy
        v[:,:-1] =  (psi[:,1:]-psi[:,:-1])/dx
        v[:,-1],v[-1,:] = v[:,0],v[0,:]
        u[-1,:],u[:,-1] = u[0,:],u[:,0]
        
        #the Courant number
        cx = J*u*dt/dx
        cy = J*v*dt/dy

        return cx, cy

    #-----------------------------
    # Initial tracer distribution
    #-----------------------------
    if mesh == 'W':
        Y = comput_deform(X,Y1,ymax)
        J = np.where(Y1>=0, ymax/(ymax-f_W(x,Lx)),ymin/(ymin-f_W(x,Lx)))
    else:
        Y = Y1
        J = np.ones_like(X)
    phi = 0.95*np.exp(- 5*((X-x0)**2 + (Y-y0)**2)) + 0.95*np.exp(- 5*((X-x1)**2+ (Y-y0)**2))

    #-----------------------------
    # Expected analytical solution
    #-----------------------------
    phiExact = phi

    #-----------------------------
    # Initial stream field
    #-----------------------------
    if mesh == 'W':
        Y = comput_deform(X-0.5*dx,Y1-0.5*dy,ymax)
    else:
        Y = Y1
    psi = (u0/nt/dt)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*((X/Lx))))**2)*((np.cos(np.pi*(Y)/Lx))**2)  

    #-----------------------------
    # Courant number
    #-----------------------------
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)
    # the velocity field
    u[:-1,:] = -(psi[1:,:]-psi[:-1,:])/dy
    v[:,:-1] =  (psi[:,1:]-psi[:,:-1])/dx
    v[:,-1],v[-1,:] = v[:,0],v[0,:]
    u[-1,:],u[:,-1] = u[0,:],u[:,0]

    #the Courant number
    cx = u*dt/dx
    cy = v*dt/dy    
    
    return phi, phiExact, u, v, cx, cy, J
        