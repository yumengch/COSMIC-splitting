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
def solid(x_edge, y_edge, x_cntr, y_cntr, t, nt, dt, mesh, change):
#---------------------------------------------------------------------------------
# Author: Yumeng Chen
# function: solid body rotation test case with angular velocity of 2A 
# and Gaussian distribution
#---------------------------------------------------------------------------------

    
    #-----------------------
    # Basic grid information
    #-----------------------
    Lx, Ly = np.max(x_edge) - np.min(x_edge),np.max(y_edge) - np.min(y_edge)
    nx, ny = len(x_cntr), len(y_cntr)
    dx, dy = Lx/nx, Ly/ny
    X_edge , Y_edge   = np.meshgrid(x_edge,y_edge)
    X_cntr, Y_cntr  = np.meshgrid(x_cntr,y_cntr)
    ymax = np.max(y_edge)
    phi = np.zeros_like(X_cntr)
    #----------------------------------------------
    # initial tracer position and angular velocity
    #----------------------------------------------
    A= np.pi*5/3/1000.
    r = 2500.
    theta0 = np.pi/2.
    x0 = 0.5*Lx + r*np.cos(theta0)
    y0 = 0.5*Ly + r*np.sin(theta0)

    #-----------------------------
    # initial tracer distribution
    #-----------------------------
    if mesh == 'quad':
        Y_C = compt_to_phys_SB(Y_cntr, Ly, f_quad(X_cntr,Lx, Ly))
    elif mesh == 'V':
        Y_C = compt_to_phys_SB(Y_cntr, Ly, f_V(X_cntr,Lx, Ly))
    else:
        Y_C = Y_cntr

    # tracer distribution:
    phi = np.exp(- 0.5*(((X_cntr-x0)/500)**2 + ((Y_C-y0)/500)**2))
    # phi[:,:] = 1.0
    #---------------------
    # Analytical solution 
    #---------------------
    phiExact = np.array([phi,phi,phi,phi,phi,phi])
    it = [int(nt/6.),int(nt/3.),int(nt/2.),int(2*nt/3.),int(5*nt/6), nt]
    for i in xrange(len(it)):
        x0 = 0.5*Lx + r*np.cos(theta0+2*A*it[i]*dt)
        y0 = 0.5*Ly + r*np.sin(theta0+2*A*it[i]*dt)
        phiExact[i] = np.exp(- 0.5*(((X_cntr - x0)/(500))**2 + ((Y_C - y0)/(500))**2))

    #----------------------
    # streamfunction field
    #----------------------
    if mesh == 'quad':
        Y = compt_to_phys_SB(Y_edge, Ly, f_quad(X_edge,Lx, Ly))
    elif mesh == 'V':
        Y = compt_to_phys_SB(Y_edge, Ly, f_V(X_edge,Lx, Ly))
    else:
        Y = Y_edge
    # streamfunction definition    
    r = np.sqrt((X_edge - 0.5*Lx)**2+(Y - 0.5*Ly)**2)
    psi = A*r**2


    #-----------------------------
    # dX and dY in physical domain
    #-----------------------------
    dX = np.zeros([ny,nx])
    dY = np.zeros([ny,nx])
    dY_J = np.zeros([ny,nx])
    J =  np.zeros([ny, nx])
    J_p =  np.zeros([ny, nx])
    dY = Y[1:, :-1] - Y[:-1, :-1]
    dX = X_edge[:-1, 1:] - X_edge[:-1, :-1]
    dY_J[1:, :] = Y_C[1:, :] - Y_C[:-1, :]
    dY_J[0, :]  = Y_C[0, :] + dY_J[-1, :]
    J = dy/dY
    J_p = dy/dY_J

    #-----------------------------
    # Courant number and velocity
    #-----------------------------
    u = np.zeros([ny, nx])
    v = np.zeros([ny, nx])  
    U = np.zeros([ny, nx])
    V = np.zeros([ny, nx])

    #-----------------------------
    # velocity in physical domain
    #-----------------------------
    u = -(psi[1:,:-1]-psi[:-1,:-1])/dY
    v =  (psi[:-1,1:]-psi[:-1,:-1])/dX
  

    #-----------------------------
    # Courant number in computational domain
    #-----------------------------
    U = -(psi[1:,:-1]-psi[:-1,:-1])/dy
    V =  (psi[:-1,1:]-psi[:-1,:-1])/dx
    # plt.figure(1)
    # plt.clf()
    # plt.quiver(X_cntr,Y_C,U,V)
    # plt.show()
    cx = U*dt/dx
    cy = V*dt/dy
    dudx = (1/J)*u - U
    dudx = ((1/J[:-1,1:])*u[:-1,1:]-(1/J[:-1,:-1])*u[:-1,:-1])/dx
    dvdy = (v[1:,:-1]-v[:-1,:-1])/dY
    # dudx= J[:-1,:-1]*(cx[:-1,1:]-cx[:-1,:-1])/dt
    # dvdy = J[:-1,:-1]*(cy[1:,:-1]-cy[:-1,:-1])/dt
    for j in xrange(ny):
        print 'divergence :', np.max(dudx[j,:])
    # print "d", np.max(dvdy)
    return phi, phiExact, cx, cy, u, v, X_cntr, X_edge[:-1,:-1], Y_edge[:-1,:-1], Y[:-1,:-1], J, Y_C, J_p

def orography(x_edge, z_edge, x_cntr, z_cntr, t, nt, dt, mesh, change):
#---------------------------------------------------------------------------------
# Author: Yumeng Chen
# function: horizontal advection over orography
# source: Schar, Christoph, et al. "A new terrain-following vertical coordinate 
#           formulation for atmospheric prediction models." Monthly Weather Review 
#                130.10 (2002): 2459-2480.
#---------------------------------------------------------------------------------

    #-----------------------
    # Basic grid information
    #-----------------------
    Lx, Lz = np.max(x_edge) - np.min(x_edge),np.max(z_edge) - np.min(z_edge)
    nx, nz = len(x_cntr), len(z_cntr)
    dx, dz = Lx/nx, Lz/nz
    X_edge , Z_edge   = np.meshgrid(x_edge,z_edge)
    X_cntr, Z_cntr  = np.meshgrid(x_cntr,z_cntr)
    zmax = np.max(z_edge)
    psi = np.zeros_like(X_edge)
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
    Z_C = compt_to_phys_oro(X_cntr,Z_cntr,zmax, h0)
    r = np.sqrt(((X_cntr - x0)/Ax)**2 + ((Z_C - z0)/Az)**2)
    phi = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 

    #-----------------------------
    # Analytical solution 
    #-----------------------------
    phiExact = np.array([phi,phi])
    it = [int(nt/2.), nt]
    for i in xrange(len(it)):
        x0 = -50e3 + it[i]*dt*u0
        r = np.sqrt(((X_cntr - x0)/Ax)**2 + ((Z_C - z0)/Az)**2)
        phiExact[i] = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 

    #-----------------------------
    # streamfunction field
    #-----------------------------
    Z = compt_to_phys_oro(X_edge,Z_edge,zmax, h0)
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
    # dX and dZ in physical domain
    #-----------------------------
    dX = np.zeros([nz,nx])
    dZ = np.zeros([nz,nx])
    dZ_J = np.zeros([nz,nx])
    J =  np.zeros([nz, nx])
    J_p =  np.zeros([nz, nx])
    dZ = Z[1:, :-1] - Z[:-1, :-1]
    dX = X_edge[:-1, 1:] - X_edge[:-1, :-1]
    dZ_J[1:, :] = Z_C[1:, :] - Z_C[:-1, :]
    dZ_J[0, :]  = Z_C[0, :] + dZ_J[-1, :]
    J = dz/dZ
    J_p = dz/dZ_J

    #-----------------------------
    # Courant number and velocity
    #-----------------------------
    u = np.zeros([nz, nx])
    v = np.zeros([nz, nx])  
    U = np.zeros([nz, nx])
    V = np.zeros([nz, nx])

    #-----------------------------
    # velocity in physical domain
    #-----------------------------
    u = -(psi[1:,:-1]-psi[:-1,:-1])/dZ
    v =  (psi[:-1,1:]-psi[:-1,:-1])/dX
  

    #-----------------------------
    # Courant number in computational domain
    #-----------------------------
    U = -(psi[1:,:-1]-psi[:-1,:-1])/dz
    V =  (psi[:-1,1:]-psi[:-1,:-1])/dx
    # plt.figure(1)
    # plt.clf()
    # plt.quiver(X_cntr,Y_C,U,V)
    # plt.show()
    cx = U*dt/dx
    cy = V*dt/dz

    return phi, phiExact, cx, cy, u, v, X_cntr, X_edge[:-1,:-1], Z_edge[:-1,:-1], Z[:-1,:-1], J, Z_C, J_p

def deform(x_edge, y_edge, x_cntr, y_cntr, t, nt, dt, mesh, change):
#----------------------------------------------------------------------------------------
# Author: Yumeng Chen
# adapted from: Lauritzen, P. H., W. C. Skamarock, M. J. Prather, and M. A. Taylor, 2012:
#  A standard test case suite for two-dimensional linear transport on the sphere.
#    Geosci. Model Dev., 5, 887–901, doi:10.5194/gmd-5-887-2012.
#----------------------------------------------------------------------------------------
        
    #-----------------------
    # Basic grid information
    #-----------------------
    Lx, Ly = np.max(x_edge) - np.min(x_edge),np.max(y_edge) - np.min(y_edge)
    nx, ny = len(x_cntr), len(y_cntr)
    dx, dy = Lx/nx, Ly/ny
    X_edge , Y_edge   = np.meshgrid(x_edge,y_edge)
    X_cntr, Y_cntr  = np.meshgrid(x_cntr,y_cntr)
    ymin = np.min(y_edge)
    ymax = np.max(y_edge)
    phi = np.zeros_like(X_cntr)

    #--------------------------------------
    # initial stream and tracer parameters
    #--------------------------------------
    u0 = 10.
    x0 = 5.*Lx/12.
    y0 = 0
    x1 = 7.*Lx/12.
    T = 5.0

    if mesh == 'W':
        Y = compt_to_phys_deform(X_edge, Y_edge, Lx, ymin, ymax)
    else:
        Y = Y_edge
    #---------------------------------
    # time varying velocity in 
    #---------------------------------
    if change:
        psi = (u0/T)*((Lx/2./np.pi)**2)*((np.sin(2*np.pi*((X_edge/Lx) - (float(t)/nt))))**2)*((np.cos(Y))**2)*np.cos(np.pi*t/nt)-Lx*Y/T
    
        #-----------------------------
        # dX and dY in physical domain
        #-----------------------------
        dX = np.zeros([ny,nx])
        dY = np.zeros([ny,nx])
        dY = Y[1:, :-1] - Y[:-1, :-1]
        dX = X_edge[:-1, 1:] - X_edge[:-1, :-1]

        #-----------------------------
        # Courant number and velocity
        #-----------------------------
        u = np.zeros([ny, nx])
        v = np.zeros([ny, nx])  
        U = np.zeros([ny, nx])
        V = np.zeros([ny, nx])

        #-----------------------------
        # velocity in physical domain
        #-----------------------------
        u = -(psi[1:,:-1]-psi[:-1,:-1])/dY
        v =  (psi[:-1,1:]-psi[:-1,:-1])/dX
      

        #-----------------------------
        # Courant number in computational domain
        #-----------------------------
        U = -(psi[1:,:-1]-psi[:-1,:-1])/dy
        V =  (psi[:-1,1:]-psi[:-1,:-1])/dx
        # plt.figure(1)
        # plt.clf()
        # plt.quiver(X_cntr,Y_C,U,V)
        # plt.show()
        cx = U*dt/dx
        cy = V*dt/dy
        return u, v, cx, cy

    #-----------------------------
    # Initial tracer distribution
    #-----------------------------
    if mesh == 'W':
        Y_C = compt_to_phys_deform(X_cntr, Y_cntr, Lx, ymin, ymax)
    else:
        Y_C = Y_cntr
    phi = 0.95*np.exp(- 5*((X_cntr - x0)**2 + (Y_C - y0)**2)) + 0.95*np.exp(- 5*((X_cntr - x1)**2+ (Y_C - y0)**2))


    psi = (u0/T)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*((X_edge/Lx) - (float(0)/nt))))**2)*((np.cos(np.pi*Y/Lx))**2)      
       
    #-----------------------------
    # Expected analytical solution
    #-----------------------------
    phiExact = phi

    # #-----------------------------
    # # Initial stream field
    # #-----------------------------
    # if mesh == 'W':
    #     Y = comput_deform(X-0.5*dx,Y1-0.5*dy,ymax)
    # else:
    #     Y = Y1
 
    #-----------------------------
    # dX and dY in physical domain
    #-----------------------------
    dX = np.zeros([ny,nx])
    dY = np.zeros([ny,nx])
    dY_J = np.zeros([ny,nx])
    J =  np.zeros([ny, nx])
    J_p =  np.zeros([ny, nx])
    dY = Y[1:, :-1] - Y[:-1, :-1]
    dX = X_edge[:-1, 1:] - X_edge[:-1, :-1]
    dY_J[1:, :] = Y_C[1:, :] - Y_C[:-1, :]
    dY_J[0, :]  = Y_C[0, :] + dY_J[-1, :]
    J = dy/dY
    J_p = dy/dY_J

    #-----------------------------
    # Courant number and velocity
    #-----------------------------
    u = np.zeros([ny, nx])
    v = np.zeros([ny, nx])  
    U = np.zeros([ny, nx])
    V = np.zeros([ny, nx])

    #-----------------------------
    # velocity in physical domain
    #-----------------------------
    u = -(psi[1:,:-1]-psi[:-1,:-1])/dY
    v =  (psi[:-1,1:]-psi[:-1,:-1])/dX
  

    #-----------------------------
    # Courant number in computational domain
    #-----------------------------
    U = -(psi[1:,:-1]-psi[:-1,:-1])/dy
    V =  (psi[:-1,1:]-psi[:-1,:-1])/dx
    # plt.figure(1)
    # plt.clf()
    # plt.quiver(X_cntr,Y_C,U,V)
    # plt.show()
    cx = U*dt/dx
    cy = V*dt/dy

    return phi, phiExact, cx, cy, u, v, X_cntr, X_edge[:-1,:-1], Y_edge[:-1,:-1], Y[:-1,:-1], J, Y_C, J_p