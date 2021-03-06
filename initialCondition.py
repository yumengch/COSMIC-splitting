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
    dX = X_edge[:-1, 1:] - X_edge[:-1, :-1]
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
        Y_C = compt_to_phys_SB(Y_cntr, Ly, 0, f_quad(X_cntr,Lx, Ly))
    elif mesh == 'V':
        Y_C = compt_to_phys_SB(Y_cntr, Ly, 0, f_V(X_cntr,Lx, Ly))
    else:
        Y_C = Y_cntr

    # tracer distribution:
    phi = np.exp(- 0.5*(((X_cntr-x0)/500)**2 + ((Y_C-y0)/500)**2))

    #---------------------
    # Analytical solution 
    #---------------------
    phiExact = np.array([phi,phi,phi,phi,phi,phi])
    it = [int(nt/6.),int(nt/3.),int(nt/2.),int(2*nt/3.),int(5*nt/6), nt]
    for i in xrange(len(it)):
        x0 = 0.5*Lx + r*np.cos(theta0+2*A*it[i]*dt)
        y0 = 0.5*Ly + r*np.sin(theta0+2*A*it[i]*dt)
        phiExact[i] = np.exp(- 0.5*(((X_cntr - x0)/(500))**2 + ((Y_C - y0)/(500))**2))
    #-----------------------------
    # Jacobian for v
    #-----------------------------
    dY_s = np.zeros([ny, nx])
    dY_s[1:,:] = Y_C[1:, :] - Y_C[:-1, :]
    dY_s[0,:] = Y_C[0, :] - y_edge[0] + ymax - Y_C[-1, :]
    J_s = dy/dY_s

    #----------------------
    # streamfunction field
    #----------------------
    if mesh == 'quad':
        Y = compt_to_phys_SB(Y_edge, Ly, 0, f_quad(X_edge,Lx, Ly))
    elif mesh == 'V':
        Y = compt_to_phys_SB(Y_edge, Ly, 0, f_V(X_edge,Lx, Ly))
    else:
        Y = Y_edge
    # streamfunction definition    
    r = np.sqrt((X_edge - 0.5*Lx)**2+(Y - 0.5*Ly)**2)
    psi = A*r**2
    #-----------------------------
    # Jacobian for u
    #-----------------------------
    dY_w = Y[1:, :-1] - Y[:-1, :-1]
    J_w = dy/dY_w

    #------------------
    # Jacobian for phi
    #------------------
    if mesh == 'quad':
        Y_J = compt_to_phys_SB(Y_edge[:,:-1], Ly, 0, f_quad(np.vstack([X_cntr, x_cntr]),Lx, Ly))
    elif mesh == 'V':
        Y_J = compt_to_phys_SB(Y_edge[:,:-1], Ly, 0, f_V(np.vstack([X_cntr, x_cntr]),Lx, Ly))
    else:
        Y_J = Y_edge[:,:-1]

    dY = Y_J[1:, :] - Y_J[:-1, :]
    J = dy/dY

    #-----------------------------
    # velocity in computational domain
    #-----------------------------
    U = -J_w*(psi[1:,:-1]-psi[:-1,:-1])/dy
    V =  J_s*(psi[:-1,1:]-psi[:-1,:-1])/dx

    cx = U*dt/dx
    cy = V*dt/dy
    dcx = np.round(np.max(cx[:-1,1:]/J_w[:-1,1:] - cx[:-1,:-1]/J_w[:-1,:-1]), 5)
    dcy = np.round(np.max(cy[1:,:-1]/J_s[1:,:-1] - cy[:-1,:-1]/J_s[:-1,:-1]), 5)
    print "the deformational Courant number in x:", dcx, "in y: ", dcy
    print "the max Courant number in x:", np.round(np.max(cx/J_w), 5), "in y: ", np.round(np.max(cy/J_s), 5)
    return phi, phiExact, cx, cy, J, J_s, J_w

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
    dX = X_edge[:-1, 1:] - X_edge[:-1, :-1]
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
    # Jacobian for v
    #-----------------------------
    dZ_s = np.zeros([nz, nx])
    dZ_s[1:,:] = Z_C[1:, :] - Z_C[:-1, :]
    dZ_s[0,:] = Z_C[0, :] - Z_C[-1, :]
    J_s = dz/dZ_s

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
    # Jacobian for u
    #-----------------------------
    dZ_w = Z[1:, :-1] - Z[:-1, :-1]
    J_w = dz/dZ_w
    #-----------------------------
    # Jacobian for phi
    #-----------------------------
    Z_J = compt_to_phys_oro(np.vstack([X_cntr, x_cntr]),Z_edge[:,:-1],zmax, h0)
    dZ = Z_J[1:, :] - Z_J[:-1, :]
    J = dz/dZ

    #-----------------------------
    # Courant number in computational domain
    #-----------------------------
    U = -J_w*(psi[1:,:-1]-psi[:-1,:-1])/dz
    V =  J_s*(psi[:-1,1:]-psi[:-1,:-1])/dx

    cx = U*dt/dx
    cy = V*dt/dz
    dcx = np.round(np.max(cx[:-1,1:]/J_w[:-1,1:] - cx[:-1,:-1]/J_w[:-1,:-1]), 5)
    dcy = np.round(np.max(cy[1:,:-1]/J_s[1:,:-1] - cy[:-1,:-1]/J_s[:-1,:-1]), 5)
    print "the deformational Courant number in x:", dcx, "in y: ", dcy
    print "the max Courant number in x:", np.round(np.max(cx/J_w), 5), "in y: ", np.round(np.max(cy/J_s), 5)
    return phi, phiExact, cx, cy, J, J_s, J_w

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
    dX = X_edge[:-1, 1:] - X_edge[:-1, :-1]
    #--------------------------------------
    # initial stream and tracer parameters
    #--------------------------------------
    u0 = 10.
    x0 = 5.*Lx/12.
    y0 = 0
    x1 = 7.*Lx/12.
    T = 5.0
    #-----------------------------
    # Initial tracer distribution
    #-----------------------------
    if mesh == 'W':
        Y_C = compt_to_phys_deform(X_cntr, Y_cntr, Lx, ymin, ymax)
    else:
        Y_C = Y_cntr
    # phi = np.zeros([ny, nx])
    phi = 0.95*np.exp(- 5*((X_cntr - x0)**2 + (Y_C - y0)**2)) + 0.95*np.exp(- 5*((X_cntr - x1)**2+ (Y_C - y0)**2))
    phiExact = phi

    #-----------------------------
    # Jacobian for v
    #-----------------------------
    dY_s = np.zeros([ny, nx])
    dY_s[1:,:] = Y_C[1:, :] - Y_C[:-1, :]
    dY_s[0,:] = Y_C[0, :] - Y_C[-1, :]
    J_s = dy/dY_s

    if mesh == 'W':
        Y_E = compt_to_phys_deform(X_edge, Y_edge, Lx, ymin, ymax)
    else:
        Y_E = Y_edge
    #-----------------------------
    # Jacobian for u
    #-----------------------------
    dY_w = Y_E[1:, :-1] - Y_E[:-1, :-1]
    J_w = dy/dY_w

    #------------------
    # Jacobian for phi
    #-----------------
    if mesh == 'W':
        Y = compt_to_phys_deform(np.vstack([X_cntr, x_cntr]), Y_edge[:,:-1], Lx, ymin, ymax)
    else:
        Y = Y_edge[:,:-1]
    dY = Y[1:, :] - Y[:-1, :]
    J = dy/dY
    #---------------------------------
    # time varying velocity in 
    #---------------------------------
    if change:
        psi = (u0/T)*((Lx/2./np.pi)**2)*((np.sin(2*np.pi*((X_edge/Lx) - (float(t)/nt))))**2)*((np.cos(Y_E))**2)*np.cos(np.pi*t/nt)-Lx*Y_E/T
        #-----------------------------
        # Courant number in computational domain
        #-----------------------------
        U = -J_w*(psi[1:,:-1]-psi[:-1,:-1])/dy
        V =  J_s*(psi[:-1,1:]-psi[:-1,:-1])/dx

        cx = U*dt/dx
        cy = V*dt/dy
        dcx = np.round(np.max(cx[:-1,1:]/J_w[:-1,1:] - cx[:-1,:-1]/J_w[:-1,:-1]), 5)
        dcy = np.round(np.max(cy[1:,:-1]/J_s[1:,:-1] - cy[:-1,:-1]/J_s[:-1,:-1]), 5)
        return cx, cy

    psi = (u0/T)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*((X_edge/Lx) - (float(1)/nt))))**2)*((np.cos(np.pi*Y_E/Ly))**2) -Lx*Y_E/T

    #-----------------------------
    # Courant number in computational domain
    #-----------------------------
    U = -J_w*(psi[1:,:-1]-psi[:-1,:-1])/dy
    V =  J_s*(psi[:-1,1:]-psi[:-1,:-1])/dx

    cx = U*dt/dx
    cy = V*dt/dy
    dcx = np.round(np.max(cx[:-1,1:]/J_w[:-1,1:] - cx[:-1,:-1]/J_w[:-1,:-1]), 5)
    dcy = np.round(np.max(cy[1:,:-1]/J_s[1:,:-1] - cy[:-1,:-1]/J_s[:-1,:-1]), 5)
    print "the max deformational Courant number in x:", dcx, "in y: ", dcy
    print "the max Courant number in x:", np.round(np.max(cx/J_w), 5), "in y: ", np.round(np.max(cy/J_s), 5)
    return phi, phiExact, cx, cy, J, J_s, J_w