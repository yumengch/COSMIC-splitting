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

    return phi, phiExact, cx, cy, u, v, X_cntr, X_edge[:-1,:-1], Y_edge[:-1,:-1], Y[:-1,:-1], J, Y_C, J_p