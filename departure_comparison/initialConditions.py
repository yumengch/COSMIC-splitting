def solid(x_edge, y_edge, x_cnter, y_cnter, t, nt, dt, mesh, change):
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

    def compt_to_phys_SB(Y, Ly, fx):
        #-------------------------------------
        # non-orthogonal computational domain 
        # for solid body rotation test case
        #-------------------------------------
        return np.where(Y>0.5*Ly, fx+(Y-0.5*Ly)*(Ly-fx)/(0.5*Ly),Y*fx/(0.5*Ly))

    def phys_to_compt_SB(Y, Ly, fx):
        #-------------------------------------
        # non-orthogonal computational domain 
        # for solid body rotation test case
        #-------------------------------------
        return np.where(Y>0.5*Ly, fx+(Y-0.5*Ly)*(Ly-fx)/(0.5*Ly),Y*fx/(0.5*Ly))
    
    #-----------------------
    # Basic grid information
    #-----------------------
    Lx, Ly = np.max(x_edge) - np.min(x_edge),np.max(y_edge) - np.min(y_edge)
    nx, ny = len(x_cnter), len(y_cnter)
    dx, dy = Lx/nx, Ly/ny
    X_edge , Y_edge   = np.meshgrid(x_edge,y_edge)
    X_cnter, Y_cnter  = np.meshgrid(x_cnter,y_cnter)
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
        Y_C = compt_to_phys_SB(Y_cnter, Ly, f_quad(X_cnter,Lx))
    elif mesh == 'V':
        Y_C = compt_to_phys_SB(Y_cnter, Ly, f_V(X_cnter,Lx))
    else:
        Y_C = Y_cnter

    # tracer distribution:
    phi = np.exp(- 0.5*(((X_cnter-x0)/500)**2 + ((Y_C-y0)/500)**2))

    #---------------------
    # Analytical solution 
    #---------------------
    phiExact = np.array([phi,phi,phi,phi,phi,phi])
    it = [int(nt/6.),int(nt/3.),int(nt/2.),int(2*nt/3.),int(5*nt/6), nt]
    for i in xrange(len(it)):
        x0 = 0.5*Lx + r*np.cos(theta0+2*A*it[i]*dt)
        y0 = 0.5*Ly + r*np.sin(theta0+2*A*it[i]*dt)
        phiExact[i] = np.exp(- 0.5*(((X_cnter - x0)/(500))**2 + ((Y_C - y0)/(500))**2))

    #----------------------
    # streamfunction field
    #----------------------
    if mesh == 'quad':
        Y = compt_to_phys_SB(Y_edge, Ly, f_quad(X_edge,Lx))
    elif mesh == 'V':
        Y = compt_to_phys_SB(Y_edge, Ly, f_V(X_edge,Lx))
    else:
        Y = Y_edge
    # streamfunction definition    
    r = np.sqrt((X_edge)**2+(Y)**2)
    psi = A*r**2


    #-----------------------------
    # dX and dY in physical domain
    #-----------------------------
    dX = np.zeros([ny,nx])
    dY = np.zeros([ny,nx])
    dY_J = np.zeros([ny+1,nx+1])
    J =  np.zeros([ny+1, nx+1])
    
    dY = Y[1:, :-1] - Y[:-1, :-1]

    dX = X_edge[:-1, 1:] - X_edge[:-1, :-1]
    # dY_J[1:-1, :-1] = Y_C[1:, :] - Y_C[:-1, :]
    # dY_J[-1, :-1] = dY_J[0, :-1]  = 0.5*(dY[0,:] + dY[-1,:])
    # dY_J[:, -1] = dY_J[:,0] 
    J[:-1,:-1] = dy/dY
    # J[-1, :] = J[0, :]
    # J[:,-1] = J[:,0]

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

    cx = U*dt/dx
    cy = V*dt/dy

    return cx, cy, u, v, X_edge, Y, J