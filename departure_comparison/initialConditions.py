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

    def comput_SB(X, Y, ymax, fx, L):
        #-------------------------------------
        # non-orthogonal computational domain 
        # for solid body rotation test case
        #-------------------------------------
        return np.where(Y>0.5*L, fx+(Y-0.5*L)*(L-fx)/(0.5*L),Y*fx/(0.5*L))
    
    #-----------------------
    # Basic grid information
    #-----------------------
    dy,dx = Ly/ny,Lx/nx

    ymax = y[-1]
    A = np.pi*5/3000

    X,Y1 = np.meshgrid(x,y)

    #-----------------------------
    # streamfunction field
    #-----------------------------
    if mesh == 'quad':
        Y_psi = comput_SB(X,Y1,ymax,f_quad(X,Lx), Ly)
    elif mesh == 'V':
        Y = comput_SB(X,Y1,ymax,f_V(X,Lx), Ly)
    else:
        Y = Y1
    # streamfunction definition    
    r = np.sqrt((X-0.5*Lx)**2+(Y_psi-0.5*Ly)**2)

    psi = A*r**2
    psi[:,-1] = psi[:,0]

    #-----------------------------
    # dX and dY in physical domain
    #-----------------------------
    dX = np.zeros([ny+1,nx])
    dY = np.zeros([ny,nx+1])
    J =  np.zeros([ny+1, nx+1])
    dY = Y_psi[1:, :] - Y_psi[:-1, :]

    dX = X[:,1:] - X[:,:-1]

    J[:-1, :] = dy/dY
    J[-1, :] = J[0,:]
    #-----------------------------
    # Courant number and velocity
    #-----------------------------
    u = np.zeros([ny, nx+1])
    v = np.zeros([ny+1, nx])  
    U = np.zeros([ny, nx+1])
    V = np.zeros([ny+1, nx])

    #-----------------------------
    # velocity in physical domain
    #-----------------------------
    u = -(psi[1:,:]-psi[:-1,:])/dY
    v =  (psi[:,1:]-psi[:,:-1])/dX
  

    #-----------------------------
    # Courant number in computational domain
    #-----------------------------
    U = -(psi[1:,:]-psi[:-1,:])/dy
    V =  (psi[:,1:]-psi[:,:-1])/dx

    cx = U*dt/dx
    cy = V*dt/dy

    return cx, cy, u, v, X, Y_psi, J