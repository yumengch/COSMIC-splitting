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

    dX = np.zeros([ny+1, nx+1])
    dY = np.zeros([ny+1, nx+1])


    dX[:,:-1] = X[:, 1:] - X[:, :-1]
    dX[:,-1] = dX[:,0]
    dY[:-1, :] = Y[1:, :] - Y[:-1, :]
    dY[-1,:] = dY[-2,:]
    # dX[:,1:-1] = 0.5*(X[:, 1:-1] - X[:, :-2]) + 0.5*(X[:, 2:] - X[:, 1:-1])
    # dX[:, 0] = 0.5*(X[:, -1] - X[:, -2]) + 0.5*(X[:, 1] - X[:, 0])
    # dX[:,-1] = dX[:,0]
    # # print dX
    # dY[1:-1, :] = 0.5*(Y[1:-1, :] - Y[:-2, :]) + 0.5*(Y[2:, :] - Y[1:-1, :])
    # dY[0,:] = 0.5*(Y[-1,:] - Y[-2,:]) + 0.5*(Y[1, :] - Y[0, :])
    # dY[-1,:] = dY[0,:]
    J = dy/dY
    # print dX
    #-----------------------------
    # Courant number
    #-----------------------------
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)    
    U = np.zeros_like(psi)
    V = np.zeros_like(psi)

    # the velocity field
    u[:-1,:] = -(psi[1:,:]-psi[:-1,:])/dY[:-1,:]
    v[:,:-1] =  (psi[:,1:]-psi[:,:-1])/dX[:,:-1]
    u[-1,:] = u[0,:]
    v[:,-1] = v[:,0]

    U[:-1,:] = -(psi[1:,:]-psi[:-1,:])/dy
    V[:,:-1] =  (psi[:,1:]-psi[:,:-1])/dx
    V[:,-1] = V[:,0]
    U[-1,:] = U[0,:]

    # print dY
    # print np.max(u - U), np.max(v-V)
    # print u
    # print U
    cx = U*dt/dx
    cy = V*dt/dy

    return phi, cx, cy, u, v, X, Y, J