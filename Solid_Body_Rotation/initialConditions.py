def solid(x,y,xmin,ymin,nx,ny,Lx,Ly,t,nt,dt):
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

    dy,dx = Ly/ny,Lx/nx
    X,Y1 = np.meshgrid(x,y)
    phi = np.zeros((ny+1, nx+1))
    ymax = ymin+Ly
    # 1d arrays (x and y locations for different variables
    # Y = Y1
    Y = computational(X,Y1,ymax)
    psi= np.zeros((ny+1,nx+1))
    # J = 1
    J= np.where(Y1>0.5*Ly, 0.5*Ly/(Ly-f(X,Lx)),0.5*Ly/f(X,Lx))
    A= np.pi*5/3/1000.
    r = 2500
    theta0 = np.pi/2.
    x0 = 0.5*nx*dx+r*np.cos(theta0)
    y0 = 0.5*ny*dy+ r*np.sin(theta0)
    phi = np.exp(- 0.5*(((X-x0)/500)**2 + ((Y-y0)/500)**2))
    Y = computational(X-0.5*dx,Y1-0.5*dy,ymax)
    r = ((X-0.5*dx)-(0.5*nx)*dx)**2+((Y-0.5*dy)-0.5*Ly)**2
    # rmax = 8600.**2
    # psi = np.where(r<=rmax,A*r, A*rmax)
    psi = A*r

    UchangewithT = False
    Y = computational(X,Y1,ymax)
    return phi,X,Y,J,psi,UchangewithT,Y1