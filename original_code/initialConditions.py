'''this file is used to implement the initial conditions to test advection schemes. The scheme uses Arakawa C 
grid. The streamfunction thus must shift half of the grid points both in x and y direction'''
import numpy as np
import matplotlib.pyplot as plt
def Smolarkiewicz(x,y,xmin,ymin,nx,ny,Lx,Ly,t,nt,dt):
    psi= np.zeros((ny+1,nx+1))
    A, k = 8, 4*np.pi/Lx
    dx,dy = Lx/nx, Ly/ny
    X,Y = np.meshgrid(x,y)
    if t <0:
        #define the initial cone phi and stream function
        J = np.ones((ny+1,nx+1))
        psi = A*np.sin(k*(X-0.5*dx))*np.cos(k*(Y-0.5*dy))
        cone = lambda x,y: (1 - np.sqrt((X - dx*((nx/2)))**2 + (Y - dy*((ny/2)))**2)/(15))
        phi = np.where((((X - dx*((nx/2)))**2 + (Y - dy*((ny/2)))**2) <= ((15)**2)), cone(X,Y),0 
              )
        UchangewithT = False
        return phi,X,Y,J,psi,UchangewithT
    else:
        #define the streamfunction
        psi = Lx*(Y-0.5*dy)/(nt*dt)+np.cos(np.pi*t/nt)*A*np.sin(k*(X-0.5*dx+(Lx*t/nt)))*np.cos(k*(Y-0.5*dy))
        psi[0,:] = psi[-1,:]
        psi[:,0] = psi[:,-1]
        return psi

def constant(x,y,xmin,ymin,nx,ny,Lx,Ly,t,nt,dt):
    psi = np.zeros((ny+1,nx+1))
    J = np.ones((ny+1,nx+1))
    X,Y = np.meshgrid(x,y)
    phi = np.zeros((ny+1,nx+1))
    dx,dy = Lx/nx, Ly/ny
    u0 = 1.
    psi = -(u0*Ly/(1*np.pi))*np.sin(2*np.pi*((X-0.5*dx)-xmin)/Lx)*np.sin(1*np.pi*((Y-0.5*dy)-ymin)/Ly)
    psi[0,:] = psi[-1,:]
    psi[:,0] = psi[:,-1]
    phi[:,:] = 1
    UchangewithT = False
    return phi,X,Y,J,psi,UchangewithT

def stirring(x,y,xmin,ymin,nx,ny,Lx,Ly,t,nt,dt):
    psi= np.zeros((ny+1,nx+1))
    X,Y = np.meshgrid(x,y)
    dx,dy = Lx/nx, Ly/ny
    u0 = 1.
    if t <0:
        J = np.ones((ny+1,nx+1))
        phi = np.zeros((ny+1,nx+1))
        rmax = 3*Lx/8.
        r = np.sqrt((X-dx*(nx/2))**2+(Y-dy*(ny/2))**2)
        # psi = -(u0*Ly/(1*np.pi))*np.sin(2*np.pi*((X-0.5*dx)-xmin)/Lx)*np.sin(1*np.pi*((Y-0.5*dy)-ymin)/Ly)
        psi = u0*(np.sin(X/Lx))**2*(np.cos(Y/Ly))**2
        phi[:,:] = 0.5*(1+np.cos(np.pi*(np.minimum(r,rmax)/rmax)))
        UchangewithT = True
        return phi,X,Y,J,psi,UchangewithT
    else:
        psi = -Lx*(Y-0.5*dy)/(nt*dt)-np.cos(np.pi*t/nt)*(u0*Ly/(1*np.pi))*np.sin(2*np.pi*((X-0.5*dx)-xmin-(Lx*t/nt))/Lx)*np.sin(1*np.pi*((Y-0.5*dy)-ymin)/Ly)
        psi[0,:] = psi[-1,:]
        psi[:,0] = psi[:,-1]
    return psi

def solid(x,y,xmin,ymin,nx,ny,Lx,Ly,t,nt,dt):
    psi= np.zeros((ny+1,nx+1))
    J = np.ones((ny+1,nx+1))
    X,Y = np.meshgrid(x,y)
    phi = np.zeros((ny+1,nx+1))
    dx,dy = Lx/nx, Ly/ny
    A= 8
    r = ((X-0.5*dx)-(0.5*nx)*dx)**2+((Y-0.5*dy)-(0.5*ny)*dy)**2
    psi = A*r
    psi[0,:] = psi[-1,:]
    psi[:,0] = psi[:,-1]
    r = np.sqrt((30)**2+(30)**2)
    theta0 = np.arctan((30)/(30))
    x0 = 0.5*nx*dx+r*np.cos(theta0)
    y0 = 0.5*ny*dy+ r*np.sin(theta0)
    phi = np.exp(- 0.5*(((X-x0) / (3))**2 + ((Y-y0) / (3))**2))
    UchangewithT = False
    return phi,X,Y,J,psi,UchangewithT 

def terrain(x,z,xmin,zmin,nx,nz,Lx,Lz,t,nt,dt):
    def h(x):
        "The mountain height as a function of x"
        h0 = 3e3
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
    z1 = 4e3
    z2 = 5e3
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
    z0 = 9e3
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
    UchangewithT = False
    return phi,X,Z,J,psi,UchangewithT