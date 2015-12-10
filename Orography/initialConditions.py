'''this file is used to implement the initial conditions to test advection schemes. The scheme uses Arakawa C 
grid. The streamfunction thus must shift half of the grid points both in x and y direction'''
import numpy as np
import matplotlib.pyplot as plt
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
        return zmax*(Z-hx)/(zmax-hx)

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
                psi[j,i] = 0.
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
    plt.figure(1)
    plt.contour(X,Z1,Z,levels = np.arange(0,zmax,1000),color = 'k')
    plt.show()
    return phi,X,Z,J,psi,UchangewithT