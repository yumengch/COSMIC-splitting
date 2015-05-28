'''this file is used to implement the initial conditions to test advection schemes. The scheme uses Arakawa C 
grid. The streamfunction thus must shift half of the grid points both in x and y direction'''
import numpy as np
def Smolarkiewicz(x,y,xmin,ymin,nx,ny,Lx,Ly):
    psi= np.zeros((ny+1,nx+1))
    A, k = 8, 4*np.pi/Lx
    dx,dy = Lx/nx, Ly/ny
    X,Y = np.meshgrid(x,y)
    #define the streamfunction
    psi = A*np.sin(k*(X-0.5*dx))*np.cos(k*(Y-0.5*dy))
    psi[0,:] = psi[-1,:]
    psi[:,0] = psi[:,-1]
    #define the initial cone phi condition
    cone = lambda x,y: (1 - np.sqrt((X - dx*((nx/2)))**2 + (Y - dy*((ny/2)))**2)/(15))
    phi = np.where((((X - dx*((nx/2)))**2 + (Y - dy*((ny/2)))**2) <= ((15)**2)), cone(X,Y),0 
          )
    return phi,psi,X,Y

def constant(x,y,xmin,ymin,nx,ny,Lx,Ly):
    psi = np.zeros((ny+1,nx+1))
    X,Y = np.meshgrid(x,y)
    phi = np.zeros((ny+1,nx+1))
    dx,dy = Lx/nx, Ly/ny
    u0 = 1.
    psi = -(u0*Ly/(ny*np.pi))*np.sin(2*np.pi*((X-0.5*dx)-xmin)/Lx)*np.sin(1*np.pi*((Y-0.5*dy)-ymin)/Ly)
    psi[0,:] = psi[-1,:]
    psi[:,0] = psi[:,-1]
    #define the velocity field
    phi[:,:] = 1.0
    return phi,psi,X,Y

def blob(x,y,xmin,ymin,nx,ny,Lx,Ly):
    psi= np.zeros((ny+1,nx+1))
    X,Y = np.meshgrid(x,y)
    phi = np.zeros((ny+1,nx+1))
    rmax = np.zeros((ny+1,nx+1))
    dx,dy = Lx/nx, Ly/ny
    rmax[:,:] = 3*Lx/4.
    u0 = 1.
    psi = -(u0*Ly/(ny*np.pi))*np.sin(2*np.pi*((X-0.5*dx)-xmin)/Lx)*np.sin(1*np.pi*((Y-0.5*dy)-ymin)/Ly)
    psi[0,:] = psi[-1,:]
    psi[:,0] = psi[:,-1]
    r = np.sqrt((X-dx*(nx/2))**2+(Y-dy*(ny/2))**2)
    phi[:,:] = 0.5*(1+np.cos(np.pi*(np.minimum(r,rmax))/rmax))
    return phi,psi,X,Y

def solid(x,y,xmin,ymin,nx,ny,Lx,Ly):
    psi= np.zeros((ny+1,nx+1))
    X,Y = np.meshgrid(x,y)
    phi = np.zeros((ny+1,nx+1))
    dx,dy = Lx/nx, Ly/ny
    A= 8
    psi = A*(0.5*((X-0.5*dx)**2+(Y-0.5*dy)**2)-((nx/2.)*dx*(X-0.5*dx)+(ny/2)*dy*(Y-0.5*dy)))
    psi[0,:] = psi[-1,:]
    psi[:,0] = psi[:,-1]
    phi = np.exp(- 0.5*(((X-75*dx) / (3*dx))**2 + ((Y-50*dy) / (3*dy))**2))
    # phi = np.exp(-4*np.log(2) * ((X-75*dx)**2 + (Y-75*dx)**2) / 3**2)
    return phi,psi,X,Y