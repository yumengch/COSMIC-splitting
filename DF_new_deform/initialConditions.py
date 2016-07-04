'''this file is used to implement the initial conditions to test advection schemes. The scheme uses Arakawa C 
grid. The streamfunction thus must shift half of the grid points both in x and y direction'''
import numpy as np
import matplotlib.pyplot as plt
def deformational(x,y,xmin,ymin,nx,ny,Lx,Ly,nt,dt, distorted, t, change):
    def f(x,Lx):
        "the function that distort the grid"
        return np.where(x<=0.25*Lx,
                        -(1./np.sqrt(3))*(x)+Lx/8./np.sqrt(3),
                        np.where((x<0.5*Lx)&(x>=0.25*Lx),
                        1./np.sqrt(3)*(x-0.25*Lx)-Lx/8./np.sqrt(3),
                        np.where((x<0.75*Lx)&(x>=0.5*Lx),
                        -(1./np.sqrt(3))*(x-0.5*Lx)+Lx/8./np.sqrt(3),1./np.sqrt(3)*(x-0.75*Lx)-Lx/8./np.sqrt(3))))
        
    def computational(X,Y,ymax):
        "The computational grids"
        fx = f(X,Lx)
        return np.where(Y>=0, fx+Y*(1-fx/ymax),fx+Y*(1-fx/ymin))
    psi = np.zeros((ny+1,nx+1))
    phi = np.zeros((ny+1,nx+1))
    X,Y1 = np.meshgrid(x,y)
    dx,dy = Lx/nx, Ly/ny
    ymax = ymin+Ly
    u0 = 10.
    if distorted:
        Y = computational(X,Y1,ymax)
        J = np.where(Y1>=0, ymax/(ymax-f(x,Lx)),ymin/(ymin-f(x,Lx)))
        x0 = 5.*Lx/12.
        y0 = 0
        x1 = 7.*Lx/12.
        Y = computational(X,Y1,ymax)
        phi = 0.95*np.exp(- 5*((X-x0)**2 + (Y-y0)**2)) + 0.95*np.exp(- 5*((X-x1)**2+ (Y-y0)**2))
        Y = computational(X-0.5*dx,Y1-0.5*dy,ymax)       
        psi = (u0/nt/dt)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*(((X-0.5*dx)/Lx) - (float(t)/nt))))**2)*((np.cos(np.pi*(Y-0.5*dy)/Lx))**2)      
        Y = computational(x,Y1,ymax)
        
        if change:
            Y = computational(X-0.5*dx,Y1-0.5*dy,ymax)
            psi = (u0/nt/dt)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*(((X-0.5*dx)/Lx) - (float(t)/nt))))**2)*((np.cos(Y-0.5*dy))**2)*np.cos(np.pi*t/nt)-Lx*(Y-0.5*dy)/(nt*dt)
            return psi
        else:
            return phi,X,Y,J,psi,Y1
    else:
        J = np.ones((ny+1,nx+1))
        Y = Y1
        x0 = 5.*Lx/12.
        y0 = 0
        x1 = 7.*Lx/12.
        phi = 0.95*np.exp(- 5*((X-x0)**2 + (Y-y0)**2)) + 0.95*np.exp(- 5*((X-x1)**2+ (Y-y0)**2))
        # psi = u0*((np.sin(2*np.pi*((X-0.5*dx)/Lx))**2))*((np.cos(np.pi*(Y-0.5*dy)/Lx))**2)
        psi = (u0/nt/dt)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*((X-0.5*dx)/Lx - float(t)/nt)))**2)*((np.cos(np.pi*(Y-0.5*dy)/Ly))**2)*np.cos(np.pi*t/nt)-Lx*(Y-0.5*dy)/(nt*dt)
        if change:
            psi = (u0/nt/dt)*((Lx/2/np.pi)**2)*((np.sin(2*np.pi*((X-0.5*dx)/Lx - float(t)/nt)))**2)*((np.cos(np.pi*(Y-0.5*dy)/Ly))**2)*np.cos(np.pi*t/nt)-Lx*(Y-0.5*dy)/(nt*dt)
            return psi
        else:
            return phi,X,Y,J,psi
        