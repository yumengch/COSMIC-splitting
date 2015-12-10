import numpy as np
import matplotlib.pyplot as plt
def solid_ana(X,Y,dx,dy,nt,dt,nx,ny):
    phi = np.zeros((ny+1,nx+1))
    A= np.pi*5/3/1000.
    # r = ((X-0.5*dx)-(0.5*nx)*dx)**2+((Y-0.5*dy)-(0.5*ny)*dy)**2
    r = 2500.
    theta0 = np.pi/2.
    x0 = 0.5*nx*dx+r*np.cos(theta0+2*A*nt*dt)
    y0 = 0.5*ny*dy+ r*np.sin(theta0+2*A*nt*dt)
    # u0 = 10.
    # x0 = u0*nt*dt
    # y0 = 0.5*ny*dy+r
    phi = np.exp(- 0.5*(((X-x0) / (500))**2 + ((Y-y0) / (500))**2))
    return phi

def norms(phiExact, phi, dx, dy):
	#remove the wrap-around points
    phi = phi[0:-1,0:-1]
    phiExact = phiExact[0:-1,0:-1]
    
    # calculate the error and the error norms
    phiError = phi - phiExact
    l1 = np.sum(np.abs(phiError))/np.sum(np.abs(phiExact))
    # correct the following two lines
    l2 = np.sqrt(np.sum(phiError**2))/np.sqrt(np.sum(phiExact**2))
    linf = np.max(np.abs(phiError))

    return [l1,l2,linf]