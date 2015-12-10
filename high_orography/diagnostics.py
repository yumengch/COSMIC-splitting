import numpy as np
import matplotlib.pyplot as plt
def analytical(X,Y,dx,dy,nt,dt,nx,ny,u0=10):
    x0 = -50e3+nt*dt*u0
    Ax = 25e3
    Az = 3e3
    z0 = 12e3
    r = np.sqrt(((X-x0)/Ax)**2 + ((Y-z0)/Az)**2)
    phi = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 
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