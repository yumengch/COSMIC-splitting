import numpy as np
def analytical(X,Y,dx,dy,nt,dt,nx,ny):
	A = 8
	r = np.sqrt(((75-0.5*nx)*dx)**2+((50-0.5*ny)*dy)**2)
	theta0 = np.arctan(((50-0.5*ny)*dy)/((75-0.5*nx)*dx))
	x0 = 0.5*nx*dx+r*np.cos(theta0+A*nt*dt)
	y0 =.5*ny*dy+ r*np.sin(theta0+A*nt*dt)
	print x0,y0
	phi = np.exp(- 0.5*(((X-x0) / (3*dx))**2 + ((Y-y0) / (3*dy))**2))
	return phi

def norms(phiExact, phi, dx, dy):
	#remove the wrap-around points
    phi = phi[0:-1]
    phiExact = phiExact[0:-1]
    
    # calculate the error and the error norms
    phiError = phi - phiExact
    l1 = np.sum(np.abs(phiError))/np.sum(np.abs(phiExact))
    # correct the following two lines
    l2 = np.sqrt(np.sum(phiError**2))/np.sqrt(np.sum(phiExact**2))
    linf = np.max(np.abs(phiError))

    return [l1,l2,linf]