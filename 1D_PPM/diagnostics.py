import numpy as np
import pylab as pl

# Various diagnostics of the outputs of linear advection

def errorNorms(phi, phiExact):
    "Calculates the l1, l2 and linfinity error norms of phi in comparison to"
    "phiExact"
    
    # input checking
    if (np.ndim(phi) != 1) or  (np.ndim(phiExact)) != 1:
        raise ValueError("In errorNorms(phi, phiExact), phi and phiExact must be one dimensional array")
    if (np.size(phi) <= 3) or (np.size(phiExact) !=  np.size(phi)):
        raise ValueError("In errorNorms(phi, phiExact), phi and phiExact must have the same number of elements and must have at least 3 elements (including 2 wrap-around)")

    #remove the wrap-around points
    phi = phi[0:-1]
    phiExact = phiExact[0:-1]
    
    # calculate the error and the error norms
    phiError = phi - phiExact
    l1 = sum(np.abs(phiError))/sum(np.abs(phiExact))
    # correct the following two lines
    l2 = np.sqrt(sum(phiError**2))/np.sqrt(sum(phiExact**2))
    linf = np.max(np.abs(phiError))

    return [l1,l2,linf]


def fourierPower(signal):
    "Calculates the Fourier power spectrum of the signal"
    
    # input checking
    if np.ndim(signal) != 1:
        raise ValueError("In fourierPower(signal), signal must be one dimensional array")

    A = pl.fft(signal)
    n = int((len(signal)+1)/2)
    power = pl.zeros(n)
    for i in range(n):
        power[i] = 4*(A[i+1].real**2 + A[i+1].imag**2)
    return power
   
def disperror(alpha,kdx, c):
    "demonstrate the dispersion error"
    
    #input checking
    if np.ndim(alpha) != 1:
        raise ValueError("In disperror(alpha,kdx, c), phi must be one dimensional array")
	
    # array for physical phase speed and computational phase speed   
    phyu = pl.zeros_like(kdx)
    comu = pl.zeros_like(kdx)
    # calculate the 
    phyu = alpha/(kdx*c)
    comu = (-np.pi+alpha)/(kdx*c)
    return phyu,comu

def orders():

	# create any array for error of different nx,dx and nt
    phierrors = pl.zeros((6,3,3),dtype='f')
    deltax = pl.zeros(3,dtype='f')
	
	#call function to get errors with respect to x
    phierrors[:,:,0],deltax[0] = advection(\
        initialProfile=cosBell, xmin = 0, xmax = 1, nx = 40, nt = 50, c = 0.2)
    phierrors[:,:,1], deltax[1]  = advection(\
        initialProfile=cosBell, xmin = 0, xmax = 1, nx = 100, nt = 125 , c = 0.2)
    phierrors[:,:,2],deltax[2] = advection(\
        initialProfile=cosBell, xmin = 0, xmax = 1, nx = 200, nt = 250, c = 0.2)
    
    #plot l1 versus dx
    pl.figure(5)
    pl.ion()
    pl.clf()
    pl.loglog(deltax, phierrors[0,0,:], 'r',label='FTBS')
    pl.loglog(deltax, phierrors[1,0,:], 'g',label='CTCS')
    pl.loglog(deltax, phierrors[2,0,:], 'k',label='diffusion')
    pl.loglog(deltax, phierrors[3,0,:], 'y',label='semi-langrangian')
    pl.loglog(deltax, phierrors[4,0,:], 'c',label='clip')
    pl.loglog(deltax, phierrors[5,0,:], 'm',label='Lax-Wendroff')
    pl.loglog([deltax[0],deltax[-1]], \
        [phierrors[0,0,0],phierrors[0,0,0]*(deltax[-1]/deltax[0])**2], 'b--', label='2nd order')
    pl.loglog([deltax[0],deltax[-1]], \
        [0.8*phierrors[0,0,0],0.8*phierrors[0,0,0]*(deltax[-1]/deltax[0])**1], 'r--', label='1st order')
    pl.axis([0.004,0.04,1e-3,2])
    pl.legend(loc='best')
    pl.xlabel('dx')
    pl.ylabel('l1 norm')
    
    #plot l2 versus dx
    pl.figure(6)
    pl.ion()
    pl.clf()
    pl.loglog(deltax, phierrors[0,1,:], 'r',label='FTBS')
    pl.loglog(deltax, phierrors[1,1,:], 'g',label='CTCS')
    pl.loglog(deltax, phierrors[2,1,:], 'k',label='diffusion')
    pl.loglog(deltax, phierrors[3,1,:], 'y',label='semi-langrangian')
    pl.loglog(deltax, phierrors[4,1,:], 'c',label='clip')
    pl.loglog(deltax, phierrors[5,1,:], 'm',label='Lax-Wendroff')
    pl.loglog([deltax[0],deltax[-1]], \
        [phierrors[0,1,0],phierrors[0,1,0]*(deltax[-1]/deltax[0])**2], 'b--',label='2nd order')
    pl.loglog([deltax[0],deltax[-1]], \
        [phierrors[0,1,0],phierrors[0,1,0]*(deltax[-1]/deltax[0])**1], 'r--',label='1st order')
    pl.axis([0.004,0.03,1e-5,1])
    pl.legend(loc='lower right')
    pl.xlabel('dx')
    pl.ylabel('l2 norm')