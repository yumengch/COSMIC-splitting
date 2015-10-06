import numpy as np
import matplotlib.pyplot as plt

xmin = -150e3
xmax = 150e3
zmin = 0
zmax = 25e3

dx = 1e3
dz = 500
xmin -= (dx/2.0)
xmin += (dx/2.0)
nx = int((xmax-xmin)/dx)
nz = int((zmax-zmin)/dz)
tEnd = 10e3
dt = 25
nt = int(tEnd/dt)

x = np.linspace(xmin,xmax,nx+1)
z = np.linspace(zmin,zmax,nz+1)
X,Z1 = np.meshgrid(x,z)
rho = np.zeros((nz+1, nx+1))
Psi= np.zeros((nz+1,nx+1))
def h(x):
    "The mountain height as a function of x"
    h0 = 3e3
    a = 25e3
    lam = 8e3
    # x = x - 500
    return np.where(abs(x) <= a,
                    h0*np.cos(np.pi*x/lam)**2 * np.cos(np.pi*x*0.5/a)**2,
                    0)

def sigmaHeights(x,Z,zmax):
    "The height of sigma coordinate Z"
    hx = h(x)
    return hx + Z/zmax*(zmax - hx)

Ax = 25e3
Az = 3e3
x0 = -50e3
z0 = 9e3
Z = sigmaHeights(X,Z1,zmax)
r = np.sqrt((((X)-x0)/Ax)**2 + ((Z-z0)/Az)**2)
rho = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0)

J = zmax/(zmax - h(X-0.5*dx))
J[:,0] = J[:,-1]
J[:,:-1] = (J[:,1:] + J[:,:-1])/2
J[:,-1] = J[:,0]

u0 = 10.
z1 = 4e3
z2 = 5e3
Z = sigmaHeights(X-0.5*dx,Z1-0.5*dz,zmax)
for i in xrange(nx+1):
    for j in xrange(nz+1):   
        if (Z[j,i] <= z1):
            Psi[j,i] = 0
        elif (Z[j,i] <= z2):
            Psi[j,i] = -0.5*u0*(Z[j,i]-z1-(z2-z1)/np.pi*np.sin(np.pi*((Z[j,i])-z1)/(z2-z1)))
        else:
            Psi[j,i] = -0.5*u0*(2*(Z[j,i]) - z1 - z2)
# Psi[:,0] = Psi[:,-1]
Psi[:,-1] = Psi[:,0]
rhoNew = rho.copy()
rhoOld = rho + 0.5*dt*u0/dx\
         *(np.roll(rho,-1,axis=-1) - np.roll(rho,1,axis=-1))

# Time steps to store
rhoInit = rho.copy()
rhoMid = rho.copy()

for it in range(0,nt):
    # loop over space (away from the boundaries, no change on boundaries)
    for i in range(1,nx):
        for k in range(1,nz):
            rhoNew[k,i] = rhoOld[k,i] - (J[k,i])*dt/(dx*dz)*\
            (\
                - (Psi[k+1,i+1]-Psi[k,i+1])*rho[k,i+1] \
                + (Psi[k+1,i]-Psi[k,i])*rho[k,i-1] \
                + (Psi[k+1,i+1]-Psi[k+1,i])*rho[k+1,i] \
                - (Psi[k,i+1]-Psi[k,i])*rho[k-1,i] \
            )
    # update old arrays
    rhoOld = rho.copy()
    rho = rhoNew.copy()
    print rho[:-1,:-1].sum()
    # store mid value if we have arrived there
    if it*dt == tEnd/2:
        rhoMid = rho.copy()
# Final plot with start, end and mid time-steps
x0 = 50e3
r = np.sqrt(((X-x0)/Ax)**2 + ((Z-z0)/Az)**2)
phiExact = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 
x0 = 0
r = np.sqrt(((X-x0)/Ax)**2 + ((Z-z0)/Az)**2)
phiMidExact = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 
print np.amax(rho),np.amin(rho)
phiError = rho[0:-1] - phiExact[0:-1]
l2 = np.sqrt(np.sum(phiError**2)/np.sum(phiExact**2))
l1 = np.sum(np.abs(phiError))/np.sum(np.abs(phiExact))
linf = np.max(np.abs(phiError))
print l1,l2,linf



plt.figure(1)
plt.clf()
ax = plt.subplot(111)
plt.contour(X/1000, Z/1000, np.round(rho,3),
           np.concatenate((np.arange(-2,0,0.1), np.arange(0.1,2,0.1))),
          colors = 'k')
plt.contour(X/1000,Z/1000, np.round(rhoMid,3),
           np.concatenate((np.arange(-2,0,0.1), np.arange(0.1,2,0.1))),
          colors = 'k', hold='on')
plt.contour(X/1000, Z/1000, np.round(rhoInit,3),
           np.concatenate((np.arange(-2,0,0.1), np.arange(0.1,2,0.1))),
          colors = 'k', hold='on')
plt.plot(x/1000,h(x)/1000,'k.')
phiError = np.round(rho - phiExact,3)
phiMidError= np.round(rhoMid-phiMidExact,3)
print np.amax(phiMidError),np.amin(phiMidError)
color = ['#400056','#2B006B','#150081','#0000FF','#0031FF','#0061FF','#0092FF','#00C2FF','#00F3FF','#6DFFFF','#FFFFFF','#FFFF6D','#FFF300',\
    '#FFC200','#FF9200','#FF6100','#FF3100','#FF0000','#CA0600','#950B00','#601100']
# levs =np.concatenate([np.arange(-0.4,-0.02,0.02), np.arange(0.02,0.4,0.02)])
levs = np.arange(-0.21,0.21,0.02)
# define the bins and normalize and forcing 0 to be part of the colorbar!
# plt.contourf(X/1000,Z/1000,np.round(phiMidError,3),levels =levs,colors = color)
cs = plt.contourf(X/1000,Z/1000,np.round(phiError+phiMidError,3),levels = levs,colors = color,extend="both")
cs.cmap.set_under('#400056')
cs.cmap.set_over('#601100')
cs.set_clim(-0.2, 0.2)
cb = plt.colorbar(cs)
cb.set_ticks(np.round(np.arange(-0.16,0.24,0.08),2))
cb.set_ticklabels(['%1.2f' % i for i in np.arange(-0.16,0.24,0.08)])
plt.figtext(.3, 0.07, "min = "+str(np.round(np.min(phiError),3))+', max = '+str(np.round(np.max(phiError),3)))
plt.xlim([-75,75])
plt.ylim([0,15])
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
# plt.savefig('CTCS.pdf')
