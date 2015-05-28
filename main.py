import numpy as np

from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
exec(open("initialConditions.py").read())
exec(open("2D_scheme.py").read())
exec(open("plot.py").read())
exec(open("diagnostics.py").read())
def advection(initialProfile, xmin = 0.0, xmax = 1.0,  ymin = 0.0, ymax = 1.0, epsilon = 0.01, dt= 200, nx = 200, ny = 200, nt =500):
    ''''Advect the given initial conditions function with domain with 
    nt time steps, Courant number c, where c is the combination of cx and cy'''
    #np.seterr(all='raise')
    u,v = np.zeros((2,ny+1,nx+1))
    Lx, Ly= xmax-xmin,ymax-ymin
    dx,dy = Lx/nx, Ly/ny
    x,y = np.linspace(0,Lx,nx+1), np.linspace(0,Ly,ny+1)
    phiOld,psi,X,Y = initialProfile(x,y,xmin,ymin,nx,ny,Lx,Ly)
    #define the velocity field
    u[:-1,:] = -(psi[1:,:]-psi[:-1,:])/dy
    v[:,:-1] =  (psi[:,1:]-psi[:,:-1])/dx
    v[:,-1],v[-1,:] = v[:,0],v[0,:]
    u[-1,:],u[:,-1] = u[0,:],u[:,0]
    #calculate the divergence
    dudx= (u[:-1,1:]-u[:-1,:-1])/dx
    dvdy = (v[1:,:-1]-v[:-1,:-1])/dy
    #calculate the Courant number
    cx = u*dt/dx
    cy = v*dt/dy
    print 'the maximum and minimum Courant number is ',np.max(cx),np.max(cy),np.max(dudx)*dt,np.max(dvdy)*dt
    phiCOS = COSMIC(phiOld,cx,cy,nt,epsilon,dx,dy)
    # phiExact = analytical(X,Y,dx,dy,nt,dt,nx,ny)
    # phirev = COSMIC(phiCOS,-cx,-cy,nt,epsilon,dx,dy)
    # l1,l2,linf = norms(phiExact, phiCOS, dx, dy)
    # print 'the l1 norm is ', l1, 'the l2 norm is ', l2,'the linf norm is ', linf
    # # return phiCOS
    hist(X,Y,phiCOS)
    # contours(X,Y,phiCOS)


# advection(initialProfile = blob, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 1.5, nx = 100, ny = 100, nt =60)
# advection(initialProfile = solid, xmin = 0., xmax = 100., ymin = 0., ymax = 100., epsilon = 0.01, dt= 100, nx = 100, ny = 100, nt =30)
# advection(initialProfile = Smolarkiewicz, xmin = 0., xmax = 100., ymin = 0., ymax = 100., epsilon = 0.01, dt= 1.4, nx = 200, ny = 200, nt =32)
advection(initialProfile = constant, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 1.5, nx = 50, ny = 50, nt =20)
