import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
exec(open("initialConditions.py").read())
exec(open("2D_scheme.py").read())
exec(open("plot.py").read())
exec(open("diagnostics.py").read())
exec(open("streamfunction.py").read())
def advection(xmin, xmax,  ymin, ymax, epsilon, dt, nx, ny, nt, distorted = True):
    ''''Advect the given initial conditions function with domain with 
    nt time steps, Courant number c, where c is the combination of cx and cy'''
    Lx, Ly= xmax-xmin,ymax-ymin

    dx,dy = Lx/nx, Ly/ny
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
    phiOld,X,Y,J,psi = deformational(x,y,xmin,ymin,nx,ny,Lx,Ly,t,nt,dt,change = False)
    u,v,cx,cy = streamfunction(psi,dx,dy,dt,initialProfile)
    #calculate the divergence
    print np.max(phiOld),np.min(phiOld)
    dudx= (u[:-1,1:]-u[:-1,:-1])/dx
    dvdy = (v[1:,:-1]-v[:-1,:-1])/dy
    print 'max u', np.max(u),'max v',np.max(v)
    print 'the maximum divergence', np.max(dudx+dvdy)
    print 'the maximum Courant number is ',np.max(cx),np.max(cy),np.max(J[:-1,:-1]*dudx)*dt,np.max(J[:-1,:-1]*dvdy)*dt
    phiCOS,phiCOSMid,phiQ1,phiQ2 = COSMIC(phiOld,nt,epsilon,dx,dy,J,initialProfile,xmin,ymin,dt,UchangewithT,cx,cy)
    np.savez(str(nt), phiCOS = phiCOS, phiCOSMid = phiCOSMid, phiQ1 = phiQ1, phiQ2 = phiQ2)

    f = np.load(str(nt)+'.npz')
    print f.files
    phiCOS = f['phiCOS']
    phiCOSMid = f['phiCOSMid']
    phiQ1 = f['phiQ1']
    phiQ2 = f['phiQ2']
    print 'the max and min at Midtime',np.max(phiCOSMid),np.min(phiCOSMid)
    print 'the max and min at final time',np.max(phiCOS),np.min(phiCOS)
    l1,l2,linf = norms(phiOld, phiCOS, dx, dy)
    print 'the l1 norm at final time is ', l1, 'the l2 norm is ', l2,'the linf norm is ', linf
    contours_blob(X,Y,phiCOS,u,v,x,nt)
    contours_blob_mid(X,Y,phiCOSMid,nt)
    contours_blob_q(X,Y,phiQ1,nt)
    contours_blob_q2(X,Y,phiQ2,nt)
advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.04, nx = 60, ny = 30, nt = 125)
advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.02, nx = 120, ny = 60, nt = 250)
advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt= 0.01, nx = 240, ny = 120, nt =500)
advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.02/3, nx = 360, ny = 180, nt = 750)
advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt= 0.005, nx = 480, ny = 240, nt =1000)