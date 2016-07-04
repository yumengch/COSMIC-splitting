# import numpy as np
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
    print Lx
    dx,dy = Lx/nx, Ly/ny
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
    phiOld,X,Y,J,psi,Y1 = deformational(x,y,xmin,ymin,nx,ny,Lx,Ly,nt,dt, distorted,t = 0,change = False)
    u,v,cx,cy = streamfunction(psi,dx,dy,dt,J)
    #calculate the divergence
    print np.max(phiOld),np.min(phiOld)
    dudx= (u[:-1,1:]-u[:-1,:-1])/dx
    dvdy = (v[1:,:-1]-v[:-1,:-1])/dy
    # print Y[1:] - Y[:-1]
    print 'max u', np.max(u),'max v',np.max(v)
    print 'the maximum divergence', np.max(dudx+dvdy)
    print 'the maximum Courant number is ',np.max(cx),np.max(cy),np.max(dudx*(J[:-1,:-1]))*dt,np.max(dvdy)*dt
    # plt.figure(1)
    # plt.clf()
    # plt.quiver(X[::3],Y1[::3],u[::3],v[::3])
    # plt.contour(X,Y,np.round(phiOld,3),np.arange(0.1,1.0,0.1),colors = 'k')
    # plt.contour(X,Y,X,np.arange(0,2*np.pi,np.pi/10),colors = "k")
    # plt.contour(X,Y,Y1,np.arange(-0.5*np.pi,0.5*np.pi,np.pi/10),colors = 'k')
    # plt.show()
    # phiT1,phiT2,phiT3,phiT4,phiCOS = COSMIC(phiOld,nt,epsilon,dx,dy,J,xmin,ymin,dt,J*cx,J*cy,distorted)
    # np.savez('nx'+str(nx)+'nt'+str(nt), phiT1 = phiT1, phiT2 = phiT2, phiT3 = phiT3, phiT4 = phiT4, phiCOS = phiCOS)
    f = np.load('nx'+str(nx)+'nt'+str(nt)+'.npz')
    print f.files
    phiCOS = f['phiCOS']
    phiT1 = f['phiT1']
    phiT2 = f['phiT2']
    phiT3 = f['phiT3']
    phiT4 = f['phiT4']
    # print 'the max and min at Midtime',np.max(phiCOSMid),np.min(phiCOSMid)
    print 'the max and min at final time',np.max(phiCOS),np.min(phiCOS)
    l1,l2,linf = norms(phiOld, phiCOS, dx, dy)
    print 'the l1 norm at final time is ', l1, 'the l2 norm is ', l2,'the linf norm is ', linf
    # contours_Old(X*180./np.pi - 180.,Y*180./np.pi,phiOld,nt)
    # contours_T4(X*180/np.pi - 180.,Y*180/np.pi,phiT4,nt)
    # contours_T3(X*180/np.pi - 180.,Y*180/np.pi,phiT3,nt)
    # contours_T2(X*180/np.pi - 180.,Y*180/np.pi,phiT2,nt)
    # contours_TT1(X*180/np.pi - 180.,Y*180/np.pi,phiT1,nt)
    # contours_END(X*180/np.pi - 180.,Y*180/np.pi,phiCOS,u,v,x,nt)

# # # dx, small c
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.02, nx = 60, ny = 30, nt = 250)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.01, nx = 120, ny = 60, nt = 500)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.005, nx = 240, ny = 120, nt = 1000)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.01/3, nx = 360, ny = 180, nt = 1500)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.0025, nx = 480, ny = 240, nt = 2000)

# dx, large c

# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.2, nx = 60, ny = 30, nt = 25)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.1, nx = 120, ny = 60, nt = 50)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.05, nx = 240, ny = 120, nt = 100)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.1/3, nx = 360, ny = 180, nt = 150)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.025, nx = 480, ny = 240, nt = 200)


# dt
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.1, nx = 120, ny = 60, nt = 50)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.05, nx = 120, ny = 60, nt = 100)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.02, nx = 120, ny = 60, nt = 250)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.01, nx = 120, ny = 60, nt = 500)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.005, nx = 120, ny = 60, nt = 1000)
# advection(xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, epsilon = 0.01, dt = 0.0025, nx = 120, ny = 60, nt = 2000)