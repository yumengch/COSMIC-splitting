import numpy as np
import time
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
exec(open("initialConditions.py").read())
exec(open("2D_scheme.py").read())
exec(open("plot.py").read())
exec(open("diagnostics.py").read())
exec(open("streamfunction.py").read())
def advection(initialProfile, xmin, xmax,  ymin, ymax, epsilon, dt, nx, ny, nt):
    ''''Advect the given initial conditions function with domain with 
    nt time steps, Courant number c, where c is the combination of cx and cy'''
    #np.seterr(all='raise')
    start_time = time.time()
    Lx, Ly= xmax-xmin,ymax-ymin
    dx,dy = Lx/nx, Ly/ny
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)+0.5*dy
    t = -1
    phiOld,X,Y,J,psi,UchangewithT = initialProfile(x,y,xmin,ymin,nx,ny,Lx,Ly,t,nt,dt)
    u,v,cx,cy = streamfunction(psi,dx,dy,dt,initialProfile)
    #calculate the divergence
    print np.max(phiOld),np.min(phiOld)
    dudx= (u[:-1,1:]-u[:-1,:-1])/dx
    dvdy = (v[1:,:-1]-v[:-1,:-1])/dy
    # print 'the maximum and minimum Courant number is ',np.max(J*cx),np.max(J*cy),np.max(J[:-1,:-1]*dudx)*dt,np.max(J[:-1,:-1]*dvdy)*dt
    print 'the maximum and minimum Courant number is ',np.max(cx),np.max(cy),np.max(dudx)*dt,np.max(dvdy)*dt
    # # # print np.max(J[:-1,:-1]*(dudx+dvdy))
    phiCOS,phiCOSMid,phiquater,phiquater2 = COSMIC(phiOld,nt,epsilon,dx,dy,J,initialProfile,xmin,ymin,dt,UchangewithT,cx,cy)
    print np.max(phiCOS),np.min(phiCOS)
    print np.max(phiCOSMid),np.min(phiCOSMid)
    print phiCOS[:-1,:-1].sum()-phiOld[:-1,:-1].sum()
    print np.max(phiCOSMid),np.min(phiCOSMid)
    phiExact = analytical(X,Y,dx,dy,nt,dt,nx,ny)
    phiMidExact = analytical(X,Y,dx,dy,int(nt/2.)+1,dt,nx,ny)
    l1,l2,linf = norms(phiExact, phiCOS, dx, dy)
    l1Mid,l2Mid,linfMid = norms(phiMidExact, phiCOSMid, dx, dy)
    print 'the l1 norm is ', l1, 'the l2 norm is ', l2,'the linf norm is ', linf
    print 'the l1 norm is ', l1Mid, 'the l2 norm is ', l2Mid,'the linf norm is ', linfMid
    # init(X,Y,u,v,phiOld,x)
    contourconst(X,Y,phiCOS)
    # contours_blob(X,Y,phiCOS,u,v,x,nt)
    # contours_blob_mid(X,Y,phiCOSMid,nt)
    # contours_blob_q(X,Y,phiquater,nt)
    # contours_blob_q2(X,Y,phiquater2,nt)
    # contours_ter(x,X,Y,phiCOS,phiOld,phiCOSMid,phiExact,phiMidExact,dx,dt,dy)
    # contours_solid(X,Y,phiCOS,phiOld,phiCOSMid,nt)
    # hist(X,Y,phiCOS)
    # print time.time()-start_time
    # return l1,l2,linf,dx,dy
# advection(initialProfile = stirring, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.07, nx = 100, ny = 100, nt =120)
# advection(initialProfile = stirring, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.004375, nx = 100, ny = 100, nt =300)
# advection(initialProfile = stirring, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.004375, nx = 100, ny = 100, nt =1920)
# advection(initialProfile = solid, xmin = 0., xmax = 100., ymin = 0., ymax = 100., epsilon = 0.01, dt= np.pi/720, nx = 180, ny = 180, nt =90)
# advection(initialProfile = solid, xmin = 0., xmax = 100., ymin = 0., ymax = 100., epsilon = 0.01, dt= np.pi/480, nx = 120, ny = 120, nt =60)
# advection(initialProfile = solid, xmin = 0., xmax = 100., ymin = 0., ymax = 100., epsilon = 0.01, dt= np.pi/240, nx = 60, ny = 60, nt =30)
# advection(initialProfile = Smolarkiewicz, xmin = 0., xmax = 100., ymin = 0., ymax =    100., epsilon = 0.01, dt= 2.63, nx = 100, ny = 100, nt =15)

#same resolution
advection(initialProfile = constant, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.15, nx = 50, ny = 50, nt =20)
advection(initialProfile = constant, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.075, nx = 50, ny = 50, nt =40)
advection(initialProfile = constant, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.0375, nx = 50, ny = 50, nt =80)
advection(initialProfile = constant, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.01875, nx = 50, ny = 50, nt =160)
#same Courant number 
# advection(initialProfile = constant, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.15, nx = 50, ny = 50, nt =20)
# advection(initialProfile = constant, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.075, nx = 100, ny = 100, nt =40)
# advection(initialProfile = constant, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.0375, nx = 200, ny = 200, nt =80)
# advection(initialProfile = constant, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.01875, nx = 400, ny = 400, nt =160)

#the order of accuracy
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 50, nx =150, ny = 25, nt =200)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 25, nx =300, ny = 50, nt =400)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 12.5, nx =600, ny = 100, nt =800)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 6.25, nx =1200, ny = 200, nt =1600)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 5, nx =1500, ny = 250, nt =2000)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 3.125, nx =2400, ny = 400, nt =3200)

#the same resolution 
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 50, nx =300, ny = 50, nt =200)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 100, nx =300, ny = 50, nt =100)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 200, nx =300, ny = 50, nt =50)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 500, nx =300, ny = 50, nt =20)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 1000, nx =300, ny = 50, nt =10)


# orders_solid()
# orders_blob()
# orders_SMO()
# orders_ter()