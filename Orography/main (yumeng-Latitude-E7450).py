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
    # print Y[1:,60:90] - Y[:-1,60:90]
    #calculate the divergence
    print np.max(phiOld),np.min(phiOld)
    dudx= (u[:-1,1:]-u[:-1,:-1])/dx
    dvdy = (v[1:,:-1]-v[:-1,:-1])/dy
    print 'the maximum and minimum Courant number is ',np.max(cx),np.max(cy),np.max(dudx)*dt,np.max(dvdy)*dt
    # phiCOS,phiCOSMid = COSMIC(phiOld,nt,epsilon,dx,dy,J,initialProfile,xmin,ymin,dt,UchangewithT,J*cx,J*cy)
    # np.savez(str(nx)+'nt'+str(nt), phiCOS = phiCOS, phiCOSMid = phiCOSMid)
    f = np.load(str(nx)+'nt'+str(nt)+'.npz')
    phiCOS = f['phiCOS']
    phiCOSMid = f['phiCOSMid']
    print 'max: ', np.max(np.round(phiCOS,5)),'min: ',np.min(np.round(phiCOS,5))
    # print np.max(np.round(phiCOSMid,5)),np.min(np.round(phiCOSMid,5))
    phiExact = analytical(X,Y,dx,dy,nt,dt,nx,ny)
    phiMidExact = analytical(X,Y,dx,dy,int(nt/2.),dt,nx,ny)
    # l1,l2,linf = norms(phiExact, phiCOS, dx, dy)
    # l1Mid,l2Mid,linfMid = norms(phiMidExact, phiCOSMid, dx, dy)
    # print 'the l2 norm is ', np.round(l2,5),'the linf norm is ', np.round(linf,5)
    # print 'the l2 norm is ', np.round(l2Mid,5),'the linf norm is ', np.round(linfMid,5)
    contours_ter(x,X,Y,phiCOS,phiOld,phiCOSMid,phiExact,phiMidExact,dx,dt,dy)   


#the order of accuracy
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 50, nx =150, ny = 25, nt =200)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 25, nx =300, ny = 50, nt =400)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 12.5, nx =600, ny = 100, nt =800)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 6.25, nx =1200, ny = 200, nt =1600)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 5, nx =1500, ny = 250, nt =2000)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 3.125, nx =2400, ny = 400, nt =3200)

# #the same resolution 
advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 50, nx =300, ny = 50, nt =200)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 100, nx =300, ny = 50, nt =100)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 200, nx =300, ny = 50, nt =50)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 500, nx =300, ny = 50, nt =20)
# advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 1000, nx =300, ny = 50, nt =10)