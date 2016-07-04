# import numpy as np
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
    np.seterr(all='raise')
    Lx, Ly= xmax-xmin,ymax-ymin
    dx,dy = Lx/nx, Ly/ny
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
    t = -1
    phiOld,X,Y,J,psi,UchangewithT,Y1,r,r1,r2 = initialProfile(x,y,xmin,ymin,nx,ny,Lx,Ly,t,nt,dt)
    u,v,cx,cy = streamfunction(psi,dx,dy,dt,initialProfile)
    # plt.figure(1)
    # plt.clf()
    # plt.contour(X,Y,X,np.arange(0,10000,200),colors = "k")
    # plt.contour(X,Y,Y1,np.arange(0,10000,200),colors = 'k')
    # # plt.show()
    # plt.savefig('streamfunction.pdf',transparent = True)
    #calculate the divergence32
    # print np.max(phiOld),np.min(phiOld)
    dudx= (u[:-1,1:]-u[:-1,:-1])/dx
    dvdy = (v[1:,:-1]-v[:-1,:-1])/dy
    print 'divergence :', np.max(dudx+dvdy)
    print 'the maximum Courant number is ',np.max(cx),np.max(cy),np.max(np.max(dudx*J[:-1,:-1],axis = 1)*dt,axis =0),np.max(dvdy*J[:-1,:-1])*dt
    # phi = COSMIC(phiOld,nt,epsilon,dx,dy,J,initialProfile,xmin,ymin,dt,UchangewithT,J*cx,J*cy)
    # # np.savez(str(nx)+'nt'+str(nt), phi[0], phi[1],phi[2],phi[3],phi[4],phi[5])
    phi = [phiOld,phiOld,phiOld,phiOld,phiOld,phiOld]
    f = np.load(str(nx)+'nt'+str(nt)+'.npz')
    for i in xrange(6):
        phi[i] = f['arr_'+str(i)]
    output = [int(nt/6.),int(nt/3.),int(nt/2.),int(2*nt/3.),int(5*nt/6)]
    phiExact = np.zeros_like(phi)
    # for i in xrange(5):
    # i = 4
    # print 'The maximum of phi at nt = '+str(output[i])+'is', np.max(np.round(phi[i],5)),np.min(np.round(phi[i],5))
    # phiExact[i] = solid_ana(X,Y,dx,dy,output[i],dt,nx,ny)
    # l1,l2,linf = norms(phiExact[i], phi[i], dx, dy)
    # print 'at nt = '+str(output[i])+'time the l2 norm is ', np.round(l2,5),'the linf norm is ', np.round(linf,5)
    phiExact[-1] = solid_ana(X,Y,dx,dy,nt,dt,nx,ny)
    print  'The maximum of phi at final time is ', np.max(np.round(phi[5],5)),np.min(np.round(phi[5],5))
    l1,l2,linf = norms(phiOld, phi[-1], dx, dy)
    print 'at final time the l2 norm is ', np.round(l2,5),'the linf norm is ', np.round(linf,5)
    contours_solid(X,Y,phi,phiExact,phiOld,nt)

advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 2., nx = 50, ny = 50, nt =300)
advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 1., nx = 100, ny = 100, nt =600)
advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 0.5, nx = 200, ny = 200, nt =1200)
advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 0.25, nx = 400, ny = 400, nt =2400)
# 
# 
# advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 20., nx = 50, ny = 50, nt =30)
# advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 10., nx = 100, ny = 100, nt =60)
# advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 5., nx = 200, ny = 200, nt =120)
# advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 2.5, nx = 400, ny = 400, nt =240)

# # #solid body rotation with dt
# advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 10., nx = 100, ny = 100, nt =60)
# advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 5., nx = 100, ny = 100, nt =120)
# advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 2., nx = 100, ny = 100, nt =300)
# advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 1., nx = 100, ny = 100, nt =600)
# advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 0.5, nx = 100, ny = 100, nt =1200)
# advection(initialProfile = solid, xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., epsilon = 0.01, dt= 0.25, nx = 100, ny = 100, nt =2400)




