import numpy as np
import matplotlib.pyplot as plt
#---------------------------------------------------------------------------------
# Author: Yumeng Chen
# The main program to call the COSMIC splitting with different test cases
#---------------------------------------------------------------------------------

exec(open("initialConditions.py").read())
exec(open("scheme1.py").read())
def advection(initialProfile, mesh, xmin, xmax,  ymin, ymax, nx, ny, dt, nt):
    np.seterr(all='raise')

    #-----------------------
    # Basic grid information
    #-----------------------
    Lx, Ly= xmax-xmin,ymax-ymin
    dx,dy = Lx/nx, Ly/ny
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)

    #------------------------------
    # initial conditions retrieval
    #------------------------------
    change = False
    t = -1
    phiOld, cx, cy, u, v, X, Y, J  = initialProfile(x, y, xmin, ymin, nx, ny, Lx, Ly ,t, nt, dt, mesh, change)

    # print np.max(J)
    #------------------------------ 
    # diagnosis information for 
    # velocity field
    #------------------------------
    # print np.max(u), np.max(v)
    # dudx= (u[:-1,1:]-u[:-1,:-1])/dx[:-1,:-1]
    # dvdy = (v[1:,:-1]-v[:-1,:-1])/dy[:-1,:-1]
    # print 'divergence :', np.max(dudx+dvdy)
    # print 'the maximum deformational cx is ', np.max(dudx),' the maximum deformational cy is ',np.max(dvdy)
    # # print 'the maximum cx is ',np.max(cx),' the maximum cx is ', np.max(cy)
    # print 'the maximum deformational cx is ', np.max(dudx)*dt,' the maximum deformational cy is ',np.max(dvdy)*dt
    
    error = COSMIC(phiOld, cx, cy, dx, dy, xmin, ymin, u, v, X, Y, dt, nt, J, initialProfile, mesh, change)
    print error


advection(initialProfile = solid, mesh ='quad', xmin = 0., xmax = 10000.,  ymin = 0., ymax = 10000., nx = 50, ny = 50, dt = 0.02, nt = 1)