import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#---------------------------------------------------------------------------------
# Author: Yumeng Chen
# The main program to call the COSMIC splitting with different test cases
#---------------------------------------------------------------------------------

exec(open("initialCondition.py").read())
exec(open("scheme.py").read())
exec(open("diagnostics.py").read())
exec(open("meshes.py").read())
def advection(initialProfile, mesh, xmin, xmax,  ymin, ymax, nx, ny, dt, nt):
    np.seterr(all='raise')

    #-----------------------
    # Basic grid information
    #-----------------------
    Lx, Ly= xmax-xmin,ymax-ymin
    dx,dy = Lx/nx, Ly/ny
    x_edge,y_edge = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
    x_cntr, y_cntr = x_edge[:-1] + 0.5*dx, y_edge[:-1] + 0.5*dy

    #------------------------------
    # initial conditions retrieval
    #------------------------------
    change = False
    t = -1
    phiOld, phiExact, cx, cy, u, v, X_cntr, X_edge, Y_edge, Y, J, Y_C = solid(x_edge, y_edge, x_cntr, y_cntr, t, nt, dt, mesh, change)
    print np.max(phiOld)

    # #------------------------------
    # # diagnosis information for 
    # # velocity field
    # #------------------------------
    # print np.max(u), np.max(v)
    # dudx= (u[:-1,1:]-u[:-1,:-1])/dx[:-1,:-1]
    # dvdy = (v[1:,:-1]-v[:-1,:-1])/dy[:-1,:-1]
    # print 'divergence :', np.max(dudx+dvdy)
    # print 'the maximum deformational cx is ', np.max(dudx),' the maximum deformational cy is ',np.max(dvdy)
    # # print 'the maximum cx is ',np.max(cx),' the maximum cx is ', np.max(cy)
    # print 'the maximum deformational cx is ', np.max(dudx)*dt,' the maximum deformational cy is ',np.max(dvdy)*dt
    
    #------------------------------
    # COSMIC splitting updates
    #------------------------------
    # if initialProfile == deform:
    #     change = True
    # else:
    #     change = False
    phi = COSMIC(phiOld, cx, cy, u, v, X_cntr, X_edge, Y_edge, Y, dt, nt, J, initialProfile, mesh, change)
    
    #--------------------------------------
    # save data for ploting and error norm
    #--------------------------------------
    if initialProfile == solid:
        it = [int(nt/6.),int(nt/3.),int(nt/2.),int(2*nt/3.),int(5*nt/6), nt]
        for i in xrange(5):
            print 'The max and min of phi at nt = '+str(it[i])+' is ', np.max(np.round(phi[i],5)),np.min(np.round(phi[i],5))
            l1,l2,linf = norms(phiExact[i], phi[i], dx, dy)
            print 'at nt = '+str(it[i])+', the l2 norm is ', np.round(l2,5),' the linf norm is ', np.round(linf,5)
        print  'The max and min of phi at final time is ', np.max(np.round(phi[-1],5)),np.min(np.round(phi[-1],5))
        l1,l2,linf = norms(phiOld, phi[-1], dx, dy)
        print 'at final time the l2 norm is ', np.round(l2,5),' the linf norm is ', np.round(linf,5)

        np.savez('solid_'+mesh+'_nx'+str(nx)+'nt'+str(nt), phiOld, phi[0], phi[1],phi[2],phi[3],phi[4],phi[5])
        np.savez('solid_'+mesh+'_nx'+str(nx)+'nt'+str(nt)+'_Exact',  phiExact[0], phiExact[1],phiExact[2],phiExact[3],phiExact[4],phiExact[5])
    hfont = {'fontname':'FreeSerif'}
    fig = plt.figure(1)
    plt.clf()
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    PPMsum = phi[0]+phi[1]+phi[2]+phi[3]+phi[4]+phi[5]
    Exactsum = phiExact[0]+phiExact[1]+phiExact[2]+phiExact[3]+phiExact[4]+phiExact[5]
    plt.contour(X_cntr,Y_C,np.round(phi[0],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X_cntr,Y_C,np.round(phi[1],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X_cntr,Y_C,np.round(phi[2],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X_cntr,Y_C,np.round(phi[3],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X_cntr,Y_C,np.round(phi[4],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X_cntr,Y_C,np.round(phi[5],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X_cntr,Y_C,np.round(phiExact[0],3),np.arange(0.1,1.0,0.1),colors = 'k',linestyles = 'dashed')
    plt.contour(X_cntr,Y_C,np.round(phiExact[1],3),np.arange(0.1,1.0,0.1),colors = 'k',linestyles = 'dashed')
    plt.contour(X_cntr,Y_C,np.round(phiExact[2],3),np.arange(0.1,1.0,0.1),colors = 'k',linestyles = 'dashed')
    plt.contour(X_cntr,Y_C,np.round(phiExact[3],3),np.arange(0.1,1.0,0.1),colors = 'k',linestyles = 'dashed')
    plt.contour(X_cntr,Y_C,np.round(phiExact[4],3),np.arange(0.1,1.0,0.1),colors = 'k',linestyles = 'dashed')
    plt.contour(X_cntr,Y_C,np.round(phiExact[5],3),np.arange(0.1,1.0,0.1),colors = 'k',linestyles = 'dashed')
    plt.xticks(np.arange(0,11000,2000),**hfont)
    plt.yticks(np.arange(0,11000,2000),**hfont)
    plt.gca().set_aspect('equal', adjustable='box')
    # levs = np.arange(-0.095,0.105,0.01)
    levs = np.arange(-0.95,1.,0.1)
    # color = ['#400056','#2B006B','#150081','#0000FF','#0061FF','#0092FF','#00C2FF','#00F3FF','#6DFFFF','#FFFFFF','#FFFF6D','#FFF300',\
    #     '#FFC200','#FF9200','#FF6100','#FF3100','#FF0000','#CA0600','#950B00','#601100']
    color = ['#400056','#2B006B','#150081','#0000FF','#0061FF','#0092FF','#00C2FF','#00F3FF','#6DFFFF','#FFFFFF','#FFFF6D','#FFF300',\
        '#FFC200','#FF9200','#FF6100','#FF0000','#CA0600','#950B00','#601100']
    cs = plt.contourf(X_cntr,Y_C,PPMsum-Exactsum,levels = levs,colors = color,extend="both")
    cs.cmap.set_under('#400056')
    cs.cmap.set_over('#601100')
    cs.set_clim(-1.0, 1.0)
    plt.figtext(.5, 0.12, "min = "+str(np.round(np.min(phi[-1]-phiOld),3))+', max = '+str(np.round(np.max(phi[-1]-phiOld),3)),horizontalalignment='center',
     verticalalignment='center',fontsize = 36, **hfont)
    # cb = plt.colorbar(cs,orientation = 'horizontal',fraction=0.046, pad=0.04)
    # cb.set_ticks(np.round(np.arange(-0.08,0.1,0.02),2))
    # cb.set_ticklabels(['%1.2f' % i for i in np.arange(-0.08,0.1,0.02)])
    # cb.set_ticks(np.round(np.arange(-0.8,0.9,0.1),1))
    # cb.set_ticklabels(['%1.1f' % i for i in np.arange(-0.8,0.9,0.1)])
    plt.tight_layout()
    matplotlib.rcParams.update({'font.size': 36})
    plt.show()
    # plt.savefig('SOLUTON_solid'+str(nt)+'.pdf', transparent=True)
    # if initialProfile == orography:
    #     print 'The max and min of phi at midtime is ', np.max(np.round(phi[0],5)),np.min(np.round(phi[0],5))
    #     l1,l2,linf = norms(phiExact[0], phi[0], dx, dy)
    #     print 'at midtime, the l2 norm is ', np.round(l2,5),' the linf norm is ', np.round(linf,5)
    #     print  'The max and min of phi at final time is ', np.max(np.round(phi[1],5)),np.min(np.round(phi[1],5))
    #     l1,l2,linf = norms(phiExact[-1], phi[1], dx, dy)
    #     print 'at final time the l2 norm is ', np.round(l2,5),' the linf norm is ', np.round(linf,5)
    #     np.savez('orography_'+mesh+'_nx'+str(nx)+'nt'+str(nt), phiOld, phi[0], phi[1])
    #     np.savez('orography__'+mesh+'_nx'+str(nx)+'nt'+str(nt)+'_Exact',  phiExact[0], phiExact[1])
    # if initialProfile == deform:
    #     print  'The max and min of phi at final time is ', np.max(np.round(phi[1],5)),np.min(np.round(phi[1],5))
    #     l1,l2,linf = norms(phiExact, phi[-1], dx, dy)
    #     print 'at final time the l2 norm is ', np.round(l2,5),' the linf norm is ', np.round(linf,5)
    #     np.savez('deform_'+mesh+'_nx'+str(nx)+'nt'+str(nt), phiOld, phi[0], phi[1],phi[2],phi[3],phi[4])
    #     np.savez('deform_'+mesh+'_nx'+str(nx)+'nt'+str(nt)+'_Exact',  phiExact[0], phiExact[1],phiExact[2],phiExact[3],phiExact[4])


advection(initialProfile = solid, mesh ='quad', xmin = 0., xmax = 10000.,  ymin = 0., ymax = 10000., nx = 100, ny = 100, dt = 10., nt = 60)
# advection(initialProfile = deform, mesh = 'orthogo', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 10, ny = 5, dt = 0.2, nt = 25)