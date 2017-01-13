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
# exec(open("plot.py").read())
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
    # initial condition retrieval
    #------------------------------
    change = False
    t = -1
    phiOld, phiExact, cx, cy, J, J_s, J_w = initialProfile(x_edge, y_edge, x_cntr, y_cntr, t, nt, dt, mesh, change)
    
    # ------------------------------
    # COSMIC splitting updates
    # ------------------------------
    if initialProfile == deform:
        change = True
    else:
        change = False
    phi = COSMIC(phiOld, cx, cy, dx, dy, xmin, ymin, dt, nt, J, J_s, J_w, initialProfile, mesh, change)

    # #--------------------------------------
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
        
    if initialProfile == orography:
        print 'The max and min of phi at midtime is ', np.max(np.round(phi[0],5)),np.min(np.round(phi[0],5))
        l1,l2,linf = norms(phiExact[0], phi[0], dx, dy)
        print 'at midtime, the l2 norm is ', np.round(l2,5),' the linf norm is ', np.round(linf,5)
        print  'The max and min of phi at final time is ', np.max(np.round(phi[1],5)),np.min(np.round(phi[1],5))
        l1,l2,linf = norms(phiExact[-1], phi[1], dx, dy)
        print 'at final time the l2 norm is ', np.round(l2,5),' the linf norm is ', np.round(linf,5)
        np.savez('orography_'+mesh+'_nx'+str(nx)+'nt'+str(nt), phiOld, phi[0], phi[1])
        np.savez('orography_'+mesh+'_nx'+str(nx)+'nt'+str(nt)+'_Exact',  phiExact[0], phiExact[1])

    if initialProfile == deform:
        print  'The max and min of phi at final time is ', np.max(np.round(phi[1],5)),np.min(np.round(phi[1],5))
        l1,l2,linf = norms(phiExact, phi[-1], dx, dy)
        print 'at final time the l2 norm is ', np.round(l2,5),' the linf norm is ', np.round(linf,5)
        np.savez('deform_'+mesh+'_nx'+str(nx)+'nt'+str(nt), phiOld, phi[0], phi[1],phi[2],phi[3],phi[4])
        np.savez('deform_'+mesh+'_nx'+str(nx)+'nt'+str(nt)+'_Exact',  phiExact[0], phiExact[1],phiExact[2],phiExact[3],phiExact[4])
