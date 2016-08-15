import numpy as np
def COSMIC(u, dt, x1,L):


    #-----------------------
    # Basic grid information
    #-----------------------
    # nx = len(phiOld)
    # xmax = xmin + (nx-1)*dx
    # print xmax
    # Lx = xmax-xmin
    # x = np.linspace(xmin,xmax,nx)

    #--------------------------------
    #the updated cell-average values
    #--------------------------------
    # phi = np.zeros_like(phiOld)

    # x_depart_phys = np.zeros_like(phiOld)

    # x_depart_compt = np.zeros_like(phiOld)
    # for t in xrange(int(nt)):
        #----------------------
        # departure point calculation
        #--------------------------

        # for i in xrange(nx):
        #     y_depart_phys[:,i],y_depart_compt[:,i] = departure(y[:-1,i], Y[:-1,i], v[:,i], dt, cy[:,i], J[:,i], Ly, dy)
        # for j in xrange(ny):
    x_depart_phys = departure(u, dt, x1, L)
    
    #----------------------
    # transform computational departure points to physical domain
    #--------------------------
    # Y_depart = comput_SB(X,y_depart_compt,ymax,f_quad(X,Lx, Ly), Ly)

    print x_depart_phys
        # error = x_depart_phys - X_depart
        # for i in xrange(nx):
        #     print i, np.max(error[:,i])
    return x_depart_phys

    

def departure(u, dt, X, L):
    #---------------------------------------
    # Integer and remnant Courant number
    #---------------------------------------
    print X, u*dt
    x_depart_phys = X - u*dt
    for i in xrange(len(x_depart_phys)):
        while x_depart_phys[i] < X[0]:
            x_depart_phys[i] = x_depart_phys[i] + L
        while x_depart_phys[i] > X[-1]:
            print x_depart_phys[i]
            x_depart_phys[i] = x_depart_phys[i] - L


    # x_depart_compt = x - J*c*dx
    
    # for i in xrange(len(x_depart_compt)):
    #     while x_depart_compt[i] < 0.:
    #         x_depart_compt[i] = 0.
    #     while x_depart_compt[i] > L:
    #         x_depart_compt[i] = L
 
    return x_depart_phys#, x_depart_compt


def f(x):
    return -0.4
def comput(x,L):
    #-------------------------------------
    # unequidistant computational domain 
    #-------------------------------------
    fx = f(x)
    # return np.where(x>0.5*L, fx+(x-0.5*L)*(L-fx)/(0.5*L),x*fx/(0.5*L))
    return np.where(x>=0, fx+x*(1-fx/xmax),fx+x*(1-fx/xmin))