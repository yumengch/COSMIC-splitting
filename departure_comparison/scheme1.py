import numpy as np
def COSMIC(phiOld, cx, cy, dx, dy, xmin, ymin, u, v, X, Y, dt, nt, J, initialProfile, mesh, change):


    #-----------------------
    # Basic grid information
    #-----------------------
    nx, ny = len(phiOld[0,:]) - 1, len(phiOld[:,0]) - 1
    xmax,ymax = xmin + nx*dx,ymin+ny*dy
    Lx, Ly= xmax-xmin,ymax-ymin 
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
    x,y = np.meshgrid(x,y)

    #--------------------------------
    #the updated cell-average values
    #--------------------------------
    phi = np.zeros_like(phiOld)
    y_depart_phys = np.zeros_like(phiOld)
    x_depart_phys = np.zeros_like(phiOld)
    y_depart_compt = np.zeros_like(phiOld)
    x_depart_compt = np.zeros_like(phiOld)

    for t in xrange(int(nt)):
        if change:
            cx, cy = initialProfile(x, y, xmin, ymin, nx, ny, Lx, Ly ,t, nt, dt, mesh, change)

        for i in xrange(nx):
            y_depart_phys[:,i],y_depart_compt[:,i] = departure(y[:,i], Y[:,i], v[:,i], dt, cy[:,i], J[:,i], Ly, dy)
        for j in xrange(ny):
            x_depart_phys[j,:],x_depart_compt[j,:] = departure(x[j,:], X[j,:], u[j,:], dt, cx[j,:], J[j,:], Lx, dx)
        
        Y_depart = comput_SB(X,y_depart_compt,ymax,f_quad(X,Lx, Ly), Ly)

        error = y_depart_phys - Y_depart
        for i in xrange(nx):
            print i, np.max(error[:,i])
    return np.max(error)

    

def departure(x, X, u, dt, c, J, L, dx):


    #---------------------------------------
    # Integer and remnant Courant number
    #---------------------------------------
    x_depart_phys = X - u*dt
    x_depart_phys = np.where(x_depart_phys < X[0], X[0], x_depart_phys)
    x_depart_phys = np.where(x_depart_phys > X[-1], X[-1], x_depart_phys)
    # while len(np.where(np.logical_or(x_depart_phys < X[0], x_depart_phys > X[-1]))[0]) > 0:
    #     x_depart_phys = np.where(x_depart_phys < X[0], x_depart_phys + L, x_depart_phys)
    #     x_depart_phys = np.where(x_depart_phys > X[-1], x_depart_phys - L, x_depart_phys)

    x_depart_compt = x - J*c*dx
    x_depart_compt = np.where(x_depart_compt < x[0], x[0], x_depart_compt)
    x_depart_compt = np.where(x_depart_compt > x[-1], x[-1], x_depart_compt)
    # while len(np.where(np.logical_or(x_depart_compt < x[0], x_depart_compt > x[-1]))[0]) > 0:
    #     x_depart_compt = np.where(x_depart_compt < x[0], x_depart_compt + L, x_depart_compt)
    #     x_depart_compt = np.where(x_depart_compt > x[-1], x_depart_compt - L, x_depart_compt)
 
    return x_depart_phys, x_depart_compt


def f_quad(x,Lx, Ly):
    #-------------------------
    # quadratic mesh function
    #-------------------------
    return 2*(x-0.5*Lx)**2/Lx/np.sqrt(3) + 0.5*Ly-0.5*0.25*Lx/np.sqrt(3)

def comput_SB(X, Y, ymax, fx, Ly):

    return np.where(Y>0.5*Ly, fx+(Y-0.5*Ly)*(Ly-fx)/(0.5*Ly),Y*fx/(0.5*Ly))

def f_V(x,Lx, Ly):
    #-------------------------
    # V shape mesh function
    #-------------------------
    return np.where(x<=0.5*Lx,
                    (-1/np.sqrt(3))*x+Lx/4./np.sqrt(3)+0.5*Lx,
                    np.where(x>0.5*Lx,
                    1/np.sqrt(3)*(x-0.5*Lx)-Lx/4./np.sqrt(3)+0.5*Lx,0))