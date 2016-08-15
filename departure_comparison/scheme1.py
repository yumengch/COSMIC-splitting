import numpy as np
def COSMIC(cx, cy, dx, dy, xmin, ymin, u, v, X, Y, dt, nt, J):


    #-----------------------
    # Basic grid information
    #-----------------------
    nx, ny = len(u[0,:-1]), len(u[:,0])
    xmax,ymax = xmin + nx*dx,ymin+ny*dy

    Lx, Ly= xmax-xmin,ymax-ymin
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
    x,y = np.meshgrid(x,y)

    #--------------------------------
    #the updated cell-average values
    #--------------------------------
    y_depart_phys = np.zeros_like(Y)
    x_depart_phys = np.zeros_like(X)
    y_depart_compt = np.zeros_like(Y)
    x_depart_compt = np.zeros_like(X)
    for t in xrange(int(nt)):
        #----------------------
        # departure point calculation
        #--------------------------

        for i in xrange(nx):
            y_depart_phys[:,i],y_depart_compt[:,i] = departure(y[:,i], Y[:,i], v[:,i], dt, cy[:,i], J[:,i], Ly, dy)
        for j in xrange(ny):
            x_depart_phys[j,:],x_depart_compt[j,:] = departure(x[j,:], X[j,:], u[j,:], dt, cx[j,:], J[j,:], Lx, dx)
        
        #----------------------
        # transform computational departure points to physical domain
        #--------------------------
        Y_depart = comput_SB(X,y_depart_compt,ymax,f_quad(X,Lx, Ly), Ly)

        
        error = y_depart_phys - Y_depart
        for j in xrange(ny):
            print j, error[50,j]
    return np.max(error)

    

def departure(x, X, u, dt, c, J, L, dx):
    #---------------------------------------
    # Integer and remnant Courant number
    #---------------------------------------
    x_depart_phys = X - u*dt
    for i in xrange(len(x_depart_phys)):
        while x_depart_phys[i] < 0.:
            x_depart_phys[i] += L
        while x_depart_phys[i] > L:
            x_depart_phys[i] -= L

    x_depart_compt = x - J*c*dx
    
    for i in xrange(len(x_depart_compt)):
        while x_depart_compt[i] < 0.:
            x_depart_compt[i] += L
        while x_depart_compt[i] > L:
            x_depart_compt[i] -= L
 
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