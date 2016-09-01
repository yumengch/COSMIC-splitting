import numpy as np
def COSMIC(phiOld, c, u, X_cntr, X_edge, dt, nt):
#---------------------------------------------------------------------------------
# Author: Yumeng Chen
# Scheme: COSMIC splitting
# Basis: Leonard, B. P., A. P. Lock, and M. K. MacVean, 1996: Conservative explicit 
#           unrestricted-time-step multidimensional constancy-preserving advection schemes. 
# Input: phiOld
# Output: phi at time nt*dt
#---------------------------------------------------------------------------------

    #-----------------------
    # Basic grid information
    #-----------------------
    nx = len(phiOld)
    xmin = np.min(X_edge)
    dx = X_edge[1] - X_edge[0]
    xmax = xmin + nx*dx
    Lx = xmax - xmin

    #-----------------------------------------------
    # midpoint values and mass flux at cell boundary
    #-----------------------------------------------
    phi_mid= np.zeros_like(phiOld)
    mass = np.zeros([nx])

    #---------------------------------
    # XA, YA advective-form operator 
    # XC, YC conservative-form operator
    # X,Y denote direction,
    # A,C means advective/conservative 
    #---------------------------------
    XA = np.zeros_like(phiOld)

    XC = np.zeros_like(phiOld)

    idx_x = np.zeros_like(phiOld).astype(int)

    r_x = np.zeros_like(phiOld)

    #-------------------------------------------------------
    # phi_AX, phi_AY:inner operator (advective-form) update
    # XC_AY, YC_AX: cross term operator updates
    #-------------------------------------------------------
    phi_AX = np.zeros_like(phiOld)

    XC_AY = np.zeros_like(phiOld)


    #--------------------------------
    #the updated cell-average values
    #--------------------------------
    phi = np.zeros_like(phiOld)


    #---------------
    # time updates
    #---------------
    for t in xrange(int(nt)):
        #---------------------------
        # find the departure points
        #---------------------------
        idx_x, r_x = departure_x(X_edge, xmin, xmax, u, dx, dt, Lx)

        #----------------------------------------------------------
        # advective operator and non-cross term conservative operator 
        # updates in x direction
        #----------------------------------------------------------
        #------------------
        # 1D PPM updates
        #------------------
        phi_mid, mass= PPM(phiOld, c, nx, dx, idx_x, r_x)
        #--------------------------------------
        # mass flux at each cell boundary
        #--------------------------------------
        OUT = flux(nx, dx, c, phi_mid, mass, idx_x, r_x)
        #----------------------------------------------
        #  adavective and conservative operator updates
        #----------------------------------------------    
        XC = conservative(nx, c, OUT)
        # XA = advective(nx, c, OUT)



        #-------------------------------------------------------
        # Final COSMIC splitting updates
        #------------------------------------------------------- 
        phi = phiOld+XC

        phiOld = phi.copy() #update the time step

        #----------------------------------------
        # print the time steps and maximum value
        #----------------------------------------
        print 'at ', t,' time step, the maximum of phi is ', np.max(phiOld)

    return phi


def departure_x(X, xmin, xmax, u, dx, dt, L):
    r = np.zeros_like(u)
    x_depart = X - u*dt

    for i in xrange(len(x_depart)):
        while x_depart[i] <= xmin:
            x_depart[i] += L
        while x_depart[i] >= xmax:
            x_depart[i] -= L

    idx = np.where(abs(u[i])<=10e-6, np.ceil((x_depart - xmin)/dx).astype(int) ,np.floor((x_depart - xmin)/dx).astype(int))
    # idx = np.floor((x_depart - xmin)/dx).astype(int)
    X = np.append(X, xmax)
    for i in xrange(len(X)-1):
        if u[i]>0:
            r[i] = (X[idx[i]+1] - x_depart[i])/dx
        else:
            r[i] = -(x_depart[i] - X[idx[i]])/dx
    
    return idx, r
    

def PPM(phiOld, c, nx, dx, idx, c_r):
#---------------------------------------------------------------------------------
# Author: Yumeng Chen
# Scheme:  Piecewise Parabolic Method in 1D
# Basis: Colella, P., and P. R. Woodward, 1984: The Piecewise Parabolic Method PPM 
#           for gas-dynamical simulations. J. Comput. Phys.,.
# Input: phiOld
# Output: phi at time nt+1
#---------------------------------------------------------------------------------

    #-----------------------
    # Basic grid information
    #-----------------------
    dmphi = np.zeros_like(phiOld)           # phi increment as in PPM paper
    phi_r = np.zeros_like(phiOld)           # phi at j+1/2 boundary
    phi_l = np.zeros_like(phiOld)           # phi at j-1/2 
    phi_6 = np.zeros_like(phiOld)           # the difference between phi,j and average of phi_l and phi_r
    daj = np.zeros_like(phiOld)             # the difference between the right and left boundary (for limiters)
    phi_mid = np.zeros_like(phiOld)         # final phi at j+1/2
    mass = np.zeros(nx)
    c = np.append(c, c[0])
    idx = np.append(idx,idx[0])
    c_r = np.append(c_r, c_r[0])
    #---------------------------------------
    # phi increment as in PPM paper
    #---------------------------------------
    dmphi[1:-1] = 0.5*(phiOld[2:]-phiOld[:-2])
    #-------------------------------------------------------
    # periodic boundary value updates
    #-------------------------------------------------------
    dmphi[0] = 0.5*(phiOld[1]-phiOld[-1])
    dmphi[-1] = 0.5*(phiOld[0]-phiOld[-2]) 



    #------------------------------------------------------------
    # phi at j-1/2 and j+1/2
    #------------------------------------------------------------
    phi_l[1:] = (0.5*(phiOld[1:]+phiOld[:-1]) + (dmphi[:-1]-dmphi[1:])/6.)
    phi_l[0] = (0.5*(phiOld[0]+phiOld[-1]) + (dmphi[-1]-dmphi[0])/6.)

    phi_r[:-1] = (0.5*(phiOld[1:]+phiOld[:-1]) + (dmphi[:-1]-dmphi[1:])/6.)
    phi_r[-1] = (0.5*(phiOld[0]+phiOld[-1]) + (dmphi[-1]-dmphi[0])/6.)

    #------------------------------------------------------------
    # piecewise parabolic subcell reconstruction
    #-----------------------------------------------------------
    daj = phi_r - phi_l
    phi_6 = 6*(phiOld-0.5*(phi_l+phi_r))

    #------------------------------------------------------------
    # PPM update to get phi at j+1/2
    #-----------------------------------------------------------
    for i in xrange(nx):
        # print idx[i+1]
        if c[i+1]>= 0:
            phi_mid[i] = phi_r[idx[i+1]]-0.5*c_r[i+1]*(daj[idx[i+1]]-(1-2*c_r[i+1]/3.)*phi_6[idx[i+1]])
        else:
            phi_mid[i] = phi_l[idx[i+1]]-0.5*c_r[i+1]*(daj[idx[i+1]]+(1+2*c_r[i+1]/3.)*phi_6[idx[i+1]])
  
    #-------------------
    # cumulative mass 
    #-------------------    
    mass[0] = phiOld[0]*dx
    for j in xrange(1,nx):
        mass[j] = mass[j-1] +phiOld[j]*dx

    return phi_mid, mass

def flux(nx,dx,c,phi_mid,mass, idx, c_r):
#-----------------------------------------------------
# function: Flux calculation at j+1/2
#-----------------------------------------------------

    #---------------------------------------
    # Flux at j+1/2
    #---------------------------------------
    OUT = np.zeros_like(phi_mid)
    c_r = np.append(c_r, c_r[0])
    c = np.append(c, c[0])
    idx = np.append(idx,idx[0])
    for i in xrange(nx):
        if c[i+1]>0:          # velocity u > 0
            if i>idx[i+1]:           # if the departure cell is at the west of predicted cell 
                OUT[i] = (phi_mid[i]*c_r[i+1]*dx+(mass[i]-mass[idx[i+1]]))/dx
            elif i==idx[i+1]:        # if the departure cell is at the position of predicted cell
                OUT[i] = phi_mid[i]*c_r[i+1]*dx/dx
            else:             # if the departure cell is at the east of predicted cell
                OUT[i] = (phi_mid[i]*c_r[i+1]*dx+mass[i]+mass[-1]-mass[idx[i+1]])/dx
        elif c[i+1]<0:  
            k = np.floor(idx[i+1]-1)%nx      # velocity u < 0
            if i> k:           # if the departure cell is at the east of predicted cell  
                OUT[i] = (-phi_mid[i]*c_r[i+1]*dx+mass[-1]-mass[i]+mass[k])/dx
            elif i ==k:     # if the departure cell is at the position of predicted cell
                OUT[i] = -phi_mid[i]*c_r[i+1]*dx/dx
            elif i<k:       # if the departure cell is at the west of predicted cell 
                OUT[i] = (-phi_mid[i]*c_r[i+1]*dx-mass[i]+mass[k])/dx
        else:               # if velcity = 0
            OUT[i] = 0
    return OUT

def advective(nx,c,OUT):
#-----------------------------------------------------
# function: advecitve operator calculation
#-----------------------------------------------------

    #---------------------------------------
    # advective operator 
    #---------------------------------------
    AX = np.zeros_like(OUT)
    c = np.append(c, c[0])
    for i in xrange(nx):
        #----------------------------------------------------
        # the if statement is to update by upwind velocity
        #----------------------------------------------------
        if c[i] >= 0 and c[i+1] >0:
            AX[i] = OUT[i-1]-c[i]*OUT[i]/c[i+1]
        elif c[i] < 0 and c[i+1] <= 0 :
            AX[i] = OUT[i]-c[i+1]*(OUT[i-1]/c[i])
        else:
            AX[i] = 0
    return AX[:]

def conservative(nx,c,OUT):
#-----------------------------------------------------
# function: conservative operator calculation
#-----------------------------------------------------
    XC = np.zeros_like(OUT)
    c = np.append(c, c[0])
    #-----------------------------------------------------
    # update by influx - outflux in each cell
    #-----------------------------------------------------
    for i in xrange(nx):
        if c[i+1]<=0 and c[i]<=0:
            XC[i] = (OUT[i]-OUT[i-1])
        elif c[i+1]>=0 and c[i]>=0:
            XC[i] = (OUT[i-1]-OUT[i])
        elif c[i+1]>=0 and c[i]<=0:
            XC[i] = -(OUT[i]+OUT[i-1])
        elif c[i+1]<=0 and c[i]>=0:
            XC[i] = (OUT[i]+OUT[i-1])
    return XC[:]


