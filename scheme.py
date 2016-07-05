import numpy as np
def COSMIC(phiOld,nt,epsilon,dx,dy,J,initialProfile,xmin,ymin,dt,UchangewithT,cx,cy):
#---------------------------------------------------------------------------------
# Author: Yumeng Chen
# Scheme: COSMIC splitting
# Basis: Leonard, B. P., A. P. Lock, and M. K. MacVean, 1996: Conservative explicit 
#           unrestricted-time-step multidimensional constancy-preserving advection schemes. 
#               Mon. Wea. Rev., 124, 2585–2606.
# Input: phiOld
# Output: phi at time nt*dt
#---------------------------------------------------------------------------------

    #-----------------------
    # Basic grid information
    #-----------------------
    nx, ny = len(phiOld[0,:]) - 1, len(phiOld[:,0]) - 1
    xmax,ymax = xmin + nx*dx,ymin+ny*dy
    Lx, Ly= xmax-xmin,ymax-ymin 
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)

    #-----------------------------------------------
    # midpoint values and mass flux at cell boundary
    #-----------------------------------------------
    phi_mid= np.zeros_like(phiOld)
    mass = np.zeros([ny,nx])

    #---------------------------------
    # XA, YA advective-form operator 
    # XC, YC conservative-form operator
    # X,Y denote direction,
    # A,C means advective/conservative 
    #---------------------------------
    XA = np.zeros_like(phiOld)
    YA = np.zeros_like(phiOld)
    XC = np.zeros_like(phiOld)
    YC = np.zeros_like(phiOld)

    #-------------------------------------------------------
    # phi_AX, phi_AY:inner operator (advective-form) update
    # XC_AY, YC_AX: cross term operator updates
    #-------------------------------------------------------
    phi_AX = np.zeros_like(phiOld)
    phi_AY = np.zeros_like(phiOld)
    XC_AY = np.zeros_like(phiOld)
    YC_AX = np.zeros_like(phiOld)

    #--------------------------------
    #the updated cell-average values
    #--------------------------------
    phi = np.zeros_like(phiOld)


    #---------------
    # time updates
    #---------------
    for t in xrange(int(nt)):
        #---------------------------------------------------
        # velocity updates for deformational flow test case
        #---------------------------------------------------
        if UchangewithT == True:
            psi = initialProfile(x,y,xmin,ymin,nx,ny,Lx,Ly,t,nt,dt)
            u,v,cx,cy = streamfunction(psi,dx,dy,dt,initialProfile)

        #-------------------------------------------------------------
        # advective operator and non-cross term conservative operator 
        # updates in y direction
        #-------------------------------------------------------------
        for i in xrange(nx):
            #------------------
            # 1D PPM updates
            #------------------
            phi_mid[:,i],mass[:,i] = PPM((1/J[:,i])*phiOld[:,i],cy[:,i],ny,epsilon,dy)
            #--------------------------------------
            # mass flux at each cell boundary
            #--------------------------------------
            OUT = flux(ny,dy,cy[:,i],phi_mid[:-1,i],mass[:,i])
            #----------------------------------------------
            #  adavective and conservative operator updates
            #----------------------------------------------
            YC[:-1,i] = conservative(ny,cy[:,i],OUT)
            YA[:-1,i] = advective(ny,cy[:,i],OUT)
        #---------------------------------------
        # periodic boundary value updates
        #---------------------------------------
        YA[-1,:],YA[:,-1] = YA[0,:],YA[:,0]
        YC[-1,:],YC[:,-1] = YC[0,:],YC[:,0]


        #----------------------------------------------------------
        # advective operator and non-cross term conservative operator 
        # updates in x direction
        #----------------------------------------------------------
        for j in xrange(ny):
            #------------------
            # 1D PPM updates
            #------------------
            phi_mid[j,:], mass[j,:]= PPM((1/J[j,:])*phiOld[j,:],cx[j,:],nx,epsilon,dx)
            #--------------------------------------
            # mass flux at each cell boundary
            #--------------------------------------
            OUT = flux(nx,dx,cx[j,:],phi_mid[j,:-1],mass[j,:])
            #----------------------------------------------
            #  adavective and conservative operator updates
            #----------------------------------------------    
            XC[j,:-1] = conservative(nx,cx[j,:],OUT)
            XA[j,:-1] = advective(nx,cx[j,:],OUT)
        #---------------------------------------------------------
        # periodic boundary value updates
        #---------------------------------------------------------            
        XA[-1,:], XA[:,-1] = XA[0,:], XA[:,0]
        XC[-1,:], XC[:,-1] = XC[0,:], XC[:,0]

        #---------------------------------------------------------
        # advective operator updates
        #---------------------------------------------------------    
        phi_AX = phiOld + J*XA
        phi_AY = phiOld + J*YA

        #---------------------------------------------------------------
        # conservative operator with cross-term updates in y direction
        #---------------------------------------------------------------
        for i in xrange(nx):
            #------------------
            # 1D PPM updates
            #------------------
            phi_mid[:,i],mass[:,i] = PPM((1/J[:,i])*phi_AX[:,i],cy[:,i],ny,epsilon,dy)
            #---------------------------------
            # mass flux at each cell boundary
            #---------------------------------
            OUT = flux(ny,dy,cy[:,i],phi_mid[:-1,i],mass[:,i])
            #-------------------------------
            # conservative operator updates
            #-------------------------------
            YC_AX[:-1,i] = conservative(ny,cy[:,i],OUT)
        #----------------------------------
        # periodic boundary value updates
        #----------------------------------
        YC_AX[-1,:], YC_AX[:,-1] = YC_AX[0,:], YC_AX[:,0]

        #---------------------------------------------------------------
        # conservative operator with cross-term updates in x direction
        #---------------------------------------------------------------
        for j in xrange(ny):
            #------------------
            # 1D PPM updates
            #------------------
            phi_mid[j,:], mass[j,:]= PPM((1/J[j,:])*phi_AY[j,:],cx[j,:],nx,epsilon,dx)
            #----------------------------------
            # mass flux at each cell boundary
            #----------------------------------
            OUT = flux(nx,dx,cx[j,:],phi_mid[j,:-1],mass[j,:])
            #----------------------------------------------
            # conservative operator updates
            #----------------------------------------------                
            XC_AY[j,:-1] = conservative(nx,cx[j,:],OUT)
        #-------------------------------------------------------
        # periodic boundary value updates
        #------------------------------------------------------- 
        XC_AY[-1,:], XC_AY[:,-1]= XC_AY[0,:], XC_AY[:,0]

        #-------------------------------------------------------
        # Final COSMIC splitting updates
        #------------------------------------------------------- 
        phi = phiOld+J*(0.5*(XC+XC_AY)+0.5*(YC+YC_AX))

        #-------------------------------------------------------
        # periodic boundary value updates
        #------------------------------------------------------- 
        phi[-1,:], phi[:,-1] = phi[0,:], phi[:,0]
        phiOld = phi.copy() #update the time step


        #-----------------------------
        # intermediate value storage
        #----------------------------- 
        if t == int(nt/6)-1:
            phi0 = phiOld.copy()
        if t == int(nt/3)-1:
            phi1 = phiOld.copy()
        if t == int(nt/2)-1:
            phi2 = phiOld.copy()
        if t == int(2*nt/3)-1:
            phi3 = phiOld.copy()
        if t == int(5*nt/6)-1:
            phi4 = phiOld.copy()

        #----------------------------------------
        # print the time steps and maximum value
        #----------------------------------------
        print t,np.max(phiOld)
    return [phi0,phi1,phi2,phi3,phi4,phi]

def PPM(phiOld, c, nx, dx, epsilon = 0.01, eta1=20, eta2 = 0.05):
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
    phiOld = np.append(phiOld, [phiOld[1]]) # one ghosted cell for PPM interpolation
    dmphi = np.zeros_like(phiOld)           # phi increment as in PPM paper
    phi_r = np.zeros_like(phiOld)           # phi at j+1/2 boundary
    phi_l = np.zeros_like(phiOld)           # phi at j-1/2 boundary
    d2phi = np.zeros_like(phiOld)           # the second derivative of phi (for limiters)
    eta = np.zeros_like(phiOld)             # the weight of discontinuity (for flux limiters)
    phi_6 = np.zeros_like(phiOld)           # the difference between phi,j and average of phi_l and phi_r
    daj = np.zeros_like(phiOld)             # the difference between the right and left boundary (for limiters)
    phi_mid = np.zeros_like(phiOld)         # final phi at j+1/2
    mass = np.zeros(nx)


    #---------------------------------------
    # Integer and remnant Courant number
    #---------------------------------------
    c_N = c.astype(int)
    c_r = c-c_N

    #---------------------------------------
    # phi increment as in PPM paper
    #---------------------------------------
    dmphi[1:-1] = 0.5*(phiOld[2:]-phiOld[:-2])
    #-------------------------------------------------------
    # periodic boundary value updates
    #-------------------------------------------------------
    dmphi[0] = dmphi[-2]
    dmphi[-1] = dmphi[1]
    
    #------------------------------------------------------------
    # #Calculate the second derivative of phi(centred in space)
    # Determine the eta to weigh the discontinuity
    # limiters in PPM
    #-----------------------------------------------------------
    # d2phi[1:-1] = (phiOld[2:] - 2*phiOld[1:-1] + phiOld[:-2])/(6*dx**2)
    # # # Updated boundary condition
    # d2phi[0] = d2phi[-2]
    # d2phi[-1] = d2phi[1]
    
    # for j in xrange(1,nx+1):
    #     if -d2phi[j+1]*d2phi[j-1]>0 and\
    #             np.abs(phiOld[j+1]-phiOld[j-1])-\
    #                 epsilon*min(np.abs(phiOld[j+1]),np.abs(phiOld[j-1])) > 0:
    #         eta[j] = -(d2phi[j+1]-d2phi[j-1])*(dx**2)/(phiOld[j+1]-phiOld[j-1])
    #     else:
    #         eta[j] = 0
    #     eta[j] = max(0,min(eta1*(eta[j]-eta2,1)))
    # eta[0] = eta[-2]
    # eta[-1] = eta[1]
    

    #--------------
    # no limiters
    #--------------
    eta[:] = 0

    #------------------------------------------------------------
    # phi at j-1/2 and j+1/2
    #------------------------------------------------------------
    phi_l[1:] = (0.5*(phiOld[1:]+phiOld[:-1]) + (dmphi[:-1]-dmphi[1:])/6.)*(1-eta[1:])+(phiOld[:-1] + 0.5*dmphi[:-1])*eta[1:]
    phi_l[0] = phi_l[-2] # boundary updates

    phi_r[:-1] = (0.5*(phiOld[1:]+phiOld[:-1]) + (dmphi[:-1]-dmphi[1:])/6.)*(1-eta[:-1])+(phiOld[1:] - 0.5*dmphi[1:])*eta[:-1]
    phi_r[-1] =phi_r[1] # boundary updates

    #------------------------------------------------------------
    # avoid extreme conditions that phi,j+1/2 is not phi_l,j and phi_r,j-1
    # Determine the eta to weigh the discontinuity
    # limiters in PPM
    #-----------------------------------------------------------
    # for j in xrange(nx+2):
    #     if (phi_r[j]-phiOld[j])*(phiOld[j]-phi_l[j])<=0:
    #         phi_l[j] = phi_r[j] = phiOld[j]
    #     if (phi_r[j]-phi_l[j])*(phiOld[j]-0.5*(phi_l[j]+phi_r[j]))>(phi_r[j]-phi_l[j])**2/6.:
    #         phi_l[j] = 3*phiOld[j]-2*phi_r[j]
    #     if -(phi_r[j]-phi_l[j])**2/6.>(phi_r[j]-phi_l[j])*(phiOld[j]-0.5*(phi_r[j]+phi_l[j])):
    #         phi_r[j] = 3*phiOld[j]-2*phi_l[j]


    #------------------------------------------------------------
    # piecewise parabolic subcell reconstruction
    #-----------------------------------------------------------
    daj = phi_r - phi_l
    phi_6 = 6*(phiOld-0.5*(phi_l+phi_r))

    #------------------------------------------------------------
    # PPM update to get phi at j+1/2
    #-----------------------------------------------------------
    for i in xrange(nx):
        k= np.floor(i-c_N[i+1])%nx                               # departure points if u >= 0 
        k1 = np.floor(i+1-c_N[i+1])%nx                           # departure points if u < 0 
        if dc[i+1]>= 0:
            phi_mid[i] = phi_r[k]-0.5*c_r[i+1]*(daj[k]-(1-2*c_r[i+1]/3.)*phi_6[k])
        else:
            phi_mid[i] = phi_l[k1]-0.5*c_r[i+1]*(daj[k1]+(1+2*c_r[i+1]/3.)*phi_6[k1])
    # boundary updates
    phi_mid[-2] = phi_mid[0]
    phi_mid[-1] = phi_mid[1]
  
    #-------------------
    # cumulative mass 
    #-------------------    
    mass[0] = phiOld[0]*dx
    for j in xrange(1,nx):
        mass[j] = mass[j-1] +phiOld[j]*dx

    return phi_mid[:-1],mass

def flux(nx,dx,c,phi_mid,mass):
#-----------------------------------------------------
# Flux calculation at j+1/2
#-----------------------------------------------------

    #---------------------------------------
    # Integer and remnant Courant number
    #---------------------------------------
    c_N = c.astype(int)
    c_r = c-c_N

    #---------------------------------------
    # Flux at j+1/2
    #---------------------------------------
    OUT = np.zeros_like(phi_mid)


    for i in xrange(nx):
        k= np.floor(i-c_N[i+1])%nx
        if c[i+1]>0:          # velocity u > 0
            if i>k:           # if the departure cell is at the west of predicted cell 
                OUT[i] = (phi_mid[i]*c_r[i+1]*dx+(mass[i]-mass[k]))/dx
            elif i==k:        # if the departure cell is at the position of predicted cell
                OUT[i] = phi_mid[i]*c_r[i+1]*dx/dx
            else:             # if the departure cell is at the east of predicted cell
                OUT[i] = (phi_mid[i]*c_r[i+1]*dx+mass[i]+mass[-1]-mass[k])/dx
        elif c[i+1]<0:        # velocity u < 0
            if i>k:           # if the departure cell is at the east of predicted cell  
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
# advecitve operator calculation
#-----------------------------------------------------

    #---------------------------------------
    # advective operator 
    #---------------------------------------
    AX = np.zeros_like(OUT)

    for i in xrange(nx):
        #----------------------------------------------------
        # the if statement is to update by upwind velocity
        #----------------------------------------------------
        if c[i] >= 0 and c[i+1] >0:
            AX[i] = c[i]*(OUT[i-1]/c[i]-OUT[i]/c[i+1])
        elif c[i] < 0 and c[i+1] <=0 :
            AX[i] = c[i+1]*(OUT[i]/c[i+1]-OUT[i-1]/c[i])
        else:
            AX[i] = 0
    return AX[:]

def conservative(nx,c,OUT):
#-----------------------------------------------------
# conservative operator calculation
#-----------------------------------------------------
    XC = np.zeros_like(OUT)

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
