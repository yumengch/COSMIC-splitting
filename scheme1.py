import numpy as np
def COSMIC(phiOld, X, Y, u, v, dt, nt, J, initialProfile, mesh, change):
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
    xmax,ymax = X[0,-1], Y[-1,0]
    xmin, ymin = X[0,0], Y[0,0]
    Lx, Ly= xmax-xmin,ymax-ymin 
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
    # print X, Y
    #-----------------------------------------------
    # midpoint values and mass flux at cell boundary
    #-----------------------------------------------
    phi_mid= np.zeros_like(phiOld)
    mass = np.zeros([ny,nx])
    dx = np.zeros([ny+1, nx+1])
    dy = np.zeros([ny+1, nx+1])
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

    # dx and dy calculation
    dx[:,:-1] = X[:, 1:] - X[:, :-1]
    dx[:,-1] = dx[:,0]
    dy[:-1, :] = Y[1:, :] - Y[:-1, :]
    dy[-1,:] = dy[0,:]
    # print dy[:,-1]
    # print dx, dy
    #---------------
    # time updates
    #---------------
    for t in xrange(int(nt)):
        #---------------------------------------------------
        # velocity updates for deformational flow test case
        #---------------------------------------------------
        # if change:
        #     cx, cy = initialProfile(x, y, xmin, ymin, nx, ny, Lx, Ly ,t, nt, dt, mesh, change)

        #-------------------------------------------------------------
        # advective operator and non-cross term conservative operator 
        # updates in y direction
        #-------------------------------------------------------------
        for i in xrange(nx):
            #------------------
            # 1D PPM updates
            #------------------

            phi_mid,mass, idx = PPM((1/J[:,i])*phiOld[:,i],v[:,i],ny,dy[:,i], Ly, Y[:,i], dt)
            #--------------------------------------
            # mass flux at each cell boundary
            #--------------------------------------
            OUT = flux(ny,dy[:,i], phi_mid[:-1],mass, idx, v[:,i])
            #----------------------------------------------
            #  adavective and conservative operator updates
            #----------------------------------------------
            YC[:-1,i] = conservative(ny,v[:,i],dy[:,i],OUT)

            YA[:-1,i] = advective(ny,v[:,i],dy[:,i],OUT)
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
            phi_mid, mass, idx = PPM((1/J[j,:])*phiOld[j,:],u[j,:],nx,dx[j,:], Lx,X[j,:], dt)
            #--------------------------------------
            # mass flux at each cell boundary
            #--------------------------------------
            OUT = flux(nx,dx[j,:],phi_mid[:-1],mass, idx, u[j,:])
            #----------------------------------------------
            #  adavective and conservative operator updates
            #----------------------------------------------    
            XC[j,:-1] = conservative(nx,u[j,:], dx[j,:], OUT)
            XA[j,:-1] = advective(nx,u[j,:], dx[j,:], OUT)
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

        # ---------------------------------------------------------------
        # conservative operator with cross-term updates in y direction
        #---------------------------------------------------------------
        for i in xrange(nx):
            #------------------
            # 1D PPM updates
            #------------------
            phi_mid, mass, idx = PPM((1/J[:,i])*phi_AX[:,i],v[:,i],ny,dy[:,i], Ly,Y[:,i], dt)
            #---------------------------------
            # mass flux at each cell boundary
            #---------------------------------
            OUT = flux(ny,dy[:,i],phi_mid[:-1],mass, idx, v[:,i])
            #-------------------------------
            # conservative operator updates
            #-------------------------------
            YC_AX[:-1,i] = conservative(ny,v[:,i], dy[:,i], OUT)
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
            phi_mid, mass, idx = PPM((1/J[j,:])*phi_AY[j,:],u[j,:],nx,dx[j,:], Lx,X[j,:], dt)
            #----------------------------------
            # mass flux at each cell boundary
            #----------------------------------
            OUT = flux(nx,dx[j,:],phi_mid[:-1],mass, idx, u[j,:])
            #----------------------------------------------
            # conservative operator updates
            #----------------------------------------------                
            XC_AY[j,:-1] = conservative(nx,u[j,:], dx[j,:], OUT)
        #-------------------------------------------------------
        # periodic boundary value updates
        #------------------------------------------------------- 
        XC_AY[-1,:], XC_AY[:,-1]= XC_AY[0,:], XC_AY[:,0]
        # print 'at ', t,' time step, the maximum of phi is ', np.max(phiOld + YC_AX)
        #-------------------------------------------------------
        # Final COSMIC splitting updates
        #------------------------------------------------------- 
        phi = phiOld+J*(0.5*(XC+XC_AY)+0.5*(YC+YC_AX))

        #-------------------------------------------------------
        # periodic boundary value updates
        #------------------------------------------------------- 
        phi[-1,:], phi[:,-1] = phi[0,:], phi[:,0]
        phiOld = phi.copy() #update the time step

        #----------------------------------------
        # print the time steps and maximum value
        #----------------------------------------
        print 'at ', t,' time step, the maximum of phi is ', np.max(phiOld)

        #-----------------------------
        # intermediate value storage
        #-----------------------------
        if initialProfile == solid: 
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
        if initialProfile == orography:
            if t == int(nt/2)-1:
                phi0 = phiOld.copy()
        if initialProfile == deform:
            if t == int(nt/5)-1:
                phi0 = phiOld.copy()
            if t == int(2*nt/5)-1:
                phi1 = phiOld.copy()
            if t == int(3*nt/5)-1:
                phi2 = phiOld.copy()
            if t == int(4*nt/5)-1:
                phi3 = phiOld.copy()

    if initialProfile == solid: 
        return [phi0, phi1, phi2, phi3, phi4, phi]
    if initialProfile == orography:
        return [phi0, phi]
    if initialProfile == deform:
        return [phi0, phi1, phi2, phi3, phi]


    

def PPM(phiOld, u, nx, dx, L, x, dt, epsilon = 0.01, eta1=20, eta2 = 0.05):
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
    dist = np.zeros_like(x)

    #---------------------------
    # departure points find
    #--------------------------
    x_depart = np.mod((np.mod(x - 0.5*dx,L) - u*dt), L)

    idx = np.rint(np.floor(np.searchsorted(x, x_depart)%nx))

    for i in xrange(nx+1):
        k = idx[i]
        if np.mod(x[k] - x_depart[i],L) > 0.5*dx[k] +10e-10:
            idx[i] = int(np.floor((idx[i] -1)%nx))


    #---------------------------------------
    # phi increment as in PPM paper
    #---------------------------------------
    dmphi[1:-1] = 0.5*(phiOld[2:]-phiOld[:-2])
    # -------------------------------------------------------
    # periodic boundary value updates
    #-------------------------------------------------------
    dmphi[0] = dmphi[-2]
    dmphi[-1] = dmphi[1]    

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
    # piecewise parabolic subcell reconstruction
    #-----------------------------------------------------------
    daj = phi_r - phi_l
    phi_6 = 6*(phiOld-0.5*(phi_l+phi_r))

    #------------------------------------------------------------
    # PPM update to get phi at j+1/2
    #-----------------------------------------------------------
    for i in xrange(nx):
        if u[i+1] >= 0.:
            dist = np.mod((x[idx[i]] + 0.5*dx[idx[i]])-x_depart[i],L)
            phi_mid[i] = (phi_r[idx[i+1]]-0.5*(dist/dx[idx[i+1]])*(daj[idx[i+1]]-(1-2*(dist/dx[idx[i+1]])/3.)*phi_6[idx[i+1]]))*dist
        else:
            dist = np.fmod(x_depart[i+1] - (x[idx[i+1]] - 0.5*dx[idx[i+1]]),L)
            phi_mid[i] = (phi_l[idx[i+1]]+0.5*(dist/dx[idx[i+1]])*(daj[idx[i+1]]+(1-2*(dist/dx[idx[i+1]])/3.)*phi_6[idx[i+1]]))*dist

    # boundary updates
    phi_mid[-2] = phi_mid[0]
    phi_mid[-1] = phi_mid[1]
    # print np.max(phi_mid)w
    #-------------------
    # cumulative mass 
    #-------------------    
    mass[0] = phiOld[0]*dx[0]
    for j in xrange(1,nx):
        mass[j] = mass[j-1] +phiOld[j]*dx[j]
    # mass[-1] = mass[-2] + mass[0]
    return phi_mid[:-1],mass, idx

def flux(nx,dx, phi_mid,mass, idx, u):
#-----------------------------------------------------
# function: Flux calculation at j+1/2
#-----------------------------------------------------

    #---------------------------------------
    # Flux at j+1/2
    #---------------------------------------
    OUT = np.zeros_like(phi_mid)


    for i in xrange(nx):
        if u[i+1]>0:          # velocity u > 0
            k = idx[i+1]
            if i>k:           # if the departure cell is at the west of predicted cell 

                OUT[i] = (phi_mid[i]+(mass[i]-mass[k]))
            elif i==k:        # if the departure cell is at the position of predicted cell
                OUT[i] = phi_mid[i]
            else:         # if the departure cell is at the east of predicted cell
                OUT[i] = (phi_mid[i]+mass[i]+mass[-1]-mass[k])
        elif u[i+1]<0:        # velocity u < 0            
            k = int(np.floor((idx[i+1]-1)%nx))
            if k>=i+1:           # if the departure cell is at the east of predicted cell  
                OUT[i] = (phi_mid[i]-mass[i]+mass[k])
            else:
                OUT[i] = (phi_mid[i]+mass[-1]-mass[i]+mass[k])

        else:               # if velcity = 0
            OUT[i] = 0

    return OUT

def advective(nx,u, dx, OUT):
#-----------------------------------------------------
# function: advecitve operator calculation
#-----------------------------------------------------

    #---------------------------------------
    # advective operator 
    #---------------------------------------
    AX = np.zeros_like(OUT)

    for i in xrange(nx):
        #----------------------------------------------------
        # the if statement is to update by upwind velocity
        #----------------------------------------------------
        if u[i] >= 0 and u[i+1] >0:
            AX[i] = (OUT[i-1]-(u[i]/dx[i])*OUT[i]/(u[i+1]/dx[i+1]))/dx[i]
        elif u[i] < 0 and u[i+1] <= 0 :
            AX[i] = (OUT[i]-(u[i+1]/dx[i+1])*(OUT[i-1]/(u[i]/dx[i])))/dx[i]
        else:
            AX[i] = 0
    return AX[:]

def conservative(nx,u, dx,OUT):
#-----------------------------------------------------
# function: conservative operator calculation
#-----------------------------------------------------
    XC = np.zeros_like(OUT)

    #-----------------------------------------------------
    # update by influx - outflux in each cell
    #-----------------------------------------------------
    for i in xrange(nx):
        if u[i+1]<=0 and u[i]<=0:
            XC[i] = (OUT[i]-OUT[i-1])/dx[i]
        elif u[i+1]>=0 and u[i]>=0:
            XC[i] = (OUT[i-1]-OUT[i])/dx[i]
        elif u[i+1]>=0 and u[i]<=0:
            XC[i] = -(OUT[i]+OUT[i-1])/dx[i]
        elif u[i+1]<=0 and u[i]>=0:
            XC[i] = (OUT[i]+OUT[i-1])/dx[i]
    return XC[:]
