import numpy as np
def COSMIC(phiOld, cx, cy, u, v, X_cntr, X_edge, Y_edge, Y, dt, nt, J, J_p, initialProfile, mesh, change):
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
    nx, ny = len(phiOld[0,:]), len(phiOld[:,0])
    xmin, ymin = np.min(X_edge[0,:]), np.min(Y_edge[:,0])
    dx, dy = X_edge[0,1] - X_edge[0,0], Y_edge[1,0] - Y_edge[0,0]
    xmax, ymax = xmin + nx*dy, ymin + ny*dy
    Lx, Ly = xmax - xmin, ymax - ymin 

    print nx, ny, xmin, ymin, dx, dy, xmax, ymax, Lx, Ly
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
    idx_x = np.zeros_like(phiOld).astype(int)
    idx_y = np.zeros_like(phiOld).astype(int)
    r_x = np.zeros_like(phiOld)
    r_y = np.zeros_like(phiOld)
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

    
    #---------------------------
    # find the departure points
    #---------------------------
    # for j in xrange(ny):
    #     idx_x[j,:], r_x[j,:] = departure_x(X_edge[j,:], xmin, xmax, u[j,:], dx, dt, Lx)
    # for i in xrange(nx):
    #     idx_y[:,i], r_y[:,i] = departure_y(X_edge[:,i], Y_edge[:,i], Y[:,i], ymin, ymax, v[:,i], dy, dt, Lx, Ly, mesh)
    #---------------
    # time updates
    #---------------
    for t in xrange(int(nt)):
        #---------------------------------------------------
        # velocity updates for deformational flow test case
        #---------------------------------------------------
        if change:
            x_edge,y_edge = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
            x_cntr, y_cntr = x_edge[:-1] + 0.5*dx, y_edge[:-1] + 0.5*dy
            u, v, cx, cy = initialProfile(x_edge, y_edge, x_cntr, y_cntr, t, nt, dt, mesh, change)
        #---------------------------
        # find the departure points
        #---------------------------
        for j in xrange(ny):
            idx_x[j,:], r_x[j,:], cx[j,:] = departure_x(X_edge[j,:], xmin, xmax, u[j,:], dx, dt, Lx)
        for i in xrange(nx):
            idx_y[:,i], r_y[:,i], cy[:,i] = departure_y(X_edge[:,i], Y_edge[:,i], Y[:,i], ymin, ymax, v[:,i], dy, dt, Lx, Ly, mesh)


        #-------------------------------------------------------------
        # advective operator and non-cross term conservative operator 
        # updates in y direction
        #-------------------------------------------------------------
        for i in xrange(nx):
            #------------------
            # 1D PPM updates
            #-----------------
            phi_mid[:,i],mass[:,i] = PPM((1/J[:,i])*phiOld[:,i],J[:,i]*cy[:,i],ny,dy, idx_y[:,i], J[:,i]*r_y[:,i])
            #--------------------------------------
            # mass flux at each cell boundary
            #--------------------------------------
            OUT = flux(ny,dy,J[:,i]*cy[:,i],phi_mid[:,i],mass[:,i], idx_y[:,i], J[:,i]*r_y[:,i])
            #----------------------------------------------
            #  adavective and conservative operator updates
            #----------------------------------------------
            YC[:,i] = conservative(ny,J[:,i]*cy[:,i],OUT)
            YA[:,i] = conservative(ny,J[:,i]*cy[:,i],OUT)


        #----------------------------------------------------------
        # advective operator and non-cross term conservative operator 
        # updates in x direction
        #----------------------------------------------------------
        for j in xrange(ny):
            #------------------
            # 1D PPM updates
            #------------------
            phi_mid[j,:], mass[j,:]= PPM((1/J[j,:])*phiOld[j,:],J[j,:]*cx[j,:],nx,dx, idx_x[j,:], J[j,:]*r_x[j,:])
            #--------------------------------------
            # mass flux at each cell boundary
            #--------------------------------------
            OUT = flux(nx,dx,J[j,:]*cx[j,:],phi_mid[j,:],mass[j,:], idx_x[j,:], J[j,:]*r_x[j,:])
            #----------------------------------------------
            #  adavective and conservative operator updates
            #----------------------------------------------    
            XC[j,:] = conservative(nx,J[j,:]*cx[j,:],OUT)
            XA[j,:] = conservative(nx,J[j,:]*cx[j,:],OUT)

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
            phi_mid[:,i],mass[:,i] = PPM((1/J[:,i])*phi_AX[:,i],J[:,i]*cy[:,i],ny,dy, idx_y[:,i], J[:,i]*r_y[:,i])
            #---------------------------------
            # mass flux at each cell boundary
            #---------------------------------
            OUT = flux(ny,dy,J[:,i]*cy[:,i],phi_mid[:,i],mass[:,i], idx_y[:,i], J[:,i]*r_y[:,i])
            #-------------------------------
            # conservative operator updates
            #-------------------------------
            YC_AX[:,i] = conservative(ny,J[:,i]*cy[:,i],OUT)


        #---------------------------------------------------------------
        # conservative operator with cross-term updates in x direction
        #---------------------------------------------------------------
        for j in xrange(ny):
            #------------------
            # 1D PPM updates
            #------------------
            phi_mid[j,:], mass[j,:]= PPM((1/J[j,:])*phi_AY[j,:],J[j,:]*cx[j,:],nx,dx, idx_x[j,:], J[j,:]*r_x[j,:])
            #----------------------------------
            # mass flux at each cell boundary
            #----------------------------------
            OUT = flux(nx,dx,J[j,:]*cx[j,:],phi_mid[j,:],mass[j,:], idx_x[j,:],J[j,:]*r_x[j,:])
            #----------------------------------------------
            # conservative operator updates
            #----------------------------------------------                
            XC_AY[j,:] = conservative(nx,J[j,:]*cx[j,:],OUT)


        #-------------------------------------------------------
        # Final COSMIC splitting updates
        #------------------------------------------------------- 
        phi = phiOld+J*(0.5*(XC+XC_AY)+0.5*(YC+YC_AX))

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
    # return phi

# def departure_y(X_cntr, X_edge, Y_edge, Y, ymin, ymax, v, dy, dt, Lx, Ly, mesh):
#     r = np.zeros_like(v)
#     k = np.zeros_like(v)
#     idx = np.zeros_like(v)
#     y_depart = Y - v*dt
#     for i in xrange(len(y_depart)):
#         while y_depart[i] <= ymin:
#             y_depart[i] += Ly
#         while y_depart[i] >= ymax or abs(y_depart[i] - ymax) <= 10e-8:
#             y_depart[i] -= Ly
#     # for i in xrange(len(y_depart)):
#     #     while abs(y_depart[i] - ymin) <= 10e-8 and (v[i] > 0.0 or abs(v[i])< 10e-8):
#     #         y_depart[i] += Ly
#     #     while  and v[i] < -10e-8:
#     #         y_depart[i] -= Ly
#     #-----------------------
#     # transform back to 
#     # computational domain
#     #-----------------------
#     if mesh == 'quad':
#         y_depart = phys_to_compt_SB(y_depart, Ly, f_quad(X_edge,Lx, Ly))
#     elif mesh == 'V':
#         y_depart = phys_to_compt_SB(y_depart, Ly, f_V(X_edge,Lx, Ly))
#     elif mesh == 'high':
#         y_depart = phys_to_compt_oro(X_edge,y_depart,ymax, 6e3)
#     elif mesh == 'low':
#         y_depart = phys_to_compt_oro(X_edge,y_depart,ymax, 3e3)
#     elif mesh == 'W':
#         y_depart = phys_to_compt_deform(X_edge, y_depart, Lx, ymin, ymax)
#     else:
#         y_depart = y_depart

#     k = y_depart - ymin
#     k = np.where(abs(k)<= 10e-6, 0., k)
#     for i in xrange(len(v)):
#         if v[i] >= 0.0:
#             idx[i] = np.ceil(k[i]/dy).astype(int)-1
#         else:
#             idx[i] = np.floor(k[i]/dy).astype(int)
#     # idx = np.where(v[i]>=0.0, np.ceil((y_depart - ymin)/dy).astype(int)-1 ,np.floor((y_depart - ymin)/dy).astype(int))
#     # print y_depart[0], idx[0], v[0], k[0], v[0]>=0.0
#     Y_edge = np.append(Y_edge, ymax)
#     # print len(Y_edge), len(Y), len(y_depart)
#     # print np.max(idx)
#     for i in xrange(len(Y)):
#         # print idx[i], y_depart[i]
#         if v[i] >= 0.:
#             r[i] = (Y_edge[idx[i]+1] - y_depart[i])/dy
#         else:
#             r[i] = -(y_depart[i] - Y_edge[idx[i]])/dy
#         # print i, r[i]
#     return idx, r

# def departure_x(X, xmin, xmax, u, dx, dt, L):
#     r = np.zeros_like(u)
#     x_depart = X - u*dt

#     for i in xrange(len(x_depart)):
#         while x_depart[i] <= xmin:
#             x_depart[i] += L
#         while x_depart[i] >= xmax:
#             x_depart[i] -= L

#     idx = np.where(u[i]>=0.0, np.ceil((x_depart - xmin)/dx).astype(int)-1 ,np.floor((x_depart - xmin)/dx).astype(int))

#     X = np.append(X, xmax)
#     for i in xrange(len(X)-1):
#         if u[i] > 0:
#             r[i] = ( X[idx[i]+1] - x_depart[i])/dx
#         else:
#             r[i] = -(x_depart[i] - X[idx[i]])/dx

#     return idx, r
    
def departure_y(X_edge, Y_edge, Y, ymin, ymax, v, dy, dt, Lx, Ly, mesh):
    r = np.zeros_like(v)
    c = np.zeros_like(v)
    idx = np.zeros_like(v)
    y_length = np.zeros_like(v)
    ny = len(Y)

    y_depart = Y - v*dt
    y_length[:] = y_depart[:]
    for i in xrange(len(y_depart)):
        while y_depart[i] < ymin:
            y_depart[i] += Ly
        while y_depart[i] > ymax or abs(y_depart[i] - ymax) <= 10e-8:
            y_depart[i] -= Ly

    #-----------------------
    # transform back to 
    # computational domain
    #-----------------------
    if mesh == 'quad':
        y_depart = phys_to_compt_SB(y_depart, Ly, f_quad(X_edge,Lx, Ly))
    elif mesh == 'V':
        y_depart = phys_to_compt_SB(y_depart, Ly, f_V(X_edge,Lx, Ly))
    elif mesh == 'high':
        y_depart = phys_to_compt_oro(X_edge,y_depart,ymax, 6e3)
    elif mesh == 'low':
        y_depart = phys_to_compt_oro(X_edge,y_depart,ymax, 3e3)
    elif mesh == 'W':
        y_depart = phys_to_compt_deform(X_edge, y_depart, Lx, ymin, ymax)
    else:
        y_depart = y_depart

    c = (Y_edge - y_depart)/dy
    c = np.where(np.absolute(c) < 10e-12, 0.0, c)
    c = np.where( (v> 0) &  ((c <0) | (y_length > ymin)), (ymax - y_depart + Y_edge - ymin)/dy, c )
    c = np.where( (v< 0) &  ((c > 0) | (y_length > ymax)), -(ymax + y_depart - Y_edge - ymin)/dy, c )
    N = c.astype(int)
    r = c - N

    for i in xrange(len(v)):
        if c[i] >= 0.0:
            idx[i] = int(np.floor(i-1-N[i])%ny)
        else:
            idx[i] = int(np.floor(i-N[i])%ny)


    return idx, r, c

def departure_x(X, xmin, xmax, u, dx, dt, L):
    r = np.zeros_like(u)
    c = np.zeros_like(u)
    idx = np.zeros_like(u)
    x_length = np.zeros_like(u)

    nx = len(X)

    x_depart = X - u*dt
    x_length[:] = x_depart[:]
    for i in xrange(len(x_depart)):
        while x_depart[i] <= xmin:
            x_depart[i] += L
        while x_depart[i] >= xmax:
            x_depart[i] -= L


    c = (X - x_depart)/dx
    c = np.where(np.absolute(c) < 10e-12, 0.0, c)
    c = np.where( (u > 0) & ((c <0) | (x_length < xmin)), (xmax - x_depart + X - xmin)/dx, c )
    c = np.where( (u < 0) & ((c >0) | (x_length > xmax)), -(xmax + x_depart - X - xmin)/dx, c )

    N = c.astype(int)
    r = c - N

    for i in xrange(len(u)):
        if c[i] >= 0.0:
            idx[i] = int(np.floor(i-1-N[i])%nx)
        else:
            idx[i] = int(np.floor(i-N[i])%nx)

    return idx, r, c

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
            elif i< k:       # if the departure cell is at the west of predicted cell 
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


