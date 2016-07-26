import numpy as np
def COSMIC(phiOld, x, u, dt, nt, J):
#------------------------------------------------------------------------------------------
# Author: Yumeng Chen
# Scheme: COSMIC splitting
# Basis: Leonard, B. P., A. P. Lock, and M. K. MacVean, 1996: Conservative explicit 
#           unrestricted-time-step multidimensional constancy-preserving advection schemes. 
# Input: phiOld
# Output: phi at time nt*dt
#------------------------------------------------------------------------------------------

    #-----------------------
    # Basic grid information
    #-----------------------
    nx = len(phiOld[:]) - 1
    L = x[-1] - x[0]
    #-----------------------------------------------
    # midpoint values and mass flux at cell boundary
    #-----------------------------------------------
    phi_mid= np.zeros_like(phiOld)
    mass = np.zeros([nx])
    dx = np.zeros([nx+1])
    YC = np.zeros_like(phiOld)

    #--------------------------------
    #the updated cell-average values
    #--------------------------------
    phi = np.zeros_like(phiOld)


    

    dx[1:-1] = 0.5*(x[1:-1] - x[:-2]) + 0.5*(x[2:] - x[1:-1])
    dx[0] = 0.5*(x[1] - x[0]) + 0.5*(x[-1] - x[-2])
    dx[-1] = dx[0]

    #---------------
    # time updates
    #---------------
    for t in xrange(int(nt)):
        #----------------
        # 1D PPM updates
        #----------------
        phi_mid,mass, idx, x_depart = PPM(J*phiOld,u,nx,dx, L, x, dt)
        #---------------------------------
        # mass flux at each cell boundary
        #---------------------------------
        OUT = flux(nx,dx, phi_mid[:-1],mass, idx, u)
        #-------------------------------
        #  conservative operator updates
        #-------------------------------
        YC[:-1] = conservative(nx,u,dx,OUT)


        #---------------------------------
        # periodic boundary value updates
        #---------------------------------
        YC[-1] = YC[0]

        #---------------------------
        # phi n+1 update
        #--------------------------- 
        phi = phiOld + (1/J)*YC
        phiOld = phi.copy()
        
    return phiOld


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
    dx = np.append(dx, [dx[1]])
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
    x_depart = np.mod((np.mod(x - 0.5*dx[:-1],L) - u*dt), L)
    # print x_depart
    # print x
    # print u*dt

    idx = np.floor(np.searchsorted(x, x_depart)%nx)


    for i in xrange(nx+1):
        k = idx[i]
        # if u[i] <0:
        if np.mod(x[k] - x_depart[i],L) > 0.5*dx[k] +10e-10:
            idx[i] = int(np.floor((idx[i] -1)%nx))
        elif u[i]>0:
            if np.mod(x[k] - x_depart[i],L)  > 0.5*dx[k]-10e-10:
                idx[i] = int(np.floor((idx[i] -1)%nx))

    # for i in xrange(nx+1):
    #     dist[i] = np.fmod((x[idx[i]] + 0.5*dx[idx[i]])-x_depart[i],L)
    # print dist
    # print idx
    #---------------------------------------
    # phi increment as in PPM paper
    #---------------------------------------
    dmphi[1:-1] = 0.5*(phiOld[2:]-phiOld[:-2])
    # # --------------------------------
    # # periodic boundary value updates
    # #---------------------------------
    dmphi[0] = dmphi[-2]
    dmphi[-1] = dmphi[1]    
    # dmq[1:-1] = (dx[1:-1]/(dx[0:-2] + dx[1:-1] +
    #         dx[2:]))*(((2*dx[0:-2] + dx[1:-1]) / (dx[2:] + dx[1:-1]) )*(phiOld[2:] -
    #             phiOld[1:-1] ) + ( (2*dx[2:] + dx[1:-1]) / (dx[:-2] + dx[1:-1]) )* 
    #                 (phiOld[1:-1] - phiOld[0:-2]))
    # dmq[0] = dmq[-2]
    # dmq[-1] = dmq[1]
    #-------------
    # no limiters
    #-------------
    eta[:] = 0

    #-------------------------
    # phi at j-1/2 and j+1/2
    #-------------------------
    phi_l[1:] = (0.5*(phiOld[1:]+phiOld[:-1]) + (dmphi[:-1]-dmphi[1:])/6.)*(1-eta[1:])+(phiOld[:-1] + 0.5*dmphi[:-1])*eta[1:]
    phi_l[0] = phi_l[-2] # boundary updates

    phi_r[:-1] = (0.5*(phiOld[1:]+phiOld[:-1]) + (dmphi[:-1]-dmphi[1:])/6.)*(1-eta[:-1])+(phiOld[1:] - 0.5*dmphi[1:])*eta[:-1]
    phi_r[-1] =phi_r[1] # boundary updates
    # phi_r[1:-2] = phiOld[1:-2] + (dx[1:-2]/(dx[1:-2] + dx[2:-1]))*(phiOld[2:-1]- phiOld[1:-2])+(1/(dx[:-3] + dx[1:-2] + dx[2:-1] + dx[3:]))*((2*dx[2:-1]*dx[1:-2]/(dx[1:-2] + dx[2:-1]))*
    #             ( (dx[:-3] + dx[1:-2])/(2*dx[1:-2] + dx[2:-1]) - (dx[3:] + dx[2:-1])/(2*dx[2:-1] + dx[1:-2]) )*(phiOld[2:-1] - phiOld[1:-2])- 
    #                 dx[1:-2]*dmq[2:-1]*(dx[:-3] + dx[1:-2])/(2*dx[1:-2] + dx[2:-1]) + 
    #                   dx[2:-1]*dmq[1:-2]*(dx[2:-1] + dx[3:])/(dx[1:-2] + 2*dx[2:-1]))
    # phi_r[0] = phiOld[0] + (dx[0]/(dx[0] + dx[1]))*(phiOld[1]- phiOld[0]) + (1/(dx[-3] + dx[0] + dx[1] + dx[2]))*((2*dx[1]*dx[0]/(dx[0] + dx[1]))* 
    #          ( (dx[-3] + dx[0])/(2*dx[0] + dx[1]) - (dx[2] + dx[1])/(2*dx[1] + dx[0]) )*(phiOld[1] - phiOld[0])- 
    #             dx[0]*dmq[1]*(dx[-3] + dx[0])/(2*dx[0] + dx[1]) + 
    #               dx[1]*dmq[0]*(dx[1] + dx[2])/(dx[0] + 2*dx[1]))
    phi_r[-1] = phi_r[1]
    phi_r[-2] = phi_r[0]

    phi_l[1:] = phi_r[0:-1]
    phi_l[0] = phi_r[-3]
    #--------------------------------------------
    # piecewise parabolic subcell reconstruction
    #--------------------------------------------
    daj = phi_r - phi_l
    phi_6 = 6*(phiOld-0.5*(phi_l+phi_r))

    #--------------------------------
    # PPM update to get phi at j+1/2
    #--------------------------------
    for i in xrange(nx):
        if u[i+1] >= 0.:
            dist = np.fmod((x[idx[i+1]] + 0.5*dx[idx[i+1]])-x_depart[i+1],L)
            phi_mid[i] = (phi_r[idx[i+1]]-0.5*(dist/dx[idx[i+1]])*(daj[idx[i+1]]-(1-2*(dist/dx[idx[i+1]])/3.)*phi_6[idx[i+1]]))*dist
        else:
            dist = np.fmod(x_depart[i+1] - (x[idx[i+1]] - 0.5*dx[idx[i+1]]),L)
            phi_mid[i] = (phi_l[idx[i+1]]+0.5*(dist/dx[idx[i+1]])*(daj[idx[i+1]]+(1-2*(dist/dx[idx[i+1]])/3.)*phi_6[idx[i+1]]))*dist
    # print phi_mid
    # boundary updates
    phi_mid[-2] = phi_mid[0]
    phi_mid[-1] = phi_mid[1]

    #-----------------
    # cumulative mass 
    #-----------------   
    mass[0] = phiOld[0]*dx[0]
    for j in xrange(1,nx):
        mass[j] = mass[j-1] +phiOld[j]*dx[j]
    # mass[-1] = mass[-2] + mass[0]
    return phi_mid[:-1],mass, idx, x_depart

def flux(nx,dx, phi_mid,mass, idx, u):
#-------------------------------------
# function: Flux calculation at j+1/2
#-------------------------------------

    #---------------
    # Flux at j+1/2
    #---------------
    OUT = np.zeros_like(phi_mid)


    for i in xrange(nx):
        if u[i+1]>0:          # velocity u > 0
            k = np.floor((idx[i+1])%nx)
            if i+1>=k:           # if the departure cell is at the west of predicted cell 
                OUT[i] = (phi_mid[i]+(mass[i]-mass[k]))
            else:         # if the departure cell is at the east of predicted cell
                OUT[i] = (phi_mid[i]+mass[i]+mass[-1]-mass[k])
        elif u[i+1]<0:        # velocity u < 0            
            k = np.floor((idx[i+1]-1)%nx)
            if k>=i+1:           # if the departure cell is at the east of predicted cell  
                OUT[i] = (phi_mid[i]-mass[i]+mass[k])
            else:
                OUT[i] = (phi_mid[i]+mass[-1]-mass[i]+mass[k])
        else:               # if velcity = 0
            OUT[i] = 0
        # print OUT[i]
    return OUT

def advective(nx,u, dx, OUT):
#------------------------------------------
# function: advecitve operator calculation
#------------------------------------------

    #--------------------
    # advective operator 
    #--------------------
    AX = np.zeros_like(OUT)

    for i in xrange(nx):
        #--------------------------------------------------
        # the if statement is to update by upwind velocity
        #--------------------------------------------------
        if u[i] >= 0 and u[i+1] >0:
            AX[i] = (OUT[i-1]-(u[i]/dx[i])*OUT[i]/(u[i+1]/dx[i+1]))/dx[i]
        elif u[i] < 0 and u[i+1] <= 0 :
            AX[i] = (OUT[i]-(u[i+1]/dx[i+1])*(OUT[i-1]/(u[i]/dx[i])))/dx[i]
        else:
            AX[i] = 0
    return AX[:]

def conservative(nx,u, dx,OUT):
#---------------------------------------------
# function: conservative operator calculation
#---------------------------------------------
    XC = np.zeros_like(OUT)

    #-----------------------------------------
    # update by influx - outflux in each cell
    #-----------------------------------------
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