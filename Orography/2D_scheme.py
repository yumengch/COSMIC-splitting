import numpy as np
import matplotlib.pyplot as plt
def COSMIC(phiOld,nt,epsilon,dx,dy,J,initialProfile,xmin,ymin,dt,UchangewithT,cx,cy):
    # plt.ion()
    #Calculate the number of grid points
    nx = len(phiOld[0,:]) - 1
    ny = len(phiOld[:,0]) - 1
    xmax,ymax = xmin + nx*dx,ymin+ny*dy
    Lx, Ly= xmax-xmin,ymax-ymin
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
    #Define the arrays used in COSMIC
    #advective-form one-dimensional flux difference of each cell
    XA = np.zeros_like(phiOld)
    YA = np.zeros_like(phiOld)
    #advective-form one-dimensional phi update
    phi_AX = np.zeros_like(phiOld)
    phi_AY = np.zeros_like(phiOld)
    #conservative-form one-dimensional flux difference of each cell
    XC = np.zeros_like(phiOld)
    XC_AY = np.zeros_like(phiOld)
    YC = np.zeros_like(phiOld)
    YC_AX = np.zeros_like(phiOld)
    #the updated phi values
    phi = np.zeros_like(phiOld)
    #calculate the integer and remnant part of the Courant number
    phi_mid= np.zeros_like(phiOld)
    mass = np.zeros((ny,nx))
    a = np.zeros((nt,nx+1,ny+1))
    for t in xrange(int(nt)):
        #the y-direction calculation
        for i in xrange(nx):
            phi_mid[:,i],mass[:,i] = PPM((1/J[:,i])*phiOld[:,i],cy[:,i],ny,epsilon,dy)
            OUT = flux(ny,dy,cy[:,i],phi_mid[:-1,i],mass[:,i])
            YC[:-1,i] = conservative(ny,cy[:,i],OUT)
            YA[:-1,i] = advective(ny,cy[:,i],OUT)
        YA[-1,:],YA[:,-1] = YA[0,:],YA[:,0]
        YC[-1,:],YC[:,-1] = YC[0,:],YC[:,0]
        #the x-direction calculation
        for j in xrange(ny):
            phi_mid[j,:], mass[j,:]= PPM((1/J[j,:])*phiOld[j,:],cx[j,:],nx,epsilon,dx)
            OUT = flux(nx,dx,cx[j,:],phi_mid[j,:-1],mass[j,:])
            XC[j,:-1] = conservative(nx,cx[j,:],OUT)
            XA[j,:-1] = advective(nx,cx[j,:],OUT)
            
        XA[-1,:], XA[:,-1] = XA[0,:], XA[:,0]
        XC[-1,:], XC[:,-1] = XC[0,:], XC[:,0]
        #intermediate calculation
        phi_AX = phiOld + J*XA
        phi_AY = phiOld + J*YA

        for i in xrange(nx):
            phi_mid[:,i],mass[:,i] = PPM((1/J[:,i])*phi_AX[:,i],cy[:,i],ny,epsilon,dy)
            OUT = flux(ny,dy,cy[:,i],phi_mid[:-1,i],mass[:,i])
            YC_AX[:-1,i] = conservative(ny,cy[:,i],OUT)
        YC_AX[-1,:], YC_AX[:,-1] = YC_AX[0,:], YC_AX[:,0]

        for j in xrange(ny):
            phi_mid[j,:], mass[j,:]= PPM((1/J[j,:])*phi_AY[j,:],cx[j,:],nx,epsilon,dx)
            OUT = flux(nx,dx,cx[j,:],phi_mid[j,:-1],mass[j,:])
            XC_AY[j,:-1] = conservative(nx,cx[j,:],OUT)
        XC_AY[-1,:], XC_AY[:,-1]= XC_AY[0,:], XC_AY[:,0]

        #the COSMIC splitting update
        phi = phiOld+J*(0.5*(XC+XC_AY)+0.5*(YC+YC_AX))
        phi[-1,:], phi[:,-1] = phi[0,:], phi[:,0]
        phiOld = phi.copy() #update the time step
        if t == int(nt/2)-1:
            phiMid = phiOld.copy()
        print t,np.max(phiOld[:-1,:-1])#,phi_AX[:-2,:-2].sum(),YC_AX[:-2,:-2].sum(),XC[:-2,:-2].sum()
    return phi,phiMid
    # return a

def PPM(phiOld,c,nx,epsilon,dx,eta1=20 ,eta2 = 0.05):
    phiOld = np.append(phiOld, [phiOld[1]])
    dmphi = np.zeros_like(phiOld) #define the phi increment
    phi_r = np.zeros_like(phiOld) #define the right boundary
    phi_l = np.zeros_like(phiOld) #define the left boundary
    d2phi = np.zeros_like(phiOld) #define the second derivative of phi
    eta = np.zeros_like(phiOld) #define the weight of discontinuity
    phi_6 = np.zeros_like(phiOld) #define the difference between phi,j and average of phi_l and phi_r
    daj = np.zeros_like(phiOld) #define the difference between the right and left boundary
    phi_mid = np.zeros_like(phiOld) #define the midpoints
    mass = np.zeros(nx)
    #the integer part of the Courant number
    #calculate the integer and remnant part of the Courant number
    N = c.astype(int)
    dc = c-N

    #Calculate dmphi
    dmphi[1:-1] = 0.5*(phiOld[2:]-phiOld[:-2])
    # Updated boundary condition
    dmphi[0] = dmphi[-2]
    dmphi[-1] = dmphi[1]
    
    eta[:] = 0
   #Calculate the right and left boundary value
    phi_l[1:] = (0.5*(phiOld[1:]+phiOld[:-1]) + (dmphi[:-1]-dmphi[1:])/6.)*(1-eta[1:])+(phiOld[:-1] + 0.5*dmphi[:-1])*eta[1:]
    phi_l[0] = phi_l[-2]

    phi_r[:-1] = (0.5*(phiOld[1:]+phiOld[:-1]) + (dmphi[:-1]-dmphi[1:])/6.)*(1-eta[:-1])+(phiOld[1:] - 0.5*dmphi[1:])*eta[:-1]
    phi_r[-1] =phi_r[1]

    #define parabolic profile
    daj = phi_r - phi_l
    phi_6 = 6*(phiOld-0.5*(phi_l+phi_r))
    for i in xrange(nx):
        k= np.floor(i-N[i+1])%nx
        k1 = np.floor(i+1-N[i+1])%nx
        if dc[i+1]>= 0:
            phi_mid[i] = phi_r[k]-0.5*dc[i+1]*(daj[k]-(1-2*dc[i+1]/3.)*phi_6[k])
        else:
            phi_mid[i] = phi_l[k1]-0.5*dc[i+1]*(daj[k1]+(1+2*dc[i+1]/3.)*phi_6[k1])
    phi_mid[-2] = phi_mid[0]
    phi_mid[-1] = phi_mid[1]
  
    mass[0] = phiOld[0]*dx
    for j in xrange(1,nx):
        mass[j] = mass[j-1] +phiOld[j]*dx    
    return phi_mid[:-1],mass

def flux(nx,dx,c,phi_mid,mass):
    N = c.astype(int)
    dc = c-N
    OUT = np.zeros_like(phi_mid)
    for i in xrange(nx):
        k= np.floor(i-N[i+1])%nx
        if c[i+1]>0:
            if i>k:
                OUT[i] = (phi_mid[i]*dc[i+1]*dx+(mass[i]-mass[k]))/dx #+ int(N[i+1]/nx)*mass[-1]/dx
            elif i==k:
                OUT[i] = phi_mid[i]*dc[i+1]*dx/dx #+ int(N[i+1]/nx)*mass[-1]/dx
            else:
                OUT[i] = (phi_mid[i]*dc[i+1]*dx+mass[i]+mass[-1]-mass[k])/dx #+ int(N[i+1]/nx)*mass[-1]/dx
        elif c[i+1]<0:
            if i>k:
                OUT[i] = (-phi_mid[i]*dc[i+1]*dx+mass[-1]-mass[i]+mass[k])/dx #+ int(N[i+1]/nx)*mass[-1]/dx
            elif i ==k:
                OUT[i] = -phi_mid[i]*dc[i+1]*dx/dx #+ int(N[i+1]/nx)*mass[-1]/dx
            elif i<k:
                OUT[i] = (-phi_mid[i]*dc[i+1]*dx-mass[i]+mass[k])/dx #+ int(N[i+1]/nx)*mass[-1]/dx
        else:
            OUT[i] = 0
    return OUT

def advective(nx,c,OUT):
    AX = np.zeros_like(OUT)
    for i in xrange(nx):
        if c[i] > 0 and c[i+1] >0:
            AX[i] = c[i]*(OUT[i-1]/c[i]-OUT[i]/c[i+1])
        elif c[i] == 0 and c[i+1] >0:
            AX[i] = 0
        elif c[i] < 0 and c[i+1] <0 :
            AX[i] = c[i+1]*(OUT[i]/c[i+1]-OUT[i-1]/c[i])
        elif c[i] < 0 and c[i+1] ==0 :
            AX[i] = 0
        elif c[i]*c[i+1]<0:
            AX[i] = 0
        else:
            AX[i] = 0
    return AX[:]

def conservative(nx,c,OUT):
    XC = np.zeros_like(OUT)
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