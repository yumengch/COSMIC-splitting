import numpy as np
import matplotlib.pyplot as plt
def COSMIC(phiOld,nt,epsilon,dx,dy,J,xmin,ymin,dt,cx,cy):
	"The calculation of COSMIC splitting"
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
    for t in xrange(int(nt)):
        psi = deformational(x,y,xmin,ymin,nx,ny,Lx,Ly,t,nt,dt)
        u,v,cx,cy = streamfunction(psi,dx,dy,dt)
        #the y-direction calculation
        for i in xrange(nx):
            phi_mid[:,i],mass[:,i] = PPM(phiOld[:,i],cy[:,i],ny,epsilon,dy)
            OUT = flux(ny,dy,cy[:,i],phi_mid[:-1,i],mass[:,i])
            YC[:-1,i] = conservative(ny,cy[:,i],OUT)
            YA[:-1,i] = advective(ny,cy[:,i],OUT)
        YA[-1,:],YA[:,-1] = YA[0,:],YA[:,0]
        YC[-1,:],YC[:,-1] = YC[0,:],YC[:,0]
        #the x-direction calculation
        for j in xrange(ny):
            phi_mid[j,:], mass[j,:]= PPM(phiOld[j,:],cx[j,:],nx,epsilon,dx)
            OUT = flux(nx,dx,cx[j,:],phi_mid[j,:-1],mass[j,:])
            XC[j,:-1] = conservative(nx,cx[j,:],OUT)
            XA[j,:-1] = advective(nx,cx[j,:],OUT)
            
        XA[-1,:], XA[:,-1] = XA[0,:], XA[:,0]
        XC[-1,:], XC[:,-1] = XC[0,:], XC[:,0]
        #intermediate calculation
        phi_AX = phiOld + J*XA
        phi_AY = phiOld + J*YA

        for i in xrange(nx):
            phi_mid[:,i],mass[:,i] = PPM(phi_AX[:,i],cy[:,i],ny,epsilon,dy)
            OUT = flux(ny,dy,cy[:,i],phi_mid[:-1,i],mass[:,i])
            YC_AX[:-1,i] = conservative(ny,cy[:,i],OUT)
        YC_AX[-1,:], YC_AX[:,-1] = YC_AX[0,:], YC_AX[:,0]

        for j in xrange(ny):
            phi_mid[j,:], mass[j,:]= PPM(phi_AY[j,:],cx[j,:],nx,epsilon,dx)
            OUT = flux(nx,dx,cx[j,:],phi_mid[j,:-1],mass[j,:])
            XC_AY[j,:-1] = conservative(nx,cx[j,:],OUT)
        XC_AY[-1,:], XC_AY[:,-1]= XC_AY[0,:], XC_AY[:,0]

        #the COSMIC splitting update
        phi = phiOld+J*(0.5*(XC+XC_AY)+0.5*(YC+YC_AX))
        phi[-1,:], phi[:,-1] = phi[0,:], phi[:,0]
        phiOld = phi.copy() #update the time step
        if t == int(nt/2.)-1:
            phiMid = phiOld.copy()
        if t == int(nt/4.)-1:
            phiQ1 = phiOld.copy()
        if t == int(nt*3/4.)-1:
            phiQ2 = phiOld.copy()
        print t,phiOld[:-1,:-1].sum()
 
    return phi,phiMid,phiQ1,phiQ2

def PPM(phiOld,c,nx,epsilon,dx,eta1=20 ,eta2 = 0.05):
	"One-dimensional PPM"
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
    
    # #Calculate the second derivative of phi(centred in space)
    # d2phi[1:-1] = (phiOld[2:] - 2*phiOld[1:-1] + phiOld[:-2])/(6*dx**2)
    # # # Updated boundary condition
    # d2phi[0] = d2phi[-2]
    # d2phi[-1] = d2phi[1]
    
    # Determine the eta to weigh the discontinuity
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
    eta[:] = 0
    # dmphi[:] = 0
   #Calculate the right and left boundary value
    phi_l[1:] = (0.5*(phiOld[1:]+phiOld[:-1]) + (dmphi[:-1]-dmphi[1:])/6.)*(1-eta[1:])+(phiOld[:-1] + 0.5*dmphi[:-1])*eta[1:]
    phi_l[0] = phi_l[-2]

    phi_r[:-1] = (0.5*(phiOld[1:]+phiOld[:-1]) + (dmphi[:-1]-dmphi[1:])/6.)*(1-eta[:-1])+(phiOld[1:] - 0.5*dmphi[1:])*eta[:-1]
    phi_r[-1] =phi_r[1]

    #To avoid some extreme conditions that phi,j+1/2 is not phi_l,j and phi_r,j-1
    # for j in xrange(nx+2):
    #     if (phi_r[j]-phiOld[j])*(phiOld[j]-phi_l[j])<=0:
    #         phi_l[j] = phi_r[j] = phiOld[j]
    #     if (phi_r[j]-phi_l[j])*(phiOld[j]-0.5*(phi_l[j]+phi_r[j]))>(phi_r[j]-phi_l[j])**2/6.:
    #         phi_l[j] = 3*phiOld[j]-2*phi_r[j]
    #     if -(phi_r[j]-phi_l[j])**2/6.>(phi_r[j]-phi_l[j])*(phiOld[j]-0.5*(phi_r[j]+phi_l[j])):
    #         phi_r[j] = 3*phiOld[j]-2*phi_l[j]

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
    "calculate the mass flux in or out each grid"
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
    'calculate advective operator'
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
    'calculate conservative operator'
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