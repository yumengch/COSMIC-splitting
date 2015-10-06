import numpy as np
import matplotlib.pyplot as plt
def analytical(X,Y,dx,dy,nt,dt,nx,ny,u0=10):
    x0 = -50e3+nt*dt*u0
    Ax = 25e3
    Az = 3e3
    z0 = 12e3
    r = np.sqrt(((X-x0)/Ax)**2 + ((Y-z0)/Az)**2)
    phi = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 
    # plt.figure(1)
    # plt.contour(X/1000,Y/1000,phi,np.concatenate((np.arange(-2,0,0.1), np.arange(0.1,2,0.1))))
    return phi

def analytical2(X,Y,dx,dy,nt,dt,nx,ny,u0=10):
    A = 8
    r = np.sqrt((30)**2+(30)**2)
    theta0 = np.arctan((30)/(30))
    x0 = 0.5*nx*dx+r*np.cos(theta0+2*A*nt*dt)
    y0 = 0.5*ny*dy+ r*np.sin(theta0+2*A*nt*dt)
    print x0,y0
    phi = np.exp(- 0.5*(((X-x0) / (3))**2 + ((Y-y0) / (3))**2))
    return phi

def norms(phiExact, phi, dx, dy):
	#remove the wrap-around points
    phi = phi[0:-1,0:-1]
    phiExact = phiExact[0:-1,0:-1]
    
    # calculate the error and the error norms
    phiError = phi - phiExact
    l1 = np.sum(np.abs(phiError))/np.sum(np.abs(phiExact))
    # correct the following two lines
    l2 = np.sqrt(np.sum(phiError**2))/np.sqrt(np.sum(phiExact**2))
    linf = np.max(np.abs(phiError))

    return [l1,l2,linf]

def orders_solid():
    # create any array for error of different nx,dx and nt
    phierrors1 = np.zeros(3,dtype='f')
    phierrors2 = np.zeros(3,dtype='f')
    phierrors3 = np.zeros(3,dtype='f')
    deltax = np.zeros(3,dtype='f')
    deltay = np.zeros(3,dtype='f')
    
    #call function to get errors with respect to x

    phierrors1[2],phierrors2[2],phierrors3[2], deltax[2],deltay[2]  = advection(initialProfile = solid, xmin = 0., xmax = 100., ymin = 0., ymax = 100., epsilon = 0.01, dt= np.pi/720, nx = 256, ny = 256, nt =90)
    phierrors1[1],phierrors2[1],phierrors3[1], deltax[1],deltay[1] = advection(initialProfile = solid, xmin = 0., xmax = 100., ymin = 0., ymax = 100., epsilon = 0.01, dt= np.pi/480, nx = 128, ny = 128, nt =60)
    
    phierrors1[0],phierrors2[0],phierrors3[0],deltax[0],deltay[0] = advection(initialProfile = solid, xmin = 0., xmax = 100., ymin = 0., ymax = 100., epsilon = 0.01, dt= np.pi/240, nx = 64, ny = 64, nt =30)
    L = 100.
    #plot l2 versus dx
    plt.figure(1)
    plt.clf()
    plt.loglog(deltax, phierrors2, 'r',label='COSMIC with PPM')
    plt.loglog([deltax[0],deltax[-1]], \
        [phierrors2[0],phierrors2[0]*(deltax[-1]/deltax[0])**2], 'b--', label='2nd order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.8*phierrors2[0],0.8*phierrors2[0]*(deltax[-1]/deltax[0])**1], 'r--', label='1st order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.6*phierrors2[0],0.6*phierrors2[0]*(deltax[-1]/deltax[0])**3], 'g--', label='3rd order')
    plt.xlim([deltax[0],deltax[-1]])
    # plt.gca().invert_xaxis()
    
    plt.xticks(deltax, [str(int(L/deltax[0]))+'X'+str(int(L/deltay[0])),str(int(L/deltax[1]))+'X'+str(int(L/deltay[1])),str(int(L/deltax[2])+1)+'X'+str(int(L/deltay[2])+1)])
    
    plt.legend(loc='best')
    plt.xlabel('Gridpoint number')
    plt.ylabel('$\ell_2$ norm')
    plt.title('Order of accuracy (Solid Body Rotation Test Case)')
    plt.savefig('order_solid.pdf')

def orders_blob():
    # create any array for error of different nx,dx and nt
    phierrors1 = np.zeros(3,dtype='f')
    phierrors2 = np.zeros(3,dtype='f')
    phierrors3 = np.zeros(3,dtype='f')
    deltax = np.zeros(3,dtype='f')
    deltay = np.zeros(3,dtype='f')
    
    #call function to get errors with respect to x

    phierrors1[0],phierrors2[0],phierrors3[0], deltax[0],deltay[0]  = advection(initialProfile = stirring, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.0021875, nx = 200, ny = 200, nt =600)
    phierrors1[1],phierrors2[1],phierrors3[1], deltax[1],deltay[1] = advection(initialProfile = stirring, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.004375, nx = 100, ny = 100, nt =300)
    phierrors1[2],phierrors2[2],phierrors3[2],deltax[2],deltay[2] = advection(initialProfile = stirring, xmin = 0., xmax = 1., ymin = 0., ymax = 1., epsilon = 0.01, dt= 0.00875, nx = 50, ny = 50, nt =150)
    L = 1.
    #plot l2 versus dx
    plt.figure(1)
    plt.clf()
    plt.loglog(deltax, phierrors2, 'r',label='COSMIC with PPM')
    plt.loglog([deltax[0],deltax[-1]], \
        [phierrors2[0],phierrors2[0]*(deltax[-1]/deltax[0])**2], 'b--', label='2nd order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.8*phierrors2[0],0.8*phierrors2[0]*(deltax[-1]/deltax[0])**1], 'r--', label='1st order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.6*phierrors2[0],0.6*phierrors2[0]*(deltax[-1]/deltax[0])**3], 'g--', label='3rd order')
    plt.xlim([deltax[0],deltax[-1]])
    xticks,xticklabels = plt.xticks(deltax, [deltax[0],deltax[1],deltax[2]])
    plt.xlim([deltax[0]/1.5,deltax[-1]*1.5])
    
    plt.legend(loc='best')
    plt.xlabel('Gridpoint number')
    plt.ylabel('$\ell_2$ norm')
    plt.title('Order of accuracy (Stirring Test Case)')
    plt.savefig('order_Z.pdf')

def orders_ter():
    # create any array for error of different nx,dx and nt
    phierrors1 = np.zeros(3,dtype='f')
    phierrors2 = np.zeros(3,dtype='f')
    phierrors3 = np.zeros(3,dtype='f')
    deltax = np.zeros(3,dtype='f')
    deltay = np.zeros(3,dtype='f')
    
    #call function to get errors with respect to x

    phierrors1[2],phierrors2[2],phierrors3[2], deltax[2],deltay[2]  = advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 10, nx =750, ny = 120, nt =1000)
    phierrors1[1],phierrors2[1],phierrors3[1], deltax[1],deltay[1] = advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 25, nx =300, ny = 50, nt =400)
    
    phierrors1[0],phierrors2[0],phierrors3[0],deltax[0],deltay[0] = advection(initialProfile = terrain, xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, epsilon = 0.01, dt= 50, nx =150, ny = 28, nt =200)
    #plot l2 versus dx
    plt.figure(1)
    plt.clf()
    plt.loglog(deltax, phierrors2, 'r',label='COSMIC with PPM')
    plt.loglog([deltax[0],deltax[-1]], \
        [phierrors2[0],phierrors2[0]*(deltax[-1]/deltax[0])**2], 'b--', label='2nd order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.8*phierrors2[0],0.8*phierrors2[0]*(deltax[-1]/deltax[0])**1], 'r--', label='1st order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.6*phierrors2[0],0.6*phierrors2[0]*(deltax[-1]/deltax[0])**3], 'g--', label='3rd order')
    plt.xlim([deltax[0],deltax[-1]])
    # plt.gca().invert_xaxis()
    
    plt.xticks(deltax, [str(150)+'X'+str(28),str(300)+'X'+str(50),str(750)+'X'+str(120)])
    
    plt.legend(loc='best')
    plt.xlabel('Gridpoint number')
    plt.ylabel('$\ell_2$ norm')
    plt.title('Order of accuracy (Horizontal Advection Over Orography Test Case)')
    plt.savefig('order_ter3.pdf')