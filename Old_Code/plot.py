import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
def h(x):
        "The mountain height as a function of x"
        h0 = 3e3
        a = 25e3
        lam = 8e3
        return np.where(np.abs(x)<=a,
                        h0*np.cos(np.pi*x/lam)**2 * np.cos(np.pi*x*0.5/a)**2,
                        0)
def hist(X,Y,phiCOS):
    fig = plt.figure(1)
    ax = plt.gca(projection='3d') 
    # for j in xrange(len(Y[:,0])):
    #     phiCOS[j,:-1] = np.roll(phiCOS[j,:-1], (len(X[0,:])-1)/2)
    # phiCOS[:,-1] = phiCOS[:,0]
    # ax._axis3don = False
    c = ax.plot_surface(X,Y,phiCOS,rstride=1, cstride=1, cmap=cm.coolwarm,
       linewidth=0, antialiased=False)
    fig.colorbar(c)
    # ax.set_xlabel('x direction')
    # ax.set_ylabel('y direction')
    # ax.set_zlabel('$\phi$')
    # ax.set_visible(False)
    # plt.gca().invert_xaxis()
    fig.tight_layout()
    # plt.title('Smolarkiewicz Test Case Solution in Midtime')
    plt.savefig('smo.pdf')
    # plt.show()
def contours_ter(x,X,Y,phiCOS,phiOld,phiCOSMid,phiExact,phiMidExact,dx,dt,dy):
    plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    plt.contour(X/1000,Y/1000,phiOld,np.concatenate((np.arange(-2,0,0.1), np.arange(0.1,2,0.1))),colors ='k')
    plt.contour(X/1000,Y/1000,phiCOSMid,np.concatenate((np.arange(-2,0,0.1), np.arange(0.1,2,0.1))),colors ='k')
    plt.contour(X/1000,Y/1000,phiCOS,np.concatenate((np.arange(-2,0,0.1), np.arange(0.1,2,0.1))),colors ='k')
    plt.plot(x/1000,h(x)/1000,'k.')
    levs = np.arange(-0.21,0.21,0.02)
    color = ['#400056','#2B006B','#150081','#0000FF','#0031FF','#0061FF','#0092FF','#00C2FF','#00F3FF','#6DFFFF','#FFFFFF','#FFFF6D','#FFF300',\
        '#FFC200','#FF9200','#FF6100','#FF3100','#FF0000','#CA0600','#950B00','#601100']
    cs = plt.contourf(X/1000,Y/1000,phiCOS-phiExact+phiCOSMid-phiMidExact,levels = levs,colors = color,extend="both")
    cs.cmap.set_under('#400056')
    cs.cmap.set_over('#601100')
    cs.set_clim(-0.2, 0.2)
    cb = plt.colorbar(cs)
    cb.set_ticks(np.round(np.arange(-0.16,0.24,0.08),2))
    cb.set_ticklabels(['%1.2f' % i for i in np.arange(-0.16,0.24,0.08)])
    plt.xlim([-75,75])
    plt.ylim([0,15])
    plt.figtext(.3, 0.07, "min = "+str(np.round(np.min(phiCOS-phiExact),3))+', max = '+str(np.round(np.max(phiCOS-phiExact),3)))
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    # plt.title('$\Delta x = $'+str(int(dx))+'m, $\Delta z = $'+str(int(dy))+'m, $\Delta t = $'+str(dt)+'s')
    plt.savefig('SOLUTON_ter'+str(dt)+'.pdf')


def contours_solid(X,Y,phiCOS,phiOld,phiCOSMid,nt):
    
    fig = plt.figure(1)
    plt.clf()
    # plt.title('Solid Body Rotation Test Case Solution')
    plt.contour(X,Y,np.round(phiCOS,3),np.arange(-1.0,1.0,0.1))
    cs = plt.contour(X,Y,np.round(phiCOSMid,3),np.arange(-1.0,1.0,0.1))
    plt.contour(X,Y,np.round(phiOld,6),np.arange(-1.0,1.0,0.1),colors = 'k',linestyles = 'dashed')
    plt.colorbar(cs, orientation = 'horizontal')
    # fig.tight_layout()
    plt.savefig('SOLUTON_solid'+str(nt)+'.pdf')






def contours_blob(X,Y,phiCOS,u,v,x,nt):
    fig = plt.figure(1)
    plt.clf()

    levs = np.arange(-0.3,1,0.05)
    col = ['#400056','#2B006B','#150081','#0000FF','#0031FF','#0061FF','#FFFFFF','#6DFFFF','#00FF99','#BFFF00','#FFF300','#00FF00','#99FF33',\
            '#FFC200','#FF9200','#FF6100','#FF3100','#FF0000','#D90000','#CA0600','#B20000','#950B00','#601100','#260000','#000000']
    cs = plt.contourf(X,Y,np.round(phiCOS,6),levels = levs,colors = col,extend="both")
    # plt.quiver(X[::5,::5],Y[::5,::5],u[::5,::5],v[::5,::5],pivot='mid',linewidths=0.001)    
    plt.title('Solution of Stirring Test Case')
    plt.colorbar(cs)
    fig.tight_layout()
    plt.savefig('SOLUTON_Z'+str(nt)+'.pdf')

def contours_blob_mid(X,Y,phiCOSMid,nt):
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(-0.3,1,0.05)
    col = ['#400056','#2B006B','#150081','#0000FF','#0031FF','#0061FF','#FFFFFF','#6DFFFF','#00FF99','#BFFF00','#FFF300','#00FF00','#99FF33',\
            '#FFC200','#FF9200','#FF6100','#FF3100','#FF0000','#D90000','#CA0600','#B20000','#950B00','#601100','#260000','#000000']
    plt.contourf(X,Y,np.round(phiCOSMid,6),levels = levs,colors = col,extend="both")
    # plt.title('Midtime Solution of Stirring Test Case')
    # plt.colorbar(orientation = 'horizontal')
    fig.tight_layout()
    plt.savefig('Z_mid'+str(nt)+'.pdf')

def contours_blob_q(X,Y,phiquater,nt):
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(-0.3,1,0.05)
    col = ['#400056','#2B006B','#150081','#0000FF','#0031FF','#0061FF','#FFFFFF','#6DFFFF','#00FF99','#BFFF00','#FFF300','#00FF00','#99FF33',\
            '#FFC200','#FF9200','#FF6100','#FF3100','#FF0000','#D90000','#CA0600','#B20000','#950B00','#601100','#260000','#000000']
    plt.contourf(X,Y,np.round(phiquater,6),levels = levs,colors = col,extend="both")
    # plt.title('Solution of Stirring Test Case at a quater of time steps')
    # plt.colorbar()
    fig.tight_layout()
    plt.savefig('Z_q'+str(nt)+'.pdf')
def contours_blob_q2(X,Y,phiquater2,nt):
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(-0.3,1,0.05)
    col = ['#400056','#2B006B','#150081','#0000FF','#0031FF','#0061FF','#FFFFFF','#6DFFFF','#00FF99','#BFFF00','#FFF300','#00FF00','#99FF33',\
            '#FFC200','#FF9200','#FF6100','#FF3100','#FF0000','#D90000','#CA0600','#B20000','#950B00','#601100','#260000','#000000']
    plt.contourf(X,Y,np.round(phiquater2,6),levels = levs,colors = col,extend="both")
    # plt.title('Solution of Stirring Test Case at three quater of time steps')
    # plt.colorbar()
    fig.tight_layout()
    plt.savefig('Z_q2'+str(nt)+'.pdf')

















def contourconst(X,Y,phiCOS):
    plt.figure(1)
    plt.clf()
    # cs = plt.contourf(X,Y,phiCOS-1.0,levels = np.linspace(-1e-14,1e-14,16))
    cs = plt.contourf(X,Y,phiCOS)
    cb = plt.colorbar(cs)

    # plt.savefig('constancy.pdf')
    # cb.set_ticks(np.linspace(-1e-14,1e-14,17))
    # cb.set_ticklabels(['%.1e' % i for i in np.linspace(-1e-14,1e-14,17)])
    # levs = np.arange(-0.3,1,0.05)
    # col = ['#400056','#2B006B','#150081','#0000FF','#0031FF','#0061FF','#FFFFFF','#6DFFFF','#00FF99','#BFFF00','#FFF300','#00FF00','#99FF33',\
    #         '#FFC200','#FF9200','#FF6100','#FF3100','#FF0000','#D90000','#CA0600','#B20000','#950B00','#601100','#260000','#000000']
    plt.savefig('constancy.pdf')



    
def contour_ani(X,Y,a,nt):
    fig = plt.figure()
    ax = fig.gca()
    # levs = 
    col = ['#400056','#2B006B','#150081','#0000FF','#0031FF','#0061FF','#FFFFFF','#6DFFFF','#00FF99','#BFFF00','#FFF300','#00FF00','#99FF33',\
            '#FFC200','#FF9200','#FF6100','#FF3100','#FF0000','#D90000','#CA0600','#B20000','#950B00','#601100','#260000','#000000']
    def update_contour_plot(i):
      
        im = ax.contourf(X,Y,a[i])
        return im
    ani = animation.FuncAnimation(fig, update_contour_plot, frames=nt, interval= 100, repeat= True)
    ani.save('blob.mp4', fps=30)

def init(X,Y,u,v,phiCOS,x):
    fig = plt.figure(1)
    plt.clf()
    plt.quiver(X[::5,::5],Y[::5,::5],u[::5,::5],v[::5,::5],pivot='mid',linewidths=0.001)
    # plt.quiver(X[::3,::3],Y[::3,::3],u[::3,::3],v[::3,::3],pivot='mid',linewidths=0.001)
    # plt.plot(x/1000,h(x)/1000)
    # levs = np.concatenate([np.arange(-0.3,-0.1,0.1),np.arange(0.1,1.0,0.1)])
    plt.contour(X,Y,phiCOS,levels = np.arange(0.1,1.0,0.1))
    plt.colorbar()
    # plt.xlim([-75,75])
    # plt.xlabel('$x (km)$')
    # plt.ylabel('$Height (km)$')
    plt.xlabel('$x $')
    plt.ylabel('$y $')
    fig.tight_layout()
    plt.title('Solid Body Rotation Test Case')    
    plt.savefig('st_init.pdf')
    # plt.show()

def mountain(X,Y,u,v,phiCOS):
    plt.figure(1)
    plt.plot(x,h(x))
    plt.ylim ([0,25e3])
    # el = Ellipse((x[nx/2], h(x[nx/2])), x[nx/2]+2500, h(x[nx/2]+2500))
    plt.annotate('Mountain', xy=(x[nx/2], h(x[nx/2])),  xycoords='data',
                xytext=(-100, 60), textcoords='offset points',
                size=20,
                #bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",
                                # patchB=el,
                                connectionstyle="angle3,angleA=0,angleB=-90"),
                )
    plt.title('The Mountain')
    plt.savefig('terrain.pdf')
def coord():
    plt.figure(1)
    plt.contour(X,Z,Z1,np.arange(0,25e3,1e3),colors = 'k')
    plt.xlabel('Horizontal $(m)$')
    plt.ylabel('Height $(m)$')
    plt.title('Terrain Following Coordinate')
    plt.savefig('terrain_following_coordinate.pdf')

def expect():
    plt.figure(1)
    plt.plot(x/1000,h(x)/1000)
    plt.contour(X/1000,Z/1000,phi)
    plt.annotate('$t = 0$', xy = (-57, 12),size = 20)
    x0 = 50e3
    r = np.sqrt(((X-x0)/Ax)**2 + ((Z1-z0)/Az)**2)
    phi = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 
    plt.contour(X/1000,Z1/1000,phi)
    plt.annotate('$t = 10000s$', xy = (40, 12),size = 20)
    x0 = 0
    r = np.sqrt(((X-x0)/Ax)**2 + ((Z1-z0)/Az)**2)
    phi = np.where(r <= 1, np.cos(np.pi*r/2)**2, 0) 
    plt.contour(X/1000,Z1/1000,phi)
    plt.annotate('$t = 5000s$', xy = (-10, 12),size = 20)
    plt.ylim ([0,25])
    plt.xlim([-75,75])
    plt.xlabel('$x (km)$')
    plt.ylabel('$Height (m)$')
    plt.title('The Expected Results of Horizontal Advection Over Orography Test Case')
    plt.savefig('terrain_expected.pdf')
    plt.show()

def col():
    # if n<0:
    #     n = 0
    cdict = {'red':  ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.8, 1.0),
                   (0.75,1.0, 1.0),
                   (1.0, 0.4, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.9, 0.9),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.4),
                   (0.25,1.0, 1.0),
                   (0.5, 1.0, 0.8),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0)),
        
         'alpha':  ((0.0, 1.0, 1.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 1.0, 1.0))
        }
    plt.register_cmap(name='BlueRedAlpha', data=cdict)
    plt.rcParams['image.cmap'] = 'BlueRedAlpha'

    # cmap = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue', 
    #                                              colors =[(0, 0, 1), 
    #                                                       (1, 1., 1), 
    #                                                       (1, 0, 0)],
    #                                              N=len(levs)-1,
    #                                              )