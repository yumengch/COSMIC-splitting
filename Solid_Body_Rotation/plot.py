import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# import matplotlib.animation as animation
def contours_solid(X,Y,phi,phiExact,phiOld,nt):
    hfont = {'fontname':'FreeSerif'}
    fig = plt.figure(1)
    plt.clf()
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    PPMsum = phi[0]+phi[1]+phi[2]+phi[3]+phi[4]+phi[5]
    Exactsum = phiExact[0]+phiExact[1]+phiExact[2]+phiExact[3]+phiExact[4]+phiExact[5]
    plt.contour(X,Y,np.round(phi[0],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X,Y,np.round(phi[1],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X,Y,np.round(phi[2],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X,Y,np.round(phi[3],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X,Y,np.round(phi[4],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X,Y,np.round(phi[5],3),np.arange(0.1,1.0,0.1),colors = 'k')
    plt.contour(X,Y,np.round(phiOld,3),np.arange(0.1,1.0,0.1),colors = 'k',linestyles = 'dashed')
    plt.xticks(np.arange(0,11000,2000),**hfont)
    plt.yticks(np.arange(0,11000,2000),**hfont)
    plt.gca().set_aspect('equal', adjustable='box')
    # levs = np.arange(-0.095,0.105,0.01)
    levs = np.arange(-0.95,1.,0.1)
    # color = ['#400056','#2B006B','#150081','#0000FF','#0061FF','#0092FF','#00C2FF','#00F3FF','#6DFFFF','#FFFFFF','#FFFF6D','#FFF300',\
    #     '#FFC200','#FF9200','#FF6100','#FF3100','#FF0000','#CA0600','#950B00','#601100']
    color = ['#400056','#2B006B','#150081','#0000FF','#0061FF','#0092FF','#00C2FF','#00F3FF','#6DFFFF','#FFFFFF','#FFFF6D','#FFF300',\
        '#FFC200','#FF9200','#FF6100','#FF0000','#CA0600','#950B00','#601100']
    cs = plt.contourf(X,Y,PPMsum-Exactsum,levels = levs,colors = color,extend="both")
    cs.cmap.set_under('#400056')
    cs.cmap.set_over('#601100')
    cs.set_clim(-1.0, 1.0)
    plt.figtext(.5, 0.12, "min = "+str(np.round(np.min(phi[-1]-phiOld),3))+', max = '+str(np.round(np.max(phi[-1]-phiOld),3)),horizontalalignment='center',
     verticalalignment='center',fontsize = 36, **hfont)
    # cb = plt.colorbar(cs,orientation = 'horizontal',fraction=0.046, pad=0.04)
    # cb.set_ticks(np.round(np.arange(-0.08,0.1,0.02),2))
    # cb.set_ticklabels(['%1.2f' % i for i in np.arange(-0.08,0.1,0.02)])
    # cb.set_ticks(np.round(np.arange(-0.8,0.9,0.1),1))
    # cb.set_ticklabels(['%1.1f' % i for i in np.arange(-0.8,0.9,0.1)])
    plt.tight_layout()
    matplotlib.rcParams.update({'font.size': 36})
    plt.savefig('SOLUTON_solid'+str(nt)+'.pdf', transparent=True)

def init(X,Y,u,v,phiCOS,x):
    fig = plt.figure(1)
    plt.clf()
    plt.quiver(X[::3,::3],Y[::3,::3],u[::3,::3],v[::3,::3],pivot='mid',linewidths=0.001)
    plt.contour(X,Y,phiCOS,levels = np.arange(0.1,1.0,0.1),colors = 'k')
    plt.xticks(np.arange(0,11000,1000))
    plt.yticks(np.arange(0,11000,1000))
    plt.gca().set_aspect('equal', adjustable='box') 
    plt.savefig('st_init.pdf')