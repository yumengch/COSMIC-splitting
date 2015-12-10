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
    plt.xlim([-75,75])
    plt.ylim([0,20])
    plt.figtext(.3, 0.07, "min = "+str(np.round(np.min(phiCOS-phiExact),3))+', max = '+str(np.round(np.max(phiCOS-phiExact),3)))
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    # plt.title('$\Delta x = $'+str(int(dx))+'m, $\Delta z = $'+str(int(dy))+'m, $\Delta t = $'+str(dt)+'s')
    plt.savefig('SOLUTON_ter'+str(dt)+'.pdf')