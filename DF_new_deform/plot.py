import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
f = open('wh-bl-gr-ye-re.rgb','r')
f.readline()
f.readline()
f.readline()
colors = []
for i in xrange(199):
    a = [int(x) for x in f.readline().split()]
    colors.append(a)
f.close()
def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap
### Create a list of RGB tuples
### Call the function make_cmap which returns your colormap
my_cmap = make_cmap(colors, bit=True)
### Use your colormap
def contours_END(X,Y,phiCOS,u,v,x,nt):
    'plot the results at the endtime'
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(0.,1.05,0.05)
    cs = plt.contourf(X,Y,phiCOS,levels = levs,cmap = my_cmap,extend="both",vmin = -0.05,vmax = 1.05)
    plt.contour(X,Y,np.round(phiCOS,3),levels = np.arange(0.05,1.05,0.05),colors ='k',linewidths=0.2)
    plt.xticks(np.arange(-180,210,30))
    plt.yticks(np.arange(-90,120,30))
    plt.figtext(.5, 0.16, "min = "+str(np.round(np.min(phiCOS),3))+', max = '+str(np.round(np.max(phiCOS),3)),horizontalalignment='center',
     verticalalignment='center')
    plt.axis('scaled')
    fig.tight_layout()
    plt.savefig('deformingAdvection_onOrthogW_c1_480x240_5_T.pdf')

def contours_TT1(X,Y,phiT1,nt):
    'plot the results at the midtime'
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(0.,1.05,0.05)
    cs = plt.contourf(X,Y,phiT1,levels = levs,cmap = my_cmap,extend="both",vmin = -0.05,vmax = 1.05)
    plt.contour(X,Y,np.round(phiT1,3),levels = np.arange(0.05,1.05,0.05),colors ='k',linewidths=0.2)
    plt.xticks(np.arange(-180,210,30))
    plt.yticks(np.arange(-90,120,30))
    plt.figtext(.5, 0.16, "min = "+str(np.round(np.min(phiT1),3))+', max = '+str(np.round(np.max(phiT1),3)),horizontalalignment='center',
     verticalalignment='center')
    plt.axis('scaled')
    fig.tight_layout()
    plt.savefig('deformingAdvection_onOrthogW_c1_480x240_1_T.pdf')

def contours_T2(X,Y,phiT2,nt):
    'plot the results at the first a quarter of total time'
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(0.,1.05,0.05)
    cs = plt.contourf(X,Y,phiT2,levels = levs,cmap = my_cmap,extend="both",vmin = -0.05,vmax = 1.05)
    plt.contour(X,Y,np.round(phiT2,3),levels = np.arange(0.05,1.05,0.05),colors ='k',linewidths=0.2)
    plt.xticks(np.arange(-180,210,30))
    plt.yticks(np.arange(-90,120,30))
    plt.figtext(.5, 0.16, "min = "+str(np.round(np.min(phiT2),3))+', max = '+str(np.round(np.max(phiT2),3)),horizontalalignment='center',
     verticalalignment='center')
    plt.axis('scaled')
    fig.tight_layout()
    plt.savefig('deformingAdvection_onOrthogW_c1_480x240_2_T.pdf')
def contours_T3(X,Y,phiT3,nt):
    'plot the results at the last a quarter of total time'
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(0.,1.05,0.05)
    cs = plt.contourf(X,Y,phiT3,levels = levs,cmap = my_cmap,extend="both",vmin = -0.05,vmax = 1.05)
    plt.contour(X,Y,np.round(phiT3,3),levels = np.arange(0.05,1.05,0.05),colors ='k',linewidths=0.2)
    plt.xticks(np.arange(-180,210,30))
    plt.yticks(np.arange(-90,120,30))
    plt.figtext(.5, 0.16, "min = "+str(np.round(np.min(phiT3),3))+', max = '+str(np.round(np.max(phiT3),3)),horizontalalignment='center',
     verticalalignment='center')
    plt.axis('scaled')
    fig.tight_layout()
    plt.savefig('deformingAdvection_onOrthogW_c1_480x240_3_T.pdf')
def contours_T4(X,Y,phiT4,nt):
    'plot the results at the last a quarter of total time'
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(0.,1.05,0.05)
    cs = plt.contourf(X,Y,phiT4,levels = levs,cmap = my_cmap,extend="both",vmin = -0.05,vmax = 1.05)
    plt.contour(X,Y,np.round(phiT4,3),levels = np.arange(0.05,1.05,0.05),colors ='k',linewidths=0.2)
    plt.xticks(np.arange(-180,210,30))
    plt.yticks(np.arange(-90,120,30))
    plt.figtext(.5, 0.16, "min = "+str(np.round(np.min(phiT4),3))+', max = '+str(np.round(np.max(phiT4),3)),horizontalalignment='center',
     verticalalignment='center')
    plt.axis('scaled')
    fig.tight_layout()
    plt.savefig('deformingAdvection_onOrthogW_c1_480x240_4_T.pdf')
def contours_Old(X,Y,phiOld,nt):
    'plot the results at the last a quarter of total time'
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(0.,1.05,0.05)
    cs = plt.contourf(X,Y,phiOld,levels = levs,cmap = my_cmap,extend="both",vmin = -0.05,vmax = 1.05)
    plt.contour(X,Y,np.round(phiOld,3),levels = np.arange(0.05,1.05,0.05),colors ='k',linewidths=0.2)
    plt.xticks(np.arange(-180,210,30))
    plt.yticks(np.arange(-90,120,30))
    plt.figtext(.5, 0.16, "min = "+str(np.round(np.min(phiOld),3))+', max = '+str(np.round(np.max(phiOld),3)),horizontalalignment='center',
     verticalalignment='center')
    plt.axis('scaled')
    fig.tight_layout()
    plt.savefig('deformingAdvection_onOrthogW_c1_480x240_0_T.pdf')
    # plt.show()