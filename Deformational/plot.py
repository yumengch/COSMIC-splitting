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
def contours_blob(X,Y,phiCOS,u,v,x,nt):
    'plot the results at the endtime'
    fig = plt.figure(1)
    plt.clf()

    levs = np.arange(-0.1,1.2,0.05)
    cs = plt.contourf(X,Y,np.round(phiCOS,6),levels = levs,cmap = my_cmap)
    # plt.colorbar(cs,orientation = 'horizontal')
    fig.tight_layout()
    plt.savefig('SOLUTON_Z'+str(nt)+'.pdf')

def contours_blob_mid(X,Y,phiCOSMid,nt):
    'plot the results at the midtime'
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(-0.1,1.2,0.05)
    plt.contourf(X,Y,np.round(phiCOSMid,6),levels = levs,cmap = my_cmap)#,colors = col,extend="both")
    # plt.colorbar(orientation = 'horizontal')
    fig.tight_layout()
    plt.savefig('Z_mid'+str(nt)+'.pdf')

def contours_blob_q(X,Y,phiquater,nt):
    'plot the results at the first a quarter of total time'
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(-0.1,1.2,0.05)
    plt.contourf(X,Y,np.round(phiquater,6),levels = levs,cmap = my_cmap)
    # plt.colorbar()
    fig.tight_layout()
    plt.savefig('Z_q'+str(nt)+'.pdf')
def contours_blob_q2(X,Y,phiquater2,nt):
    'plot the results at the last a quarter of total time'
    fig = plt.figure(1)
    plt.clf()
    levs = np.arange(-0.1,1.2,0.05)
    plt.contourf(X,Y,np.round(phiquater2,6),levels = levs,cmap = my_cmap)
    # plt.colorbar(orientation = 'horizontal')
    fig.tight_layout()
    plt.savefig('Z_q2'+str(nt)+'.pdf')