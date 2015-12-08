import numpy as np
import matplotlib.pyplot as plt
def orders_solid():
    # create any array for error of different nx,dx and nt
    phierrors1 = np.zeros(3,dtype='f')
    phierrors2 = np.zeros(3,dtype='f')
    phierrors3 = np.zeros(3,dtype='f')
    deltax = np.zeros(3,dtype='f')
    deltay = np.zeros(3,dtype='f')
    
    #call function to get errors with respect to x
    # phierrors1[0],phierrors2[0],phierrors3[0], deltax[0],deltay[0]  = [0.0843,0.0641,0.0728,0.56,0.56]
    # phierrors1[1],phierrors2[1],phierrors3[1], deltax[1],deltay[1] = [0.1855,0.1560,0.1825,0.83,0.83]
    
    # phierrors1[2],phierrors2[2],phierrors3[2],deltax[2],deltay[2] = [0.7529,0.5617,0.6132,1.67,1.67]
    phierrors1[0],phierrors2[0],phierrors3[0], deltax[0],deltay[0]  = [0.086539004,0.065184053,0.072010247,0.56,0.56]
    phierrors1[1],phierrors2[1],phierrors3[1], deltax[1],deltay[1] = [0.189461052,0.159654295,0.191852967,0.83,0.83]
    
    phierrors1[2],phierrors2[2],phierrors3[2],deltax[2],deltay[2] = [0.761352988, 0.566369695, 0.61749484,1.67,1.67]
    L = 100.
    #plot l2 versus dx
    plt.figure(1)
    plt.clf()
    # plt.xlim([-1,2])
    plt.loglog(deltax, phierrors2, 'r',label='PPM with COSMIC')
    plt.loglog(deltax, phierrors2, 'r*')
    plt.loglog([deltax[0],deltax[-1]], \
        [phierrors2[0],phierrors2[0]*(deltax[-1]/deltax[0])**2], 'b-.', label='2nd order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.8*phierrors2[0],0.8*phierrors2[0]*(deltax[-1]/deltax[0])**1], 'r--', label='1st order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.6*phierrors2[0],0.6*phierrors2[0]*(deltax[-1]/deltax[0])**3], 'g:', label='3rd order')
        # plt.gca().invert_xaxis()
     # plt.xlim([-2e10,2])
    xticks,xticklabels = plt.xticks(deltax, [deltax[0],deltax[1],deltax[2]])
    # plt.xticks(deltax, [str(int(L/deltax[0]))+'X'+str(int(L/deltay[0])),str(int(L/deltax[1]))+'X'+str(int(L/deltay[1])),str(int(L/deltax[2])+1)+'X'+str(int(L/deltay[2])+1)])
    # plt.autoscale(enable=False, axis='x', tight=None)
    # xmin = (3*xticks[0] - xticks[1])/2.
# shaft half a step to the right
    # xmax = (3*xticks[-1] - xticks[-2])/2.
    plt.xlim([0.4,deltax[-1]*1.5])
    plt.xlabel('$\Delta x$')
    plt.ylabel('$\ell_2$ norm')
    plt.legend(loc= 'best')

    plt.title('Order of accuracy without limiters (Solid Body Rotation Test Case)')
    plt.savefig('order_solid.pdf')

def orders_blob():
    # create any array for error of different nx,dx and nt
    phierrors1 = np.zeros(3,dtype='f')
    phierrors2 = np.zeros(3,dtype='f')
    phierrors3 = np.zeros(3,dtype='f')
    deltax = np.zeros(3,dtype='f')
    deltay = np.zeros(3,dtype='f')
    
    #call function to get errors with respect to x

    phierrors1[0],phierrors2[0],phierrors3[0], deltax[0],deltay[0]  = [0.018648196,0.016980369,0.016491908,0.005,0.005]
    phierrors1[1],phierrors2[1],phierrors3[1], deltax[1],deltay[1] = [0.040316844,0.036124967,0.043649361,0.01,0.01]
    phierrors1[2],phierrors2[2],phierrors3[2],deltax[2],deltay[2] = [0.096870074,0.084497244,0.121724888,0.02,0.02]
    L = 1.
    #plot l2 versus dx
    plt.figure(1)
    plt.clf()
    plt.loglog(deltax, phierrors2, 'r',label='PPM with COSMIC')
    plt.loglog(deltax, phierrors2, 'r*')
    plt.loglog([deltax[0],deltax[-1]], \
        [phierrors2[0],phierrors2[0]*(deltax[-1]/deltax[0])**2], 'b-.', label='2nd order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.8*phierrors2[0],0.8*phierrors2[0]*(deltax[-1]/deltax[0])**1], 'r--', label='1st order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.6*phierrors2[0],0.6*phierrors2[0]*(deltax[-1]/deltax[0])**3], 'g:', label='3rd order')
    xticks,xticklabels = plt.xticks(deltax, [deltax[0],deltax[1],deltax[2]])
    plt.xlim([0.35*1e-2,deltax[-1]*1.5])
    # plt.gca().invert_xaxis()
    
    # plt.xticks(deltax, [str(int(L/deltax[0]))+'X'+str(int(L/deltay[0])),str(int(L/deltax[1]))+'X'+str(int(L/deltay[1])),str(int(L/deltax[2]))+'X'+str(int(L/deltay[2]))])
    
    plt.legend(loc='best')
    plt.xlabel('$\Delta x$')
    plt.ylabel('$\ell_2$ norm')
    plt.title('Order of accuracy without limiters (Stirring Test Case)')
    plt.savefig('order_Z.pdf')

def orders_terg():
    # create any array for error of different nx,dx and nt
    phierrors2 = np.zeros(6,dtype='f')
    phierrors3 = np.zeros(6,dtype='f')
    deltax = np.zeros(6,dtype='f')
    deltay = np.zeros(6,dtype='f')
    
    #call function to get errors with respect to x
    phierrors2[0],phierrors3[0], deltax[0],deltay[0]  = [0.00134041,0.00222191,125,62.5]
    phierrors2[1],phierrors3[1], deltax[1],deltay[1]  = [0.00330812,0.0045832,200,100]
    phierrors2[2],phierrors3[2], deltax[2],deltay[2]  = [0.00374586,0.00640013250,250,125]
    phierrors2[3],phierrors3[3], deltax[3],deltay[3]  = [0.0161193,0.0205361,500,250]
    phierrors2[4],phierrors3[4], deltax[4],deltay[4] = [0.0759938,0.0741328,1000,500]
    phierrors2[5],phierrors3[5],deltax[5],deltay[5] = [0.255043,0.237488,2000,1000]    #plot l2 versus dx
    plt.figure(1)
    plt.clf()
    plt.loglog(deltax, phierrors2, 'r',label='Genuinely multidimensional scheme')
    plt.loglog(deltax, phierrors2, 'r*')
    plt.loglog([deltax[0],deltax[-1]], \
        [phierrors2[0],phierrors2[0]*(deltax[-1]/deltax[0])**2], 'b-.', label='2nd order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.8*phierrors2[0],0.8*phierrors2[0]*(deltax[-1]/deltax[0])**1], 'r--', label='1st order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.6*phierrors2[0],0.6*phierrors2[0]*(deltax[-1]/deltax[0])**3], 'g:', label='3rd order')
    
    xticks,xticklabels = plt.xticks(deltax, [int(deltax[0]),int(deltax[1]),int(deltax[2]),int(deltax[3]),int(deltax[4]),int(deltax[5])])
    plt.xlim([95,deltax[-1]*1.5])
    # plt.gca().invert_xaxis()
    
    # plt.xticks(deltax, [str(150)+'X'+str(28),str(300)+'X'+str(50),str(750)+'X'+str(120)])
    
    plt.legend(loc='best')
    plt.xlabel('$\Delta x$ (m)')
    plt.ylabel('$\ell_2$ norm')
    plt.title('Order of accuracy (Genuinely multidimensional scheme)')
    plt.savefig('order_terg.pdf')

def orders_terp():
    # create any array for error of different nx,dx and nt
    phierrors2 = np.zeros(5,dtype='f')
    phierrors3 = np.zeros(5,dtype='f')
    deltax = np.zeros(5,dtype='f')
    deltay = np.zeros(5,dtype='f')
    
    #call function to get errors with respect to x
    # phierrors2[0],phierrors3[0],deltax[0],deltay[0]  = [0.00783108906154,0.0135762045364,125,62.5]
    phierrors2[0],phierrors3[0], deltax[0],deltay[0]  = [0.009132149,0.008262844,200,100]
    phierrors2[1],phierrors3[1], deltax[1],deltay[1]  = [0.010670699,0.010598636,250,125]
    phierrors2[2],phierrors3[2], deltax[2],deltay[2]  = [0.033958722,0.033875867,500,250]
    phierrors2[3],phierrors3[3], deltax[3],deltay[3] = [0.199368351,0.172507739,1000,500]
    phierrors2[4],phierrors3[4], deltax[4],deltay[4] = [0.44989719,0.429787943,2000,1000]    #plot l2 versus dx
    fig = plt.figure(1)
    plt.clf()
    plt.loglog(deltax, phierrors2, 'r',label='PPM with COSMIC')
    plt.loglog(deltax, phierrors2, 'r*')
    plt.loglog([deltax[0],deltax[-1]], \
        [phierrors2[0],phierrors2[0]*(deltax[-1]/deltax[0])**2], 'b-.', label='2nd order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.8*phierrors2[0],0.8*phierrors2[0]*(deltax[-1]/deltax[0])**1], 'r--', label='1st order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.6*phierrors2[0],0.6*phierrors2[0]*(deltax[-1]/deltax[0])**3], 'g:', label='3rd order')
    xticks,xticklabels = plt.xticks(deltax, [int(deltax[0]),int(deltax[1]),int(deltax[2]),int(deltax[3]),int(deltax[4])])
    plt.xlim([deltax[0]/1.5,deltax[-1]*1.5])
    # plt.gca().invert_xaxis()
    # fig.tight_layout()
    # plt.xticks(deltax, [str(150)+'X'+str(28),str(300)+'X'+str(50),str(750)+'X'+str(120)])
    
    plt.legend(loc='best')
    plt.xlabel('$\Delta x$ (m)')
    plt.ylabel('$\ell_2$ norm')
    plt.title('Order of accuracy (PPM with COSMIC)')
    plt.savefig('order_terp.pdf')

def orders_asym():
    # create any array for error of different nx,dx and nt
    phierrors2 = np.zeros(5,dtype='f')
    phierrors3 = np.zeros(5,dtype='f')
    deltax = np.zeros(5,dtype='f')
    deltay = np.zeros(5,dtype='f')
    
    #call function to get errors with respect to x
    # phierrors2[0],phierrors3[0],deltax[0],deltay[0]  = [0.00783108906154,0.0135762045364,125,62.5]
    phierrors2[0],phierrors3[0], deltax[0],deltay[0]  = [0.006986196,0.009865418,200,100]
    phierrors2[1],phierrors3[1], deltax[1],deltay[1]  = [0.009360691,0.012756212,250,125]
    phierrors2[2],phierrors3[2], deltax[2],deltay[2]  = [0.031425731,0.02591136,500,250]
    phierrors2[3],phierrors3[3], deltax[3],deltay[3] = [0.140452866,0.114490655,1000,500]
    phierrors2[4],phierrors3[4], deltax[4],deltay[4] = [0.317680935,0.269987787,2000,1000]    #plot l2 versus dx
    fig = plt.figure(1)
    plt.clf()
    plt.loglog(deltax, phierrors2, 'r',label='PPM with COSMIC')
    plt.loglog(deltax, phierrors2, 'r*')
    plt.loglog([deltax[0],deltax[-1]], \
        [phierrors2[0],phierrors2[0]*(deltax[-1]/deltax[0])**2], 'b-.', label='2nd order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.8*phierrors2[0],0.8*phierrors2[0]*(deltax[-1]/deltax[0])**1], 'r--', label='1st order')
    plt.loglog([deltax[0],deltax[-1]], \
        [0.6*phierrors2[0],0.6*phierrors2[0]*(deltax[-1]/deltax[0])**3], 'g:', label='3rd order')
    xticks,xticklabels = plt.xticks(deltax, [int(deltax[0]),int(deltax[1]),int(deltax[2]),int(deltax[3]),int(deltax[4])])
    plt.xlim([deltax[0]/1.5,deltax[-1]*1.5])
    # plt.gca().invert_xaxis()
    # fig.tight_layout()
    # plt.xticks(deltax, [str(150)+'X'+str(28),str(300)+'X'+str(50),str(750)+'X'+str(120)])
    
    plt.legend(loc='best')
    plt.xlabel('$\Delta x$ (m)')
    plt.ylabel('$\ell_2$ norm')
    plt.title('Order of accuracy (PPM with COSMIC)')
    plt.savefig('order_terp.pdf')
# orders_solid()
# orders_blob()
# orders_terg()
# orders_terp()
orders_asym()