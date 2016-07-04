import numpy as np
import matplotlib.pyplot as plt
def orders_asym():
    # create any array for error of different nx,dx and nt
    phierrors2 = np.zeros(5,dtype='f')
    phierrors3 = np.zeros(5,dtype='f')
    deltax = np.zeros(5,dtype='f')
    deltay = np.zeros(5,dtype='f')
    
    #call function to get errors with respect to x
    # phierrors2[0],phierrors3[0],deltax[0],deltay[0]  = [0.00783108906154,0.0135762045364,125,62.5]
    phierrors2[0],phierrors3[0], deltax[0],deltay[0]  = [0.0730965364189,0.009865418,200,100]
    phierrors2[1],phierrors3[1], deltax[1],deltay[1]  = [0.0883140020345,0.012756212,250,125]
    phierrors2[2],phierrors3[2], deltax[2],deltay[2]  = [0.119114265008,0.02591136,500,250]
    phierrors2[3],phierrors3[3], deltax[3],deltay[3] = [0.2099130116382,0.114490655,1000,500]
    phierrors2[4],phierrors3[4], deltax[4],deltay[4] = [0.445442673596,0.269987787,2000,1000]    #plot l2 versus dx
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