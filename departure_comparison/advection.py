import numpy as np
import matplotlib.pyplot as plt
#---------------------------------------------------------------------------------
# Author: Yumeng Chen
# The main program to call the COSMIC splitting with different test cases
#---------------------------------------------------------------------------------

exec(open("initialConditions.py").read())
exec(open("scheme1.py").read())
def advection(initialProfile, mesh, xmin, xmax,  ymin, ymax, nx, ny, dt, nt):
    np.seterr(all='raise')

    #-----------------------
    # Basic grid information
    #-----------------------
    Lx, Ly= xmax-xmin,ymax-ymin
    dx,dy = Lx/(nx-1), Ly/ny
    x,y = np.linspace(xmin,xmax,nx), np.linspace(ymin,ymax,ny+1)
    print x[1]- x[0], y[1]- y[0]

    #------------------------------
    # initial conditions retrieval
    #------------------------------
    change = False
    t = -1
    phiOld, cx, cy, u, v, X, Y, J  = initialProfile(x, y, xmin, ymin, nx, ny, Lx, Ly,t, nt, dt, mesh, change)
    
    error = COSMIC(phiOld, cx, cy, dx, dy, xmin, ymin, u, v, X, Y, dt, nt, J, initialProfile, mesh, change)
    print error


advection(initialProfile = solid, mesh ='quad', xmin = 0., xmax = 10000.,  ymin = 0., ymax = 10000., nx = 100, ny = 100, dt = 2, nt = 1)