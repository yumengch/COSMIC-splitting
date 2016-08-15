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
    dx,dy = Lx/nx, Ly/ny
    x,y = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
    print y[50]

    #------------------------------
    # initial conditions retrieval
    #------------------------------
    change = False
    t = -1
    cx, cy, u, v, X, Y, J  = initialProfile(x, y, xmin, ymin, nx, ny, Lx, Ly,t, nt, dt, mesh, change)
    
    error = COSMIC(cx, cy, dx, dy, xmin, ymin, u, v, X, Y, dt, nt, J)
    print error


advection(initialProfile = solid, mesh ='quad', xmin = 0., xmax = 10000.,  ymin = 0., ymax = 10000., nx = 50, ny = 50, dt = 0.1, nt = 1)