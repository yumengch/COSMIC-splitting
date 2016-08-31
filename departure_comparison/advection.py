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
    x_edge,y_edge = np.linspace(xmin,xmax,nx+1), np.linspace(ymin,ymax,ny+1)
    x_cnter, y_cnter = x_edge[:-1] + 0.5*dx, y_edge[:-1] + 0.5*dy

    #------------------------------
    # initial conditions retrieval
    #------------------------------
    change = False
    t = -1
    cx, cy, u, v, X, Y, J  = solid(x_edge, y_edge, x_cnter, y_cnter, t, nt, dt, mesh, change)
    
    error = COSMIC(cx, cy, dx, dy, xmin, ymin, u, v, X, Y, dt, nt, J)
    print error


advection(initialProfile = solid, mesh ='quad', xmin = 0., xmax = 10000.,  ymin = 0., ymax = 10000., nx = 50, ny = 50, dt = 1., nt = 1)