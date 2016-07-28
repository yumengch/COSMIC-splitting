# #!/usr/bin/python

# Outer code for setting up the diffusion problem and calling the
# function to diffusion.
# from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# all the diffusion schemes and initial conditions
execfile("Schemes.py")
execfile("initialConditions.py")
execfile("diagnostics.py")

def advection(initialProfile=topHat, xmin = 0, xmax = 1, nx = 100, nt = 125, dt = 0.2):
    "Diffuse between x = xmin and x = xmax split over 40 spatial steps"
    "with diffusion coeffient K, time step dt for nt time steps"


    # spatial points for plotting and for defining initial conditions
    x = np.linspace(xmin, xmax, nx+1)
    def f(x):
        return -0.4
    def comput(x,L):
        #-------------------------------------
        # unequidistant computational domain 
        #-------------------------------------
        fx = f(x)
        # return np.where(x>0.5*L, fx+(x-0.5*L)*(L-fx)/(0.5*L),x*fx/(0.5*L))
        return np.where(x>=0, fx+x*(1-fx/xmax),fx+x*(1-fx/xmin))

    L = xmax - xmin

    #----------------------
    # computational domain 
    #----------------------
    x1 = comput(x, L)
    J = np.where(x>0.5*L, 2*(L-f(x))/L,2*f(x)/L)
    # print J
    y = pl.linspace(xmin, xmax, nx+1)
   
    # print x1
    # x1 = x

    dx = np.zeros([nx+1])
    # J = np.zeros_like(dx)
    dx[1:-1] = 0.5*(x1[2:] - x1[1:-1]) + 0.5*(x1[1:-1] - x1[:-2])
    dx[0] = 0.5*(x1[1] - x1[0]) + 0.5*(x1[-1] - x1[-2])
    dx[-1] = dx[0]

    # initial conditions
    u = np.zeros(nx+1)
    u[:] = 2.337								# velocity (u) change at this line!!!!

    print "Courant number:", np.max(u*dt/dx)
    # print np.max(J[1:]*u[1:] - J[:-1]*u[:-1]) 
    phiOld= initialProfile(x1)

    # analytical solution
    distanceTravelled = u*dt*nt
    # print distanceTravelled
    phiExact = initialProfile((x1 - distanceTravelled)%(xmax - xmin))

    # 1D PPM
    phi = COSMIC(phiOld.copy(), x1, u, dt, dx, nt,J)

    # # calculate the error norms 
    errors = errorNorms(phi.copy(), phiExact)
    print errors

    # # # plot the results in comparison to analytic
    plt.figure(1)
    plt.clf()
    plt.plot(x1,phi)
    plt.plot(x1,phiExact, 'r')

    plt.ylim(-0.1, 1.1)
    plt.show()


advection(initialProfile=cosBell, xmin = -1., xmax = 1., nx = 100, nt = 100,  dt = 5)