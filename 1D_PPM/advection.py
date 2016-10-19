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
    L = xmax - xmin
    dx = L/nx
    x_edge = np.linspace(xmin, xmax, nx+1)
    x_cntr = x_edge[:-1] + 0.5*dx

    # initial conditions
    u = np.zeros_like(x_cntr)
    u[:] = 38.224
    c = u*dt/dx
    print "Courant number:", np.max(u*dt/dx)
    # print np.max(J[1:]*u[1:] - J[:-1]*u[:-1]) 
    phiOld= initialProfile(x_cntr)

    # # analytical solution
    distanceTravelled = u*dt*nt
    # # print distanceTravelled
    phiExact = initialProfile((x_cntr - distanceTravelled)%(xmax - xmin))

    # 1D PPM
    phi = COSMIC(phiOld, c, u, x_cntr, x_edge[:-1], dt, nt)

    # # calculate the error norms 
    # errors = errorNorms(phi.copy(), phiExact)
    # print errors

    # # # # plot the results in comparison to analytic
    plt.figure(1)
    plt.clf()
    plt.plot(x_cntr,phi)
    plt.plot(x_cntr,phiExact, 'r')

    # plt.ylim(-0.1, 1.1)
    plt.show()


advection(initialProfile=cosBell, xmin = 0., xmax = 1., nx = 50, nt = 5000,  dt = 0.0025)