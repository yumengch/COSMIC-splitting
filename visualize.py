from math import radians, sin, cos, sqrt, asin
import numpy as np
import matplotlib.pyplot as plt
def initial(x_edge, y_edge, x_cntr, y_cntr, t, nt, dt, mesh, change):
  dlat = np.zeros([nlat+1])
  dlat[:] = np.pi/nlat
  dlon = 2*np.pi/nlon
  slat = np.zeros([nlat+1])
  # nlat = nlat +1
  # nlon = nlon +1
  r = np.zeros([nlon, nlat])
  R = 6.37122e6
  for i in xrange(nlon):
      for j in xrange(nlat):
        r[i,j] = haversine(0, 1.5*np.pi, (-0.5*np.pi+j*dlat[j] + 0.5*dlat[j+1]), i*dlon + 0.5*dlon)

  R0 = R/3
  q0 = 1000.
  q = np.zeros([nlon, nlat])
  for i in xrange(nlon):
      for j in xrange(nlat):
      
          if r[i,j]<R0:
              q[i, j] = 0.5*q0*(1+np.cos(np.pi*r[i,j]/R0))
          else:
              q[i, j] = 0.

  beta = 0.0
  u0 = 2*np.pi*R/(12*3600*24)
  psi = np.zeros([nlon, nlat+1])
  slat[0] = -0.5*np.pi
  for j in xrange(1,nlat+1):
      slat[j] = slat[j-1] + dlat[j-1]
  print np.round(slat*180/np.pi,3)
  for j in xrange(nlat+1):
      for i in xrange(nlon):

          psi[i, j] = -R*u0*(np.sin(slat[j])*np.cos(beta) - np.cos(i*dlon)*np.cos(slat[j])*np.sin(beta))

  u = np.zeros([nlon, nlat])
  v = np.zeros([nlon, nlat])
  for j in xrange(nlat):
      for i in xrange(nlon-1):
          u[i, j] = -(psi[i,j+1] - psi[i,j])/(R*dlat[j+1])
          v[i, j] = (psi[i+1,j] - psi[i,j])/(R*np.cos((-0.5*np.pi+j*dlat[j+1]))*dlon)

  v[nlon-1,:] = v[0, :]
  u[nlon-1,:] = u[0, :]

  return u,v,q


 
def haversine(lat1, lon1, lat2, lon2):
 
  R = 6.37122e6 # Earth radius in kilometers
 
  dLat = lat2 - lat1
  dLon = lon2 - lon1
  lat1 = lat1
  lat2 = lat2
 
  a = sin(dLat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
  c = 2*asin(sqrt(a))
 
  return R * c

u, v, q = initial(nlon = 10, nlat = 10)
print q[:,0]
print np.max(u), np.max(v)
lons = np.linspace(0,2*np.pi, 10)
lats = np.linspace(-0.5*np.pi, 0.5*np.pi, 10)
plt.figure(1)
plt.clf()
plt.contourf(lons, lats, np.transpose(q))
# CS = plt.contour(lons, lats, np.transpose(u))
# plt.clabel(CS, inline=1, fontsize=10)
# plt.quiver(lons, lats, np.transpose(u), np.transpose(v))
# plt.colorbar()
plt.show()