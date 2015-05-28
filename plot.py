import numpy as np
import matplotlib.pyplot as plt

def hist(X,Y,phiCOS):
	fig = plt.figure(1)
	ax = fig.gca(projection = '3d')
	c = ax.plot_surface(X,Y,phiCOS,rstride=1, cstride=1, cmap=cm.coolwarm,
       linewidth=0, antialiased=False)
	fig.colorbar(c)
	plt.show()

def contours(X,Y,phiCOS):
	plt.figure(1)
	plt.contour(X,Y,phiCOS)
	plt.colorbar()
	plt.show()

def init(X,Y,u,v,phiCOS):
	plt.figure(1)
	plt.streamplot(X,Y,u,v)
	plt.contour(X,Y,phiCOS)
	plt.show()
