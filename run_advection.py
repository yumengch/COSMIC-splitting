# ---------------------------------------------------------------------------------
# Author: Yumeng Chen
# Run advection test casess
#---------------------------------------------------------------------------------
exec(open("advection.py").read())
#---------------------------------------------------------------------------------
# Solid body rotation on orthogonal grids
#---------------------------------------------------------------------------------
# small c
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 50, ny = 50, dt= 2., nt =300)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 1., nt =600)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 200, ny = 200, dt= 0.5, nt =1200)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 400, ny = 400, dt= 0.25, nt =2400)

# large c
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 50, ny = 50, dt= 20., nt =30)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 10, nt =60)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 200, ny = 200, dt= 5., nt =120)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 400, ny = 400, dt= 2.5, nt =240)

# same spatial resolution with varying dt
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 10., nt =60)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 5., nt =120)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 2., nt =300)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 1., nt =600)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 0.5, nt =1200)
advection(initialProfile = solid, mesh = 'orthog', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 0.25, nt =2400)

#---------------------------------------------------------------------------------
# Solid body rotation on V shape grids
#---------------------------------------------------------------------------------
# small c
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 50, ny = 50, dt= 2., nt =300)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 1., nt =600)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 200, ny = 200, dt= 0.5, nt =1200)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 400, ny = 400, dt= 0.25, nt =2400)

# large c
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 50, ny = 50, dt= 20., nt =30)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 10, nt =60)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 200, ny = 200, dt= 5., nt =120)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 400, ny = 400, dt= 2.5, nt =240)

# same spatial resolution with varying dt`
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 10., nt =60)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 5., nt =120)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 2., nt =300)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 1., nt =600)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 0.5, nt =1200)
advection(initialProfile = solid, mesh = 'V', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 0.25, nt =2400)

#---------------------------------------------------------------------------------
# Solid body rotation on quadratic grids
# ---------------------------------------------------------------------------------
# small c
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 50, ny = 50, dt= 2., nt =300)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 1., nt =600)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 200, ny = 200, dt= 0.5, nt =1200)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 400, ny = 400, dt= 0.25, nt =2400)

# large c
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 50, ny = 50, dt= 20., nt = 30)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 10, nt =60)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 200, ny = 200, dt= 5., nt =120)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 400, ny = 400, dt= 2.5, nt =240)

# same spatial resolution with varying dt
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 10., nt =60)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 5., nt =120)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 2., nt =300)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 1., nt =600)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 0.5, nt =1200)
advection(initialProfile = solid, mesh = 'quad', xmin = 0., xmax = 10000., ymin = 0., ymax = 10000., nx = 100, ny = 100, dt= 0.25, nt =2400)

# ---------------------------------------------------------------------------------
# Horizontal advection over orography h0 = 3km
# ---------------------------------------------------------------------------------
# small c
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =150, ny = 25, dt= 50, nt = 200)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 25, nt =400)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =600, ny = 100, dt= 12.5, nt =800)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =1200, ny = 200, dt= 6.25, nt =1600)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =1500, ny = 250, dt= 5, nt =2000)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =2400, ny = 400, dt= 3.125, nt =3200)

# small c
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =150, ny = 25, dt= 50, nt = 20)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 25, nt =40)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =600, ny = 100, dt= 12.5, nt =80)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =1200, ny = 200, dt= 6.25, nt =160)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =1500, ny = 250, dt= 5, nt =200)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =2400, ny = 400, dt= 3.125, nt =320)

# same spatial resolution with varying dt
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 50, nt =200)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 100, nt =100)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 200, nt =50)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 500, nt =20)
advection(initialProfile = orography, mesh = 'low', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 1000, nt =10)

#---------------------------------------------------------------------------------
# Horizontal advection over orography h0 = 6km
#---------------------------------------------------------------------------------
# small c
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =150, ny = 25, dt= 50, nt =200)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 25, nt =400)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =600, ny = 100, dt= 12.5, nt =800)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =1200, ny = 200, dt= 6.25, nt =1600)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =1500, ny = 250, dt= 5, nt =2000)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =2400, ny = 400, dt= 3.125, nt =3200)

# large c
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =150, ny = 25, dt= 50, nt =20)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 25, nt =40)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =600, ny = 100, dt= 12.5, nt =80)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =1200, ny = 200, dt= 6.25, nt =160)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =1500, ny = 250, dt= 5, nt =200)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =2400, ny = 400, dt= 3.125, nt =320)

# same spatial resolution with varying dt
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 50, nt =200)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 100, nt =100)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 200, nt =50)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 500, nt =20)
advection(initialProfile = orography, mesh = 'high', xmin = -150e3, xmax = 150e3, ymin = 0., ymax = 25e3, nx =300, ny = 50, dt= 1000, nt =10)

#---------------------------------------------------------------------------------
# deformational flow on orthogonal mesh
#---------------------------------------------------------------------------------
# small c
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 60, ny = 30, dt = 0.02, nt = 250)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.01, nt = 500)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 240, ny = 120, dt = 0.005, nt = 1000)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 360, ny = 180, dt = 0.01/3, nt = 1500)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 480, ny = 240, dt = 0.0025, nt = 2000)

# large c
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 60, ny = 30, dt = 0.2, nt = 25)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.1, nt = 50)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 240, ny = 120, dt = 0.05, nt = 1)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 360, ny = 180, dt = 0.1/3, nt = 150)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 480, ny = 240, dt = 0.025, nt = 200)

# same spatial resolution with varying dt
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.1, nt = 50)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.05, nt = 100)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.02, nt = 250)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.01, nt = 500)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.005, nt = 1000)
advection(initialProfile = deform, mesh = 'orthog', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.0025, nt = 2000)

# ---------------------------------------------------------------------------------
# deformational flow on W shape mesh
# ---------------------------------------------------------------------------------
# small c
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 60, ny = 30, dt = 0.02, nt = 250)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.01, nt = 500)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 240, ny = 120, dt = 0.005, nt = 1000)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 360, ny = 180, dt = 0.01/3, nt = 1500)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 480, ny = 240, dt = 0.0025, nt = 2000)

# large c
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 60, ny = 30, dt = 0.2, nt = 25)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.1, nt = 50)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 240, ny = 120, dt = 0.05, nt = 100)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 360, ny = 180, dt = 0.1/3, nt = 150)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 480, ny = 240, dt = 0.025, nt = 200)

# same spatial resolution with varying dt
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.1, nt = 50)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.05, nt = 100)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.02, nt = 250)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.01, nt = 500)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.005, nt = 1000)
advection(initialProfile = deform, mesh = 'W', xmin = 0., xmax = 2.*np.pi, ymin = -0.5*np.pi, ymax = 0.5*np.pi, nx = 120, ny = 60, dt = 0.0025, nt = 2000)