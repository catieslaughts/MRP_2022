from kepler3 import *
from kepler3_addon import *

#setup pulled from test code in kepler3:

P     = 15.9 * u.year # years
eperi = 2002.3293 * u.year
a     = 0.1246 # should be in u.arcsec but Quantity barfs on this...
e     = 0.8831
i     = 134.87 * u.deg
anode = 226.53 * u.deg
w     = 64.98 * u.deg

t_orbit = np.linspace(1999,1999+P.value,300) * u.year
X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P,eperi,a,e,i,w,anode)

save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'data1.csv', overwrite=True)

#slightly different orbit for visualization:
P     = 15.9 * u.year # years
eperi = 2002.3293 * u.year
a     = 0.1246 # should be in u.arcsec but Quantity barfs on this...
e     = 0.95
i     = 134.87 * u.deg
anode = 226.53 * u.deg
w     = 64.98 * u.deg

t_orbit = np.linspace(1999,1999+P.value,300) * u.year
X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P,eperi,a,e,i,w,anode)

save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'data2.csv', overwrite=True)

P     = 15.9 * u.year # years
eperi = 2002.3293 * u.year
a     = 0.1246 # should be in u.arcsec but Quantity barfs on this...
e     = 0.8
i     = 134.87 * u.deg
anode = 226.53 * u.deg
w     = 64.98 * u.deg

t_orbit = np.linspace(1999,1999+P.value,300) * u.year
X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P,eperi,a,e,i,w,anode)

save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'data3.csv', overwrite=True)

#animate_2d('./data/', True)
animate_3d('./data/')
