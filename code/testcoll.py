from collision_funcs import *
import numpy as np
from kepler3_helperfunctions import *
from kepler3 import *

#setup pulled from test code in kepler3:

P     = 1 * u.year # years
tperi = 1999.35 * u.year
a     = 1 *u.AU#arcsec
e     = 0.0167133
inc   = 23.4406  * u.deg
omega = 64.98 * u.deg
anode = .1 * u.deg

delta_v = 10000*u.m/u.s
theta_i = np.pi/3.6
phi_i = np.pi/2

M = c.M_sun
m = c.M_earth

t_init = 1999 * u.year

P_prime, tperi_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime = kick_kep_elements(delta_v, theta_i, phi_i, t_init, P, tperi, a, e, inc, omega, anode, M, m)

print('P_prime: '+str(P_prime))
print('tperi_prime: '+str(tperi_prime))
print('a_prime: '+str(a_prime))
print('e_prime: '+str(e_prime))
print('inc_prime: '+str(inc_prime.to_value(u.deg)))
print('omega_prime: '+str(omega_prime.to_value(u.deg)))
print('anode_prime: '+str(anode_prime.to_value(u.deg)))

#Animate:

t_orbit = np.linspace(1999, 1999+P_prime.value, 300) * u.year


X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P,tperi,a,e,inc,omega,anode)
save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'earthsun_prekick.csv', overwrite=True)


X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P_prime,tperi_prime,a_prime,e_prime,inc_prime,omega_prime,anode_prime)
save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'earthsun_postkick.csv', overwrite=True)

animate_3d('./earthsuntest/')