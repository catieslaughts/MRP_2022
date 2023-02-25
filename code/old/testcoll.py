from collision_funcs import *
import numpy as np
from kepler3_addon import *
from kepler3 import *

threeD = False
#setup pulled from test code in kepler3:

P     = 1. * u.year # years
tperi = 0. * u.year
a     = 1. *u.AU#linear distance
e     = 0.
inc   = 0.001 *np.pi * u.rad
omega = 0. *np.pi * u.rad #
anode = 0. *np.pi * u.rad

num_steps = 600

t_max = (P * 1)
t_orbit = np.linspace(0,0+t_max.value,num_steps) * u.year

delta_v = .01*u.m/u.s
theta_i = (0.5)*np.pi*u.rad
phi_i =  (0.)*np.pi*u.rad

M = c.M_sun
m = c.M_earth

E, f = tperi_to_Tanom(tperi, t_orbit, P, e) #units of radians

#P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime = kick_kep_elements(delta_v, theta_i, phi_i, P, f[0], a, e, inc, omega, anode, M, m)
if threeD:
	P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime = kick_kep_elements(delta_v, theta_i, phi_i, P, f[0], a, e, inc, omega, anode, M, m)
else:
	P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime = kick_kep_elements_2d(delta_v, theta_i, phi_i, P, f[0], a, e, omega, anode, M, m)

#convert f after collision to a new tperi:
tperi_prime = Tanom_to_tperi(f_prime, e_prime, P_prime)


# print('theta: '+str(theta_i.value/np.pi)+ 'pi')
# print()
print('e_i: '+str(e))
print('e_prime: '+str(e_prime))
print()
print('inc_i: '+str(inc.to_value(u.rad)))
print('inc_prime: '+str(inc_prime.to_value(u.rad)))
print()
print('anode_i: '+str(anode.to_value(u.rad)))
print('anode_prime: '+str(anode_prime.to_value(u.rad)))
print()
print('omega_i: '+str(omega.to_value(u.rad)))
print('omega_prime: '+str(omega_prime.to_value(u.rad)))
print()
print('f_i: '+str(f[0].to_value(u.rad)))
print('f_prime: '+str(f_prime.to_value(u.rad)))
# print()
# print('e_i: '+str(e))
# print('e_prime: '+str(e_prime))
# print('anode_prime: '+str(anode_prime.to_value(u.deg)))
print()

#original particle
X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P,tperi,a,e,inc,omega,anode)
save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'earthsun_prekick.csv', overwrite=True)

X_0 = Xs.copy()
Y_0 = Ys.copy()
Z_0 = Zs.copy()

#kicked particle
X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P_prime,tperi_prime,a_prime,e_prime,inc_prime, omega_prime, anode_prime)
save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'earthsun_postkick.csv', overwrite=True)

dist = np.sqrt((X_0[0] - Xs[0])**2 + (Y_0[0] - Ys[0])**2 + (Z_0[0] - Zs[0])**2)

print(dist)


# 
animate_3d('./earthsuntest/')