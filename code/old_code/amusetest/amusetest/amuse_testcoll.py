import numpy as np
import astropy.units as u
import astropy.constants as c
from kick_amuse import *

P     = 1. *u.yr # years
tperi = 0. *u.yr
a     = 1. *u.AU #linear distance
e     = 0.5
inc   = 0. *np.pi *u.rad
omega = 0. *np.pi *u.rad
anode = 0. *np.pi *u.rad

num_particles = 1

num_steps = 600

t_max = (P * 1)
t_orbit = np.linspace(0,0+t_max.value,num_steps) *u.yr

delta_v_arr = np.ones(num_particles)*(.01)*u.m/u.s
theta_arr = np.ones(num_particles)*(0.5)*np.pi*u.rad
phi_arr =  np.ones(num_particles)*(0.)*np.pi*u.rad

stellar_mass = c.M_sun
m = c.M_earth
part_masses = np.ones(num_particles)*(m/1000)

kick_orbit_elements(t_orbit, stellar_mass, part_masses, P, a, e, tperi, inc, anode, omega, delta_v_arr, theta_arr, phi_arr)

#original particle
# X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P,tperi,a,e,inc,omega,anode)
# save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'earthsun_prekick.csv', overwrite=True)
# 
# # X_0 = Xs.copy()
# # Y_0 = Ys.copy()
# # Z_0 = Zs.copy()
# 
# #kicked particle
# X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P_prime,tperi_prime,a_prime,e_prime,inc_prime, omega_prime, anode_prime)
# save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'earthsun_postkick.csv', overwrite=True)
# 
# # dist = np.sqrt((X_0[0] - Xs[0])**2 + (Y_0[0] - Ys[0])**2 + (Z_0[0] - Zs[0])**2)
# # 
# # print(dist)
# 
# 
# # 
# animate_3d('./earthsuntest/')