from collision_funcs import *
import numpy as np
from kepler3_addon import *
from kepler3 import *

#setup pulled from test code in kepler3:

P     = 1 * u.year # years
tperi = 0.7 * u.year
a     = 1 *u.AU#linear distance
e     = 0.8
inc   = 23.4406 * u.deg #23.4406
omega = 64.98 * u.deg #64.98
anode = 15 * u.deg #15

num_steps = 300

t_max = (P * 1.5)
t_orbit = np.linspace(0,0+t_max.value,num_steps) * u.year

E, f = tperi_to_Tanom(tperi, t_orbit, P, e)*u.rad
# print()

delta_v = .00000001*u.km/u.s
theta_i = 2*np.pi/4*u.rad
phi_i =  0*np.pi*u.rad

M = c.M_sun
m = c.M_earth

P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime = kick_kep_elements(delta_v, theta_i, phi_i, P, f[0], a, e, inc, omega, anode, M, m)

#convert f at collision to a new tperi:

tperi_prime = Tanom_to_tperi(f_prime, e_prime, P_prime)


print('theta: '+str(theta_i.value/np.pi)+ 'pi')
print()
print('inc_i: '+str(inc.to_value(u.rad)))
print('inc_prime: '+str(inc_prime.to_value(u.rad)))
print()
print('anode_i: '+str(anode.to_value(u.rad)))
print('anode_prime: '+str(anode_prime.to_value(u.rad)))
print()
print('omega_i: '+str(omega.to_value(u.rad)))
print('omega_prime: '+str(omega_prime.to_value(u.rad)))
# print()
# print('e_i: '+str(e))
# print('e_prime: '+str(e_prime))
# print('anode_prime: '+str(anode_prime.to_value(u.deg)))
print()

#original particle
X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P,tperi,a,e,inc,omega,anode)
save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'earthsun_prekick.csv', overwrite=True)

#kicked particle
#X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d_trueanom(t_orbit,P_prime,f_prime_array,a_prime,e_prime,inc_prime, omega_prime, anode_prime)
X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P_prime,tperi_prime,a_prime,e_prime,inc_prime, omega_prime, anode_prime)
save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, filename = 'earthsun_postkick.csv', overwrite=True)

# 
animate_3d('./earthsuntest/')