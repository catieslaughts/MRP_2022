from collision_funcs import *
import numpy as np
from kepler3_addon import *
from kepler3 import *

from tqdm import tqdm

## This code is to attempt to systematically evaluate where kick_kep_elements breaks down for debugging purposes

threeD = False

#input parameters to cycle over:

tperi_arr = np.asarray([0.001,.25,.5,.75]) * u.year #[0.,.25,.5,.75]
e_arr     = np.asarray([0.001,.25,.5,.75])

if threeD:
	inc_arr   = np.asarray([0.001,0.5,1,1.5]) *np.pi * u.rad
else:
	inc_arr   = np.asarray([0.001])*np.pi * u.rad
	
omega_arr = np.asarray([0.001,0.5,1,1.5]) *np.pi * u.rad # 
anode_arr = np.asarray([0.001,0.5,1,1.5]) *np.pi * u.rad #

delta_v_arr = [0.001,5]*u.km/u.s
theta_arr = np.asarray([0.001,0.5,1,1.5])*np.pi*u.rad
phi_arr =  np.asarray([0.001,0.5,1,1.5])*np.pi*u.rad


#constants:
M = c.M_sun
m = c.M_earth

P_i = 1. * u.year
a_i = 1. * u.AU

#Simulation parameters:
num_steps = 300
t_max = (P_i * 1)
t_orbit = np.linspace(0,0+t_max.value,num_steps) * u.year

tperi_grid, e_grid, inc_grid, omega_grid, anode_grid, delta_v_grid, theta_grid, phi_grid = np.meshgrid(tperi_arr, e_arr, inc_arr, omega_arr, anode_arr, delta_v_arr, theta_arr, phi_arr)

tperi_flat, e_flat, inc_flat, omega_flat, anode_flat, delta_v_flat, theta_flat, phi_flat = tperi_grid.flatten(), e_grid.flatten(), inc_grid.flatten(), omega_grid.flatten(), anode_grid.flatten(), delta_v_grid.flatten(), theta_grid.flatten(), phi_grid.flatten()

#because file i/o is so time consuming, I'll run everything first, then save it

dist_flat = np.zeros_like(tperi_flat.value) * u.AU

print('Running...')
for idx in tqdm(range(tperi_flat.size)):
	#print(idx)
	
	#pull out values from grids
	tperi_i, delta_v, theta, phi, e_i, inc_i, omega_i, anode_i = tperi_flat[idx], delta_v_flat[idx], theta_flat[idx], phi_flat[idx], e_flat[idx], inc_flat[idx], omega_flat[idx], anode_flat[idx]
	
	#calculate f
	E, f = tperi_to_Tanom(tperi_flat[idx], t_orbit, P_i, e_i)
	
	#find kicked elements
	if threeD:
		P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime = kick_kep_elements(delta_v, theta, phi, P_i, f[0], a_i, e_i, inc_i, omega_i, anode_i, M, m)
	else:
		P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime = kick_kep_elements_2d(delta_v, theta, phi, P_i, f[0], a_i, e_i, omega_i, anode_i, M, m)
	
	#calculate tperi
	tperi_prime = Tanom_to_tperi(f_prime, e_prime, P_prime)
	
	#simulate original particle
	X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit, P_i, tperi_i, a_i, e_i, inc_i, omega_i, anode_i)
	
	X_0 = Xs.copy()
	Y_0 = Ys.copy()
	Z_0 = Zs.copy()
	
	#simulate kicked particle
	X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P_prime,tperi_prime,a_prime,e_prime,inc_prime, omega_prime, anode_prime)
	
	dist = np.sqrt((X_0[0] - Xs[0])**2 + (Y_0[0] - Ys[0])**2 + (Z_0[0] - Zs[0])**2)
	dist_flat[idx] = dist
	
	
worked = dist_flat < 1E-5*u.AU
worked_out = np.array(worked, dtype=str)
#print(worked_out)

#to save to file:
header = 'tperi,e,inc,omega,anode,theta,phi,delta_v,distance,worked'

if threeD:
	filename = 'gridtest_output.txt'
else:
	filename = 'gridtest_output_2d.txt'

np.savetxt(filename, np.c_[tperi_flat.to(u.yr).value, e_flat, inc_flat.to(u.rad).value/np.pi, omega_flat.to(u.rad).value/np.pi, anode_flat.to(u.rad).value/np.pi, theta_flat.to(u.rad).value/np.pi, phi_flat.to(u.rad).value/np.pi, delta_v_flat.to(u.km/u.s).value, dist_flat.to(u.AU).value, worked_out], header = header, delimiter=',', fmt='%s')










