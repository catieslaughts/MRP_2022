import numpy as np
import amuse.ext.orbital_elements as am
from amuse.units import units
import astropy.units as u
from kepler3_addon import *

@u.quantity_input
def kick_orbit_elements(t_orbit:u.year, stellar_mass:u.kg, part_masses:u.kg, P_i:u.year, a_i:u.AU, e_i, tperi_i:u.year, inc_i:u.rad, anode_i:u.rad, omega_i:u.rad, delta_v_arr:u.m/u.s, theta_arr:u.rad, phi_arr:u.rad):
	''' Takes in initial orbital elements, creates a binary, kicks the secondary, and spits 
	out the new orbital elements'''
	
	##Create initial value arrays
	
	a_arr = np.ones_like(part_masses).value*a_i
	e_arr = np.ones_like(part_masses).value*e_i
	inc_arr = np.ones_like(part_masses).value*inc_i
	anode_arr = np.ones_like(part_masses).value*anode_i
	omega_arr = np.ones_like(part_masses).value*omega_i
	
	##calculate true anomaly and save to array:
	E, f_i = tperi_to_Tanom(tperi_i, t_orbit, P_i, e_i)
	f_arr = np.ones_like(part_masses).value * f_i[0]
	
	##change units to AMUSE units:
	stellar_mass = stellar_mass.to(u.kg).value | units.kg
	part_masses = part_masses.to(u.kg).value | units.kg
	a_arr = a_arr.to(u.AU).value | units.AU
	f_arr = f_arr.to(u.rad).value | units.rad
	inc_arr = inc_arr.to(u.rad).value | units.rad
	anode_arr = anode_arr.to(u.rad).value | units.rad
	omega_arr = omega_arr.to(u.rad).value | units.rad
	
	##calculate relative positions and velocities
	pos, vels = am.rel_posvel_arrays_from_orbital_elements(stellar_mass, part_masses, a_arr, e_arr, f_arr, inc_arr, anode_arr, omega_arr, G=None)
	
	##calculate delta v in cartesian coordinates
	dv = covert_delta_v_cartesian(delta_v_arr, theta_arr, phi_arr)
	dv = dv.astype(float).reshape(1,3) | units.m/units.s
	
	##add velocities
	kicked_vels = np.zeros_like(vels)
	##TODO: is there an element-wise way to do this?
	for i in range(vels.shape[0]):
		kicked_vels[i] = (vels[i]+dv[i]).value_in(units.m/units.s)
	
	kicked_vels = kicked_vels | units.m/units.s
	
	#print(kicked_vels)
	##calculate total mass
	total_masses = part_masses + stellar_mass
	
	
	##calculate orbital elements for new velocities
	am.orbital_elements_for_rel_posvel_arrays(pos, vels, total_masses)
	
	return

@u.quantity_input	
def covert_delta_v_cartesian(delta_v_arr:u.m/u.s, theta_arr:u.rad, phi_arr:u.rad):
	
	cart_dv = [delta_v_arr*np.cos(phi_arr)*np.sin(theta_arr), delta_v_arr*np.sin(phi_arr)*np.sin(theta_arr), delta_v_arr*np.cos(theta_arr)]
	cart_dv = np.asarray(cart_dv)
	
	return cart_dv