from kepler3_addon import *
from fileio import save_kep3d, save_kep3d_short, empty_data_directory
from spheres import fibonacci_sphere, draw_spherical
from kepler3_catherine import kep3d
import astropy.units as u
import astropy.constants as c
import os
from tqdm import tqdm

@u.quantity_input
def simulate_spherical(n_shells = 3, kick_vel:u.m/u.s = 1.*u.km/u.s, P_i:u.year = 165.*u.year, 
		tperi_i:u.year = 2000 * u.year, e_i = 0., omega_i:u.rad = 0*u.rad, anode_i:u.rad = 0*u.rad, 
		M:u.kg = 1.0*u.Msun, m:u.kg = 1.0*u.Mearth, num_per_shell = 400, t_start:u.year = 2000*u.year, 
		delta_t:u.year = 50*u.year, n_steps = 500, show_pre = False, save_dir = './data/', 
		overwrite = False):
	'''creates a simulation for a spherically symmetric explosion
		
		inputs:
			rel_mass:    (numpy array) The relative masses of each spherical shell. The number of shells
					     is determined by the length of the array
			vel: 	     (numpy array or float) Determines the kick velocities for each shell. 
					     If an array, must be of same length as rel_mass, if an int, taken to
					     be the max velocity, and a linspace function is used to determine the
					     others
			P_i:	     Period of system pre-explosion
			f_i:	     true anomaly of system pre-explosion
			omega_i:     argument of periapsis of system pre-explosion
			anode_i:     longitude of ascending node of system pre-explosion
			M:		     primary (stellar) mass
			m:		     secondary (planitesimal) mass pre-explosion
			n_particles: the number of particles to simulate
			t_start:	 the time of the explosion
			delta_t:	 how much time to simulate after the explosion
			show_pre:	 Boolean or time, indicates whether you want to simulate the 
						 progenitor for some amount of time before the explosion, if a
						 number, that much time pre-explosion, if "True", defaults to 
						 five years
			save_dir:    directory to save output files to
			overwrite:   if it's ok for the program to overwrite the current directory at save_dir
			
	'''
	
	#calculate step rate (steps/year)
	step_rate = n_steps/delta_t
	
	#initial i is 0:
	inc_i = 0 *u.rad
	
	# given the masses and period, we really have to calculate a
	a_i = np.power(c.G*(M+m)*P_i*P_i / (4*np.pi*np.pi),1/3.).to(u.au)
	
	#calculate orbital velocity
	vorb=np.sqrt(c.G * (M+m)/ a_i).to(u.km/u.s)
	print('Initial orbital velocity: '+str(vorb))
	
	#create array of shell numbers:
	shellnums = np.zeros(n_shells*num_per_shell)
	for num in range(n_shells):
		shellnums[num*num_per_shell:(num+1)*num_per_shell] = num
	
	#print(shellnums)
	
	#prepare to save data:
	if os.path.exists(save_dir):
		if len(os.listdir(save_dir)) != 0: #if the directory is not empty
			if overwrite:
				empty_data_directory(save_dir)
			else:
				print('Save directory at '+save_dir+' is not empty. Select a new directory or set overwrite = True')
				return
	else:
		print('Save directory does not exist, creating new directory at '+save_dir)
		os.mkdir(save_dir)
		
	###PROGENITOR:
	#create single orbit:
	t_pro_orbit = np.linspace(t_start,t_start+P_i,(step_rate*P_i).astype(int))
	_, _, Xa, Ya, Za, _, _, _ = kep3d(t_pro_orbit, P_i, tperi_i, a_i, e_i, inc_i, omega_i, anode_i)
	
	save_kep3d_short(Xa, Ya, Za, filename = 'orbit.csv', overwrite=overwrite, filepath = save_dir)
	
	#to visualize progenitor:
	prog_steps = 0
	if isinstance(show_pre, bool):
		if show_pre:
			prog_steps = int(step_rate*30*u.year)
			prog_start = t_start-30*u.year
			t_prog = np.linspace(prog_start, t_start, prog_steps)
			
			_, _, Xs, Ys, Zs, _, _, _ = kep3d(t_prog, P_i, tperi_i, a_i, e_i, inc_i, omega_i, anode_i)
			save_kep3d_short(Xs, Ys, Zs, filename = 'progenitor.csv', overwrite=overwrite, post_steps = n_steps, filepath = save_dir)
			
	else:
		try:
			show_pre.to(u.year)
		except:
			print('show_pre must be a boolean or quantity with units of time')
			return
		
		prog_steps = step_rate*show_pre
		prog_start = t_start+show_pre
		t_prog = np.linspace(prog_start, t_start, prog_steps)
		
		_, _, Xs, Ys, Zs, _, _, _ = kep3d(t_prog, P_i, tperi_i, a_i, e_i, inc_i, omega_i, anode_i)
		save_kep3d_short(Xs, Ys, Zs, filename = 'progenitor.csv', overwrite=overwrite, post_steps = num_steps, filepath = save_dir)
	
	## Calculate shell velocities:
	if isinstance(kick_vel.value,np.ndarray): #if vel is passed in as an array
		if kick_vel.size != n_shells:
			print('error, kick_vel must be array of same size as rel_mass OR float')
			return
		else:
			vel_arr = kick_vel.copy()
	else:
		vel_arr = np.linspace(0, kick_vel, n_shells+1)[1:]
		
	##Assuming each particle is equal mass, calc number of particles per shell:
	#normalize massses
# 	norm_mass = rel_mass/(rel_mass.sum())
# 	num_per_shell = (norm_mass*n_particles).astype(int)
	
#	num_per_shell = (n_particles/n_shells)
		
	##Create shells:
	r, theta, phi = fibonacci_sphere(n_shells, num_per_shell, vel_arr)	
	
	#calculate f_i# calculate mean anomaly for t_explode
	M_explode = u.rad * 2*np.pi*(t_start-tperi_i)/P_i
	# calculate the true anomaly for this mean anomaly
	(E_explode, f_i) = kepler_solve(e_i, M_explode)
	
	#calculate new orbital elements
	P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime = kick_kep_elements(r, theta, phi, P_i, f_i, a_i, e_i, omega_i, anode_i, M, m)
	
	#calculate orbits, save data:
	t_post = np.linspace(t_start, t_start+delta_t, n_steps)
	
	# now we want to get the periastron times for all the exploded particle orbits...
	# calculate eccentric anomaly from true anomaly
	E_prime = np.arctan2(np.sqrt(1 - e_prime**2) * np.sin(f_prime), e_prime + np.cos(f_prime))
	
	# calculate the mean anomaly fom the eccentric anomaly
	M_prime = E_prime - e_prime*np.sin(E_prime)*u.rad
	
	tperi_prime = t_start - P_prime*M_prime/(2*np.pi*u.rad)
	
	for idx in tqdm(range(P_prime.size)):
		shellnum = shellnums[idx]
		X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_post,P_prime[idx],tperi_prime[idx],a_prime[idx],e_prime[idx],inc_prime[idx],omega_prime[idx],anode_prime[idx])
		
		shellarr = np.ones_like(Xs)*shellnum
		
		filename = 'cloud'+str(idx)+'.csv'
		save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, shellarr, filename = filename, overwrite=True, pre_steps = prog_steps, filepath = save_dir)
 		
 		
	return 

simulate_spherical(n_shells = 4, num_per_shell = 400, tperi_i = 2025*u.year, show_pre = False, overwrite = True, save_dir = './test_subset/')


	