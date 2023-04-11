from kepler3_addon import *
from fileio import *
from spheres import fibonacci_sphere, draw_spherical
from kepler3_catherine import kep3d
import astropy.units as u
import astropy.constants as c
import os
from tqdm import tqdm
from huge_sphere import *
	

#use this one to make cute animations with the progenitor
@u.quantity_input
def simulate_spherical(n_shells = 3, kick_vel:u.m/u.s = 1.*u.km/u.s, v_i:u.km/u.s = 6.*u.km/u.s, 
		t_firstpass:u.day = 750 * u.day, e_i = 0., omega_i:u.rad = 0*u.rad, anode_i:u.rad = 0*u.rad, 
		M:u.kg = 1.0*u.Msun, m:u.kg = 1.0*u.Mearth, num_per_shell = 400, t_start:u.year = 2000*u.year, 
		delta_t:u.year = 50*u.year, n_steps = 100, R_star:u.au = 1*u.Rsun, show_pre = False, save_dir = './data/', 
		overwrite = False, paramfile = 'params.csv'):
		
	'''creates a simulation for a spherically symmetric explosion, with initial orbital velocity input instead of period
		
		inputs:
			rel_mass:    (numpy array) The relative masses of each spherical shell. The number of shells
					     is determined by the length of the array
			vel: 	     (numpy array or float) Determines the kick velocities for each shell. 
					     If an array, must be of same length as rel_mass, if an int, taken to
					     be the max velocity, and a linspace function is used to determine the
					     others
			v_i:	     Orbital velocity at explosion?
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
	
	tperi_i = t_start + (kick_vel*0.8/(u.km/u.s))*t_firstpass
	
	#calculate step rate (steps/year)
	step_rate = n_steps/delta_t
	
	#initial i is 0:
	inc_i = 0 *u.rad
	
	# given the masses and orbital velocity, we really have to calculate a
	a_i = (c.G*(M+m)/(v_i**2)).to(u.au)
	
	#calculate Period
	P_i=np.sqrt((a_i**3)*4*np.pi*np.pi/(c.G*(M+m))).to(u.yr)
	print('Initial Period: '+str(P_i))
	
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
	shellnum = np.ones_like(Xa)*1000
	
	save_kep3d_npz(Xa, Ya, Za, shellnum, filename = 'orbit', overwrite=overwrite, filepath = save_dir)
	
	#to visualize progenitor:
	prog_steps = 0
	if isinstance(show_pre, bool):
		if show_pre:
			prog_steps = int(step_rate*1*u.year)
			prog_start = t_start-1*u.year
			t_prog = np.linspace(prog_start, t_start, prog_steps)
			
			_, _, Xs, Ys, Zs, _, _, _ = kep3d(t_prog, P_i, tperi_i, a_i, e_i, inc_i, omega_i, anode_i)
			
			post_steps = np.full(n_steps, np.nan)
			Xs = np.concatenate((Xs,post_steps))
			Ys = np.concatenate((Ys,post_steps))
			Zs = np.concatenate((Zs, post_steps))
			
			
			shellnum = np.ones_like(Xs)*1000
			
			save_kep3d_npz(Xs, Ys, Zs, shellnum, filename = 'progenitor', overwrite=overwrite, filepath = save_dir)
			#save_kep3d_npz(Xs, Ys, Zs, shellnum, filename = 'progenitor.csv', overwrite=overwrite, post_steps = n_steps, filepath = save_dir)
			
	else:
		try:
			show_pre.to(u.year)
		except:
			print('show_pre must be a boolean or quantity with units of time')
			return
		
		prog_steps = int(step_rate*show_pre)
		prog_start = t_start-show_pre
		t_prog = np.linspace(prog_start, t_start, prog_steps)
		
		_, _, Xs, Ys, Zs, _, _, _ = kep3d(t_prog, P_i, tperi_i, a_i, e_i, inc_i, omega_i, anode_i)
			
		post_steps = np.full(n_steps, np.nan)
		Xs = np.concatenate((post_steps,Xs))
		Ys = np.concatenate((post_steps,Ys))
		Zs = np.concatenate((post_steps,Zs))
		
		
		shellnum = np.ones_like(Xs)*1000
		
		save_kep3d_npz(Xs, Ys, Zs, shellnum, filename = 'progenitor', overwrite=overwrite, filepath = save_dir)
	
	
	#save input params to file:
	write_param_file(kick_vel, v_i, t_firstpass, t_start, delta_t, e_i, omega_i, anode_i, M, m, num_per_shell, n_steps, prog_steps, n_shells, R_star, paramfile)
	
	## Calculate shell velocities:
	if isinstance(kick_vel.value,np.ndarray): #if vel is passed in as an array
		if kick_vel.size != n_shells:
			print('error, kick_vel must be array of same size as rel_mass OR float')
			return
		else:
			vel_arr = kick_vel.copy()
	else:
		vel_arr = np.linspace(0, kick_vel, n_shells+1)[1:]
	
		
	##Create shells:
	r, theta, phi = fibonacci_sphere(n_shells, num_per_shell, vel_arr)	
	
	
	# calculate mean anomaly for t_explode
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
	
	files = np.zeros(P_prime.size, dtype = str)
	
	for idx in tqdm(range(P_prime.size)):
		shellnum = shellnums[idx]
		currtheta = theta[idx]
		currphi = phi[idx]
		X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_post,P_prime[idx],tperi_prime[idx],a_prime[idx],e_prime[idx],inc_prime[idx],omega_prime[idx],anode_prime[idx])
		
		pre_steps = np.full(prog_steps, np.nan)
		Xs = np.concatenate((pre_steps,Xs))
		Ys = np.concatenate((pre_steps,Ys))
		Zs = np.concatenate((pre_steps,Zs))
		
		filename = 'cloud'+str(idx)
		files[idx] = filename
		save_kep3d_npz(Xs, Ys, Zs, shellnum, currtheta, currphi, filename = filename, overwrite=True, filepath = save_dir)
 		
 	
	return step_rate

#use this one to actually run stuff
@u.quantity_input
def simulate_spherical_docinput(paramfile = 'params.csv', save_dir = './data/', overwrite = False):
		
	'''creates a simulation for a spherically symmetric explosion, with initial orbital velocity input instead of period
		
		inputs:
			save_dir:    directory to save output files to
			overwrite:   if it's ok for the program to overwrite the current directory at save_dir
			
	'''
	
	#read in data
	
	kick_vel,v_i,t_firstpass,t_start,delta_t,e_i,omega_i,anode_i,M,m,num_per_shell,n_steps,n_pre_steps,n_shells,R_star = read_param_file(paramfile)
	
	
	tperi_i = t_start + (kick_vel*0.8/(u.km/u.s))*t_firstpass
	
	#calculate step rate (steps/year)
	step_rate = n_steps/delta_t
	
	#initial i is 0:
	inc_i = 0 *u.rad
	
	# given the masses and orbital velocity, we really have to calculate a
	a_i = (c.G*(M+m)/(v_i**2)).to(u.au)
	
	#calculate Period
	P_i=np.sqrt((a_i**3)*4*np.pi*np.pi/(c.G*(M+m))).to(u.yr)
	print('Initial Period: '+str(P_i))
	
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
		
	
	## Calculate shell velocities:
	if isinstance(kick_vel.value,np.ndarray): #if vel is passed in as an array
		if kick_vel.size != n_shells:
			print('error, kick_vel must be array of same size as rel_mass OR float')
			return
		else:
			vel_arr = kick_vel.copy()
	else:
		vel_arr = np.linspace(0, kick_vel, n_shells+1)[1:]
		
	##Create shells:
	r, theta, phi = fibonacci_sphere(n_shells, num_per_shell, vel_arr)	
	

	# calculate mean anomaly for t_explode
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
	
	files = np.zeros(P_prime.size, dtype = str)
	
	for idx in tqdm(range(P_prime.size)):
		shellnum = shellnums[idx]
		currtheta = theta[idx]
		currphi = phi[idx]
		_,_, Xs, Ys, Zs, _,_,_ = kep3d(t_post,P_prime[idx],tperi_prime[idx],a_prime[idx],e_prime[idx],inc_prime[idx],omega_prime[idx],anode_prime[idx])
		
		filename = 'cloud'+str(idx)
		files[idx] = filename
		save_kep3d_npz(Xs, Ys, Zs, shellnum, currtheta, currphi, filename = filename, overwrite=True, filepath = save_dir)
 		
 		
	return step_rate

def simulate_postprune(paramfile = 'params.csv', foi = '', extended_set = '', read_dir = './data/', save_dir = './extended_data/', overwrite = False):
		
	'''
			
	'''
	
	#read in data
	
	kick_vel,v_i,t_firstpass,t_start,delta_t,e_i,omega_i,anode_i,M,m,num_per_shell,n_steps,n_pre_steps,n_shells,R_star = read_param_file(paramfile)
	
	
	tperi_i = t_start + (kick_vel*0.8/(u.km/u.s))*t_firstpass
	
	#calculate step rate (steps/year)
	step_rate = n_steps/delta_t
	
	#initial i is 0:
	inc_i = 0 *u.rad
	
	# given the masses and orbital velocity, we really have to calculate a
	a_i = (c.G*(M+m)/(v_i**2)).to(u.au)
	
	#calculate Period
	P_i=np.sqrt((a_i**3)*4*np.pi*np.pi/(c.G*(M+m))).to(u.yr)
	print('Initial Period: '+str(P_i))
	
	#create array of shell numbers:
# 	shellnums = np.zeros(n_shells*num_per_shell)
# 	for num in range(n_shells):
# 		shellnums[num*num_per_shell:(num+1)*num_per_shell] = num
	
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
		
	
	## Calculate shell velocities:
	if isinstance(kick_vel.value,np.ndarray): #if vel is passed in as an array
		if kick_vel.size != n_shells:
			print('error, kick_vel must be array of same size as rel_mass OR float')
			return
		else:
			vel_arr = kick_vel.copy()
	else:
		vel_arr = np.linspace(0, kick_vel, n_shells+1)[1:]
		
	##TODO: read in foi particles, find new particle locations, recreate r, theta, phi:
	
	#new particles
	new_theta, new_phi = get_giant_sphere()
	
	#foi particles
	files = np.loadtxt('foi.csv', dtype = 'str', delimiter=',')
	if files.size == 1:
		files = np.asarray([files])
	
	og_theta = np.zeros(files.size) * u.rad
	og_phi = np.zeros(files.size) * u.rad
	og_shellnum = np.zeros(files.size, dtype = int)
	
	print('Reading in og data...')
	for idx in tqdm(range(files.size)):
		filename = files[idx]
		if not filename.startswith('.'):
			try:	
				_, *_, og_shellnum[idx], og_theta[idx], og_phi[idx] = read_kep3d_npz(read_dir, filename)
			except:
				print('Error: '+filename)
	
	#Calc radius:
	ang_rad = np.sqrt(4*np.pi / num_per_shell) *u.rad#4pi steradians per sphere, divided by the number of particles per shell
	search_rad = haversin(ang_rad) #ASK MATTHEW ABOUT THIS
	#print(ang_rad)
	
	#find nearby:
	print('Finding nearby particles..')
	p_bar = tqdm(range(og_theta.size))
	for curr_shell in np.unique(og_shellnum):
	
		select_arr = np.zeros(new_theta.size, dtype=int)
		theta = np.array([])*u.rad
		phi = np.array([])*u.rad
		shellnums = np.array([], dtype = int)
		vels = np.array([])*u.km/u.s
		
		rel_theta = og_theta[og_shellnum == curr_shell] #relevant list of thetas
		rel_phi = og_phi[og_shellnum == curr_shell]
		
		for part_idx, curr_theta in enumerate(rel_theta):
			curr_phi = rel_phi[part_idx]
		
			dist_arr=haversin(new_theta-curr_theta)+np.cos(new_theta)*np.cos(curr_theta)*haversin(new_phi-curr_phi)
			
			temp_arr = dist_arr<search_rad
			
			select_arr = select_arr+temp_arr
			
			p_bar.update(1)
			p_bar.refresh()
		
		select_arr = select_arr.astype(np.bool)
		
		shellnums_toadd = np.ones(new_theta[select_arr].size, dtype = int)*curr_shell
		vels_toadd = vel_arr[shellnums_toadd]
		
		theta = np.append(theta, new_theta[select_arr])
		phi = np.append(phi, new_phi[select_arr])
		shellnums = np.append(shellnums, shellnums_toadd)
		vels = np.append(vels, vels_toadd)
	
	p_bar.close()

	# # calculate mean anomaly for t_explode
	M_explode = u.rad * 2*np.pi*(t_start-tperi_i)/P_i
	# calculate the true anomaly for this mean anomaly
	(E_explode, f_i) = kepler_solve(e_i, M_explode)
	
	#calculate new orbital elements
	P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime = kick_kep_elements(vels, theta, phi, P_i, f_i, a_i, e_i, omega_i, anode_i, M, m)
	
	#calculate orbits, save data:
	t_post = np.linspace(t_start, t_start+delta_t, n_steps)
	
	# now we want to get the periastron times for all the exploded particle orbits...
	# calculate eccentric anomaly from true anomaly
	E_prime = np.arctan2(np.sqrt(1 - e_prime**2) * np.sin(f_prime), e_prime + np.cos(f_prime))
	
	# calculate the mean anomaly fom the eccentric anomaly
	M_prime = E_prime - e_prime*np.sin(E_prime)*u.rad
	
	tperi_prime = t_start - P_prime*M_prime/(2*np.pi*u.rad)
	
	files = np.zeros(P_prime.size, dtype = str)
	
	print('Simulating...')
	for idx in tqdm(range(P_prime.size)):
		shellnum = shellnums[idx]
		currtheta = theta[idx]
		currphi = phi[idx]
		_,_, Xs, Ys, Zs, _,_,_ = kep3d(t_post,P_prime[idx],tperi_prime[idx],a_prime[idx],e_prime[idx],inc_prime[idx],omega_prime[idx],anode_prime[idx])
		
		filename = 'cloud'+str(idx)
		files[idx] = filename
		save_kep3d_npz(Xs, Ys, Zs, shellnum, currtheta, currphi, filename = filename, overwrite=True, filepath = save_dir)
 		
 		
	return step_rate
##generally depreciated

# def simulate_spherical_old(n_shells = 3, kick_vel:u.m/u.s = 1.*u.km/u.s, P_i:u.year = 165.*u.year, 
# 		tperi_i:u.year = 2000 * u.year, e_i = 0., omega_i:u.rad = 0*u.rad, anode_i:u.rad = 0*u.rad, 
# 		M:u.kg = 1.0*u.Msun, m:u.kg = 1.0*u.Mearth, num_per_shell = 400, t_start:u.year = 2000*u.year, 
# 		delta_t:u.year = 50*u.year, n_steps = 500, show_pre = False, save_dir = './data/', 
# 		overwrite = False):
# 	'''creates a simulation for a spherically symmetric explosion
# 		
# 		inputs:
# 			rel_mass:    (numpy array) The relative masses of each spherical shell. The number of shells
# 					     is determined by the length of the array
# 			vel: 	     (numpy array or float) Determines the kick velocities for each shell. 
# 					     If an array, must be of same length as rel_mass, if an int, taken to
# 					     be the max velocity, and a linspace function is used to determine the
# 					     others
# 			P_i:	     Period of system pre-explosion
# 			f_i:	     true anomaly of system pre-explosion
# 			omega_i:     argument of periapsis of system pre-explosion
# 			anode_i:     longitude of ascending node of system pre-explosion
# 			M:		     primary (stellar) mass
# 			m:		     secondary (planitesimal) mass pre-explosion
# 			n_particles: the number of particles to simulate
# 			t_start:	 the time of the explosion
# 			delta_t:	 how much time to simulate after the explosion
# 			show_pre:	 Boolean or time, indicates whether you want to simulate the 
# 						 progenitor for some amount of time before the explosion, if a
# 						 number, that much time pre-explosion, if "True", defaults to 
# 						 five years
# 			save_dir:    directory to save output files to
# 			overwrite:   if it's ok for the program to overwrite the current directory at save_dir
# 			
# 	'''
# 	
# 	print('Warning: the function simulate_spherical is considered depreciated, proceed with caution')
# 	
# 	#calculate step rate (steps/year)
# 	step_rate = n_steps/delta_t
# 	
# 	#initial i is 0:
# 	inc_i = 0 *u.rad
# 	
# 	# given the masses and period, we really have to calculate a
# 	a_i = np.power(c.G*(M+m)*P_i*P_i / (4*np.pi*np.pi),1/3.).to(u.au)
# 	
# 	#calculate orbital velocity
# 	vorb=np.sqrt(c.G * (M+m)/ a_i).to(u.km/u.s)
# 	print('Initial orbital velocity: '+str(vorb))
# 	
# 	#create array of shell numbers:
# 	shellnums = np.zeros(n_shells*num_per_shell)
# 	for num in range(n_shells):
# 		shellnums[num*num_per_shell:(num+1)*num_per_shell] = num
# 	
# 	#print(shellnums)
# 	
# 	#prepare to save data:
# 	if os.path.exists(save_dir):
# 		if len(os.listdir(save_dir)) != 0: #if the directory is not empty
# 			if overwrite:
# 				empty_data_directory(save_dir)
# 			else:
# 				print('Save directory at '+save_dir+' is not empty. Select a new directory or set overwrite = True')
# 				return
# 	else:
# 		print('Save directory does not exist, creating new directory at '+save_dir)
# 		os.mkdir(save_dir)
# 		
# 	###PROGENITOR:
# 	#create single orbit:
# 	t_pro_orbit = np.linspace(t_start,t_start+P_i,(step_rate*P_i).astype(int))
# 	_, _, Xa, Ya, Za, _, _, _ = kep3d(t_pro_orbit, P_i, tperi_i, a_i, e_i, inc_i, omega_i, anode_i)
# 	
# 	save_kep3d_short(Xa, Ya, Za, filename = 'orbit.csv', overwrite=overwrite, filepath = save_dir)
# 	
# 	#to visualize progenitor:
# 	prog_steps = 0
# 	if isinstance(show_pre, bool):
# 		if show_pre:
# 			prog_steps = int(step_rate*7*u.year)
# 			prog_start = t_start-7*u.year
# 			t_prog = np.linspace(prog_start, t_start, prog_steps)
# 			
# 			_, _, Xs, Ys, Zs, _, _, _ = kep3d(t_prog, P_i, tperi_i, a_i, e_i, inc_i, omega_i, anode_i)
# 			save_kep3d_short(Xs, Ys, Zs, filename = 'progenitor.csv', overwrite=overwrite, post_steps = n_steps, filepath = save_dir)
# 			
# 	else:
# 		try:
# 			show_pre.to(u.year)
# 		except:
# 			print('show_pre must be a boolean or quantity with units of time')
# 			return
# 		
# 		prog_steps = step_rate*show_pre
# 		prog_start = t_start+show_pre
# 		t_prog = np.linspace(prog_start, t_start, prog_steps)
# 		
# 		_, _, Xs, Ys, Zs, _, _, _ = kep3d(t_prog, P_i, tperi_i, a_i, e_i, inc_i, omega_i, anode_i)
# 		save_kep3d_short(Xs, Ys, Zs, filename = 'progenitor.csv', overwrite=overwrite, post_steps = n_steps, filepath = save_dir)
# 	
# 	## Calculate shell velocities:
# 	if isinstance(kick_vel.value,np.ndarray): #if vel is passed in as an array
# 		if kick_vel.size != n_shells:
# 			print('error, kick_vel must be array of same size as rel_mass OR float')
# 			return
# 		else:
# 			vel_arr = kick_vel.copy()
# 	else:
# 		vel_arr = np.linspace(0, kick_vel, n_shells+1)[1:]
# 		
# 	##Assuming each particle is equal mass, calc number of particles per shell:
# 	#normalize massses
# # 	norm_mass = rel_mass/(rel_mass.sum())
# # 	num_per_shell = (norm_mass*n_particles).astype(int)
# 	
# #	num_per_shell = (n_particles/n_shells)
# 		
# 	##Create shells:
# 	r, theta, phi = fibonacci_sphere(n_shells, num_per_shell, vel_arr)	
# 	
# 	#calculate f_i# calculate mean anomaly for t_explode
# 	M_explode = u.rad * 2*np.pi*(t_start-tperi_i)/P_i
# 	# calculate the true anomaly for this mean anomaly
# 	(E_explode, f_i) = kepler_solve(e_i, M_explode)
# 	
# 	#calculate new orbital elements
# 	P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime = kick_kep_elements(r, theta, phi, P_i, f_i, a_i, e_i, omega_i, anode_i, M, m)
# 	
# 	#calculate orbits, save data:
# 	t_post = np.linspace(t_start, t_start+delta_t, n_steps)
# 	
# 	# now we want to get the periastron times for all the exploded particle orbits...
# 	# calculate eccentric anomaly from true anomaly
# 	E_prime = np.arctan2(np.sqrt(1 - e_prime**2) * np.sin(f_prime), e_prime + np.cos(f_prime))
# 	
# 	# calculate the mean anomaly fom the eccentric anomaly
# 	M_prime = E_prime - e_prime*np.sin(E_prime)*u.rad
# 	
# 	tperi_prime = t_start - P_prime*M_prime/(2*np.pi*u.rad)
# 	
# 	for idx in tqdm(range(P_prime.size)):
# 		shellnum = shellnums[idx]
# 		X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_post,P_prime[idx],tperi_prime[idx],a_prime[idx],e_prime[idx],inc_prime[idx],omega_prime[idx],anode_prime[idx])
# 		
# 		shellarr = np.ones_like(Xs)*shellnum
# 		
# 		filename = 'cloud'+str(idx)+'.csv'
# 		save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, shellarr, filename = filename, overwrite=True, pre_steps = prog_steps, filepath = save_dir)
#  		
#  		
# 	return 