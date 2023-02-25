import numpy as np 
import astropy.units as u
import matplotlib.pyplot as plt

#tools for making and visualizing spherical particle distributions, as well as generic
#conversion function from spherical to cartesian coordinates

##create point spheres:
@u.quantity_input
def random_sphere(n_layers, n_particles, vels:u.m/u.s):
	'''creates concentric spheres of points with randomly selected locations
		
		Inputs:
			n_layers (int): how many concentric spherical shells to makes
			n_particles (int): how many particles PER SHELL
			max_vel (float): velocity (radius) of the outermost shell
			
		Outputs:
			r (array): the radii (velocities) for each particle
			theta (array): the theta coordinate of each particle
			phi (array): the phi coordinate of each particle
	
	'''
	r = np.zeros([n_layers, n_particles], dtype = float)
	
	vels = np.linspace(max_vel, 0, n_layers, endpoint=False)[::-1] #avoids a layer with 0 velocity
	for idx, vel in enumerate(vels):
		r[idx] = np.ones(n_particles)*vel
	
	theta = np.random.rand(n_layers, n_particles)*np.pi #randomly selected from 0 to 2pi
	phi = np.random.rand(n_layers, n_particles)*2*np.pi #randomly selected from 0 to pi
	
	#units
	r=r*max_vel.unit
	theta = theta*u.rad
	phi = phi*u.rad
	
	return r, theta, phi

def fibonacci_sphere(n_layers, n_particles, vels:u.m/u.s):
	'''makes use of the fibonacci sphere algorithm to create several circular layers of points. 
	The radius output by this function will be used as our velocities in our simulation
	Sources:
		http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
		https://arxiv.org/pdf/0912.4540.pdf
		
		Inputs:
			n_layers (int): how many concentric spherical shells to makes
			n_particles (int or array of ints): how many particles PER SHELL
			vels (array of floats): velocities of each shell
			
		Outputs:
			r (array): the radii (velocities) for each particle
			theta (array): the theta coordinate of each particle
			phi (array): the phi coordinate of each particle

	'''
	
	if not isinstance(n_particles,(np.ndarray)):
		n_particles = np.ones(n_layers, dtype = int)*n_particles
	
# 	r = np.zeros(n_layers*n_particles.sum(), dtype = float)
# 	theta = np.zeros(n_layers*n_particles.sum(), dtype = float)
# 	phi = np.zeros(n_layers*n_particles.sum(), dtype = float)
	
	gr = (1 + np.sqrt(5))/2
	
	#vels = np.linspace(max_vel, 0, n_layers, endpoint=False)[::-1] #avoids a layer with 0 velocity
	
	for idx, vel in enumerate(vels):
		parts = np.arange(0,n_particles[idx])
		
		if idx == 0:
			r = np.ones(n_particles[idx])*vel
			phi = 2*np.pi * parts / gr
			theta = np.arccos(1 - 2*(parts+0.5)/n_particles[idx])
			#theta and phi are flipped from the given sources because the conventions are flipped in math and physics
		else:
			r_temp = np.ones(n_particles[idx])*vel
			phi_temp = 2*np.pi * parts / gr
			theta_temp = np.arccos(1 - 2*(parts+0.5)/n_particles[idx])
			
			r = np.append(r, r_temp)
			phi = np.append(phi, phi_temp)
			theta = np.append(theta, theta_temp)
		
		#print(r.size)
	
	#units
	r=r#*vels.unit
	theta = theta*u.rad
	phi = phi*u.rad
	
	return r, theta, phi
	
def latlong_sphere(n_layers, n_particles, max_vel:u.m/u.s):
	'''ummmm'''
	
	r = np.zeros([n_layers, n_particles], dtype = float)
	theta = np.zeros([n_layers, n_particles], dtype = float)
	phi = np.zeros([n_layers, n_particles], dtype = float)
	pass

#coordinate conversion
def sph2cart(r, theta, phi):
	'''eqs 4.4.1-3 in 
	https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Book%3A_Electromagnetics_I_(Ellingson)/04%3A_Vector_Analysis/4.04%3A_Spherical_Coordinates'''
	
	x = r * np.cos(phi) * np.sin(theta)
	y = r * np.sin(phi) * np.sin(theta)
	z = r * np.cos(theta)
	
	return x, y, z
	
#plotting	
def draw_spherical(r,theta,phi):
	'''makes a 3d plot of the points stored in r,theta,phi
	for our purposes, r really dictates velocity, not radius, but we plot it here as radius for viewing convenience'''
	
	d3 = True
	d2 = False 
	projections = False
	
	pointsize = 15
	
	x,y,z = sph2cart(r.value, theta.value, phi.value)
	
	if d3:
		ax = plt.axes(projection='3d')
		
		if len(x.shape) > 1:
			for idx in range(x.shape[0]):
				ax.scatter3D(x[idx], y[idx], z[idx], s=pointsize//3)
		else:
			ax.scatter3D(x, y, z, s=pointsize//3)
		
		#set_axes_equal(ax)
		plt.show()
	
	if d2:
		labels = [['x','y'],['y','z'],['x','z']]
		
		for lab_idx, pair in enumerate([[x,y],[y,z],[x,z]]):
			ax= plt.axes()
			if len(x.shape) > 1:
				for idx in range(x.shape[0]):
					ax.scatter(pair[0][idx], pair[1][idx], s=pointsize)
			else:
				ax.scatter(pair[0], pair[1], s=pointsize)
			
			plt.title('Data Projection')
			plt.xlabel(labels[lab_idx][0])
			plt.ylabel(labels[lab_idx][1])
			ax.set_aspect('equal')
			plt.show()
			
	if projections:
		labels = [['x','y'],['y','z'],['x','z']]
		
		for lab_idx, pair in enumerate([[x,y,z],[y,z,x],[x,z,y]]):
			ax= plt.axes()
			
			max = x.shape[1]//2
			print(max)
			
			idx_array = np.argsort(pair[2])
			idx_array = idx_array[0:max]
			
			print(idx_array.shape)
			
			sub_first = pair[0][idx_array]
			sub_second= pair[1][idx_array]
			
			sub_first = np.reshape(sub_first, [x.shape[0], sub_first.size//x.shape[0]])
			sub_second = np.reshape(sub_second, [x.shape[0], sub_second.size//x.shape[0]])
			
			if len(x.shape) > 1:
				for idx in range(x.shape[0]):
					ax.scatter(sub_first[idx], sub_second[idx], s=pointsize)
			else:
				ax.scatter(sub_first, sub_second, s=pointsize)
			
			plt.title('Orthographic Projection')
			plt.xlabel(labels[lab_idx][0])
			plt.ylabel(labels[lab_idx][1])
			ax.set_aspect('equal')
			plt.show()
	
	pass	
 