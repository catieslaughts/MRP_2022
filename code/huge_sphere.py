import numpy as np
import astropy.units as u
from spheres import fibonacci_sphere, draw_spherical
import time


def make_huge_sphere(num_parts = 1e6, file = 'giantsphere', largenumberoverride = False):
	
	if num_parts > 1e8:
		if not largenumberoverride:
			print('Warning, number of particles exceeds the safety cutoff (100 million), if this is intentional, set largenumberoverride = True')
			return
	
	num_parts = int(num_parts)
	
# 	print('Running...')
# 	start = time.time()
	r, theta, phi = fibonacci_sphere(1, num_parts, [1]*u.km/u.s)
	
	theta = theta.to(u.rad).value
	phi = phi.to(u.rad).value
	
# 	print('Run time: '+str(time.time()-start)+' seconds')
	
# 	print('Saving...')
	np.savez(file, theta=theta, phi=phi)
	
# 	print('Total time: '+str(time.time()-start)+' seconds')
	
def get_giant_sphere(file = 'giantsphere.npz'):
	
	data = np.load(file)
	
	theta, phi = data['theta'], data['phi']
	data.close()
	
	theta = theta* u.rad
	phi = phi* u.rad
	
	return(theta, phi)


def haversin(angle):
	
	return np.power(np.sin(angle/2.),2)