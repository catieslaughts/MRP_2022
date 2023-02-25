import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import os
from fileio import *

@u.quantity_input
def single_lightcurve(directory = './data/', R:u.AU = .1 * u.AU ):
	'''determines the Earth-view lightcurve for a set of simulations (stored in the directory pointed 
		to by directory) around a star with radius R.
		
		Calculates a single light curve for all particles, regardless of which shell they are in
		
		Assumes Earth sits on the x axis in the positive direction
	'''
	
	#get file names:
	files = os.listdir(directory)
	
	#remove orbit and progenitor from list
	if 'orbit.csv' in files:
		files.remove('orbit.csv')
	
	#read in first file to get number of steps
	df = read_kep3d_pd(directory, files[0])
	
	#create array to store whether the given particle is crossing at the time in question
	crossing_array = np.zeros((df.shape[0], len(files)))
	
	for idx, file in enumerate(files):
		df = read_kep3d_pd(directory, file)
		
		X = np.asarray(df['Xs']) * u.au
		Y = np.asarray(df['Ys']) * u.au
		Z = np.asarray(df['Zs']) * u.au
		
		radii = np.sqrt(Y*Y + Z*Z)
		
		#get rid of particles behind the star 
		#(possibly redundant for our cases, but wanted to be safe)
		radii[X<0] = 1000000000*R #arbitrarily large
	
		crossing_array[:,idx] = (radii<R)
	
	lightcurve = np.count_nonzero(crossing_array, 1)
	
	plt.plot(1 - lightcurve)
	plt.yticks([])
	plt.xlabel('Time Step')
	
	plt.title('Lightcurve with shells of equal weight')
	
	plt.show()
	#print(np.count_nonzero(crossing_array, 1))
	
	return


def weighted_lightcurve(directory = './data/', R:u.AU = .1 * u.AU, num_shells = 4, shell_weights = [1, 1, 1, 1], plt_subs = True, plt_total = True, save_plt = False, savefile = './lightcurve.pdf'):
	'''determines the Earth-view lightcurve for a set of simulations (stored in the directory pointed 
		to by directory) around a star with radius R.
		
		Calculates a single light curve for all particles, regardless of which shell they are in
		
		Assumes Earth sits on the x axis in the positive direction
	'''
	
	# 
	shell_weights = np.asarray(shell_weights)
	shellnums = np.arange(num_shells)
	
	if shell_weights.size != num_shells:
		print('Error: the length of the shell weight array should be the same as the number of shells')
		return
	
	#get file names:
	files = os.listdir(directory)
	
	#remove orbit and progenitor from list
	if 'orbit.csv' in files:
		files.remove('orbit.csv')
	if 'progenitor.csv' in files:
		files.remove('progenitor.csv')
	
	#read in first file to get number of steps
	df = read_kep3d_pd(directory, files[0])
	
	#create array to store whether the given particle is crossing at the time in question
	crossing_array = np.zeros((df.shape[0], len(files), num_shells))
	
	for idx, file in enumerate(files):
		df = read_kep3d_pd(directory, file)
		
		shell_num = int(df['shellnum'][0])
		
		X = np.asarray(df['Xs']) * u.au
		Y = np.asarray(df['Ys']) * u.au
		Z = np.asarray(df['Zs']) * u.au
		
		radii = np.sqrt(Y*Y + Z*Z)
		
		#get rid of particles behind the star 
		#(possibly redundant for our cases, but wanted to be safe)
		radii[X<0] = 1000000000*R #arbitrarily large
	
		crossing_array[:,idx,shell_num] = (radii<R)
	
	sub_lcs = np.count_nonzero(crossing_array[:,:,:], 1)
	
	total_lc = np.dot(sub_lcs,shell_weights)
	
	plt.plot(-total_lc, label = 'Total Lightcurve')
	
	#print(total_lc)
	if plt_subs:
		plt.plot(-sub_lcs, label=shellnums)
	
	plt.legend()
	
	plt.yticks([])
	plt.xlabel('Time Step')
	plt.xlim(150,350)
	
	plt.title('Light Curve -- [2,1,1,1]')
	
	if save_plt:
		#print('Saving...')
		plt.savefig(savefile, bbox_inches = 'tight')
	
	plt.show()
	#print(np.count_nonzero(crossing_array, 1))
	
	return


weighted_lightcurve(directory = './test_subset/', shell_weights = [2,1,1,1], plt_subs = False, save_plt = True)




	