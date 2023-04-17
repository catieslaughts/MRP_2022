import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import os
from fileio import *
from tqdm import tqdm

@u.quantity_input
def save_lightcurve_data(directory = './data/', paramfile = 'params.csv', cut_run = False, foi_file = 'foi.csv', lc_file = 'lc_data'):
	
	save_foi = True
	save_countarr = True
	
	_,*_,n_steps,n_pre_steps,num_shells,R_star = read_param_file(paramfile)
	
	if cut_run:
		R_star = 3*R_star
	
	shellnums = np.arange(num_shells)
	
	#get file names:
	files = np.asarray(os.listdir(directory))
	
	#remove orbit and progenitor from list
	if 'orbit.npz' in files:
		files.remove('orbit.npz')
	if 'progenitor.npz' in files:
		files.remove('progenitor.npz')
	
	#create array to store whether the given particle is crossing at the time in question
	crossing_array = np.zeros((n_steps+n_pre_steps, len(files), num_shells))
	
	if save_foi:
		files_of_interest = []
	
	print('Creating lightcurves...')
	for idx in tqdm(range(files.size)):
		file = files[idx]
		X, Y, Z, shell_num, *_ = read_kep3d_npz(directory, file)
		
		radii = np.sqrt(Y*Y + Z*Z)
		
		#get rid of particles behind the star 
		#(possibly redundant for our cases, but wanted to be safe)
		radii[X<0] = 1000000000*R_star #arbitrarily large
		
		crossing_array[:,idx,shell_num] = (radii<R_star)
		
		#print(radii.shape)
		
		if save_foi:
			if np.count_nonzero(radii<R_star)>0:
				files_of_interest.append(file)
			
	if save_foi:
		foi = np.asarray(files_of_interest)
		np.savetxt(foi_file, foi, delimiter=',', fmt='%s')
		print('foi saved')
	
	lcs = np.count_nonzero(crossing_array[:,:,:], 1)
	
	if save_countarr:
		np.savez(lc_file, lcs = lcs)
		print('sub lightcurves saved')

@u.quantity_input
def plot_from_saved(readfile = 'lc_data.npz', shell_weights = np.nan, steprate:1/u.year = np.nan/u.year, plt_subs = True, plt_total = True, save_plt = False, savefile = './lightcurve.pdf'):
	
	##Read in from file:
	
	sub_lcs = np.load(readfile)['lcs']
	
	if np.isnan(shell_weights):
		shell_weights = np.ones(sub_lcs.shape[1])
	else:
		shell_weights = np.asarray(shell_weights)
	shellnums = np.arange(shell_weights.size)
	
	if sub_lcs.shape[1] != shell_weights.size:
		print('Error: the length of the shell weight array should be the same as the number of shells')
		return
	
	total_lc = np.dot(sub_lcs,shell_weights)
	
	xaxis = np.arange(total_lc.size)
	
	if not np.isnan(steprate):
		xaxis = xaxis * (1/steprate).to(u.day)
	
	if plt_total:
		plt.plot(xaxis, -total_lc, label = 'Total Lightcurve')
	
	#print(total_lc)
	if plt_subs:
		labels = np.ones(shellnums.size, dtype = str)
		for idx, label in enumerate(labels):
			label = 'shell: '+str(shellnums[idx])+' weight: '+str(shell_weights[idx])
		plt.plot(xaxis, -sub_lcs, label=labels)
	
	plt.legend()
	
	plt.yticks([])
	plt.xlabel('Time Step')
	#plt.xlim(150,350)
	
	plt.title('Light Curve')
	
	if save_plt:
		#print('Saving...')
		plt.savefig(savefile, bbox_inches = 'tight')
	
	plt.show()
	#print(np.count_nonzero(crossing_array, 1))
	
	return
	
	
#depreciated

# def single_lightcurve(directory = './data/', R:u.AU = .1 * u.AU ):
# 	'''determines the Earth-view lightcurve for a set of simulations (stored in the directory pointed 
# 		to by directory) around a star with radius R.
# 		
# 		Calculates a single light curve for all particles, regardless of which shell they are in
# 		
# 		Assumes Earth sits on the x axis in the positive direction
# 	'''
# 	
# 	print('Warning: the function single_lightcurve is considered depreciated, proceed with caution')
# 	#get file names:
# 	files = os.listdir(directory)
# 	
# 	#remove orbit and progenitor from list
# 	if 'orbit.csv' in files:
# 		files.remove('orbit.csv')
# 	
# 	#read in first file to get number of steps
# 	df = read_kep3d_pd(directory, files[0])
# 	
# 	#create array to store whether the given particle is crossing at the time in question
# 	crossing_array = np.zeros((df.shape[0], len(files)))
# 	
# 	for idx, file in enumerate(files):
# 		df = read_kep3d_pd(directory, file)
# 		
# 		X = np.asarray(df['Xs']) * u.au
# 		Y = np.asarray(df['Ys']) * u.au
# 		Z = np.asarray(df['Zs']) * u.au
# 		
# 		radii = np.sqrt(Y*Y + Z*Z)
# 		
# 		#get rid of particles behind the star 
# 		#(possibly redundant for our cases, but wanted to be safe)
# 		radii[X<0] = 1000000000*R #arbitrarily large
# 	
# 		crossing_array[:,idx] = (radii<R)
# 	
# 	lightcurve = np.count_nonzero(crossing_array, 1)
# 	
# 	plt.plot(1 - lightcurve)
# 	plt.yticks([])
# 	plt.xlabel('Time Step')
# 	
# 	plt.title('Lightcurve with shells of equal weight')
# 	
# 	plt.show()
# 	#print(np.count_nonzero(crossing_array, 1))
# 	
# 	return

# def weighted_lightcurve(directory = './data/', R:u.AU = .1 * u.AU, num_shells = 4, shell_weights = [1, 1, 1, 1], steprate:1/u.year = np.nan/u.year, plt_subs = True, plt_total = True, save_plt = False, savefile = './lightcurve.pdf'):
# 	'''determines the Earth-view lightcurve for a set of simulations (stored in the directory pointed 
# 		to by directory) around a star with radius R.
# 		
# 		Calculates a single light curve for all particles, regardless of which shell they are in
# 		
# 		Assumes Earth sits on the x axis in the positive direction
# 	'''
# 	
# 	print('Warning: the function weighted_lightcurve is considered depreciated, proceed with caution')
# 	
# 	save_foi = True
# 	save_countarr = True
# 	
# 	# 
# 	shell_weights = np.asarray(shell_weights)
# 	shellnums = np.arange(num_shells)
# 	
# 	if shell_weights.size != num_shells:
# 		print('Error: the length of the shell weight array should be the same as the number of shells')
# 		return
# 	
# 	#get file names:
# 	files = os.listdir(directory)
# 	
# 	#remove orbit and progenitor from list
# 	if 'orbit.csv' in files:
# 		files.remove('orbit.csv')
# 	if 'progenitor.csv' in files:
# 		files.remove('progenitor.csv')
# 	
# 	#read in first file to get number of steps
# 	df = read_kep3d_pd(directory, files[0])
# 	
# 	#create array to store whether the given particle is crossing at the time in question
# 	crossing_array = np.zeros((df.shape[0], len(files), num_shells))
# 	
# 	if save_foi:
# 		files_of_interest = []
# 	
# 	for idx, file in enumerate(files):
# 		df = read_kep3d_pd(directory, file)
# 		
# 		shell_num = int(df['shellnum'][0])
# 		
# 		X = np.asarray(df['Xs']) * u.au
# 		Y = np.asarray(df['Ys']) * u.au
# 		Z = np.asarray(df['Zs']) * u.au
# 		
# 		radii = np.sqrt(Y*Y + Z*Z)
# 		
# 		#get rid of particles behind the star 
# 		#(possibly redundant for our cases, but wanted to be safe)
# 		radii[X<0] = 1000000000*R #arbitrarily large
# 	
# 		crossing_array[:,idx,shell_num] = (radii<R)
# 		
# 		if save_foi:
# 			if np.count_nonzero(radii<R)>0:
# 				files_of_interest.append(file)
# 			
# 	if save_foi:
# 		foi = np.asarray(files_of_interest)
# 		np.savetxt('foi.csv', foi, delimiter=',', fmt='%s')
# 		print('foi saved')
# 	
# 	sub_lcs = np.count_nonzero(crossing_array[:,:,:], 1)
# 	
# 	if save_countarr:
# 		np.savetxt('lightcurves.csv', sub_lcs, delimiter=',', fmt = '%d')
# 		print('sub lightcurves saved')
# 	
# 	total_lc = np.dot(sub_lcs,shell_weights)
# 	
# 	xaxis = np.arange(total_lc.size)
# 	
# 	if not np.isnan(steprate):
# 		xaxis = xaxis * (1/steprate).to(u.day)
# 	
# 	if plt_total:
# 		plt.plot(xaxis, -total_lc, label = 'Total Lightcurve')
# 	
# 	#print(total_lc)
# 	if plt_subs:
# 		plt.plot(xaxis, -sub_lcs, label=shellnums)
# 	
# 	plt.legend()
# 	
# 	plt.yticks([])
# 	plt.xlabel('Time Step')
# 	#plt.xlim(150,350)
# 	
# 	plt.title('Light Curve')
# 	
# 	if save_plt:
# 		#print('Saving...')
# 		plt.savefig(savefile, bbox_inches = 'tight')
# 	
# 	plt.show()
# 	#print(np.count_nonzero(crossing_array, 1))
# 	
# 	return






	