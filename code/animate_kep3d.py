import numpy as np
from fileio import *
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

from spheres import draw_spherical
from tqdm import tqdm

#Functions for animating the outputs of kep3d, position and velocity data can be saved from kep3d outputs
#using file write functions in fileio.py. Each steps through all the inidvidual files in a given directory,
#reading each and adding the corresponding point to the animation
		
def animate_los(directory = './data/', save = False, save_file = './LOSanimation.gif', paramfile = 'params.csv'):

	colors = ['orange', 'green', 'red', 'purple', 'blue'] #outermost shell is colored by the last one
	writer = animation.PillowWriter()
	
	_,*_,R_star = read_param_file(paramfile)
	R_star = R_star.to(u.au).value

	Ys = []#setting up orbit data, using regular lists because they append more nicely for our purposes
	Zs = []
	vel = []
	shell_arr = []
	
	showpre = False
	
	files = os.listdir(directory)
	
	#files.reverse()
	for filename in tqdm(files):
		#for making animations which show the original progenitor and/or orbit
		if not filename.startswith('.'):
			if 'orbit' in filename:
				continue
			
			try:	
				Xs_temp, Ys_temp, Zs_temp, shellnum, *_ = read_kep3d_npz(directory, filename)
			except:
				print(filename)
			
			Ys.append(Ys_temp)
			Zs.append(Zs_temp)
			
			shell_arr.append(shellnum)
	
	Ys = np.asarray(Ys)#convert regular lists to np arrays
	Zs = np.asarray(Zs)
	
	shell_arr = np.asarray(shell_arr,dtype=int)
	
	##sort so shells are plotted in order
	idxs = np.arange(shell_arr.size)
	idxs = [q for _, q in sorted(zip(shell_arr,idxs))]
	
	Ys = Ys[idxs]
	Zs = Zs[idxs]
	
	shell_arr = shell_arr[idxs]
	
	#plot
	fig, ax = plt.subplots()#set up figure
	
	ax.set_xlabel('Y')
	ax.set_ylabel('Z')
	ax.set_title('RIP Alderaan')
	
	
# 	offset = .1
# 	ax.set_xlim((Ys.min()-offset,Ys.max()+offset)) #axis shape
# 	ax.set_ylim((Zs.min()-offset, Zs.max()+offset))
	
	circle1 = plt.Circle((0, 0), R_star, color='black')
	ax.add_patch(circle1) #primary
	
	ax.axis('equal')
	ax.set_xlim((-5,5)) #axis shape
	
# 	for i in range(Xs.shape[0]):
# 		ax.plot(Xs[i],Ys[i], Zs[i], linestyle='--', alpha=.5) #if you want a reference orbit, that can be plotted
# 	ax.set_prop_cycle(None) #resets color cycler so orbits match up with particles
	
	points = np.ndarray(Ys.shape[0], dtype=object) #arbitrary length based on # of filenames
	
	for idx1 in range(Ys.shape[0]):
		point, = ax.plot(Ys[idx1][0], Zs[idx1][0]) #sets up 2Dline objects
		points[idx1] = point
	
	def animate(frame):
		'''animation function for an arbitrary number of particles'''
		returns = []
		
		for idx2, point in enumerate(points):#steps through list of particles
			point.set_xdata(Ys[idx2][frame:frame+1])
			point.set_ydata(Zs[idx2][frame:frame+1])
			
			if shell_arr[idx2] == 1000: #progenitor
				point.set_color('grey')
				point.set_markersize(2)
			else:
				point.set_color(colors[shell_arr[idx2]])
				point.set_markersize(1)
			
			point.set_marker('o')
			
			returns.append(point)
		
		#print(returns)
		return returns
	
	ani = animation.FuncAnimation(fig, animate, frames = np.arange(1, Ys.shape[1]), interval = 1, blit=True, repeat = True)
	plt.show()
	
	if save:
		ani.save(save_file)

def animate_los_subset(foi_list = 'foi.csv',directory = './data/', save = False, save_file = './LOSanimation.gif', paramfile = 'params.csv'):

	colors = ['orange', 'green', 'red', 'purple', 'blue'] #outermost shell is colored by the last one
	writer = animation.PillowWriter()
	
	_,*_,n_pre_steps,_,R_star = read_param_file(paramfile)
	
	R_star = R_star.to(u.au).value

	Ys = []#setting up orbit data, using regular lists because they append more nicely for our purposes
	Zs = []
	vel = []
	shell_arr = []
	
	files = np.loadtxt('foi.csv', dtype = 'str', delimiter=',')
	if files.size == 1:
		files = np.asarray([files])
	
	#files.reverse()
	for filename in tqdm(files):
		if not filename.startswith('.'):
			try:	
				Xs_temp, Ys_temp, Zs_temp, shellnum, *_ = read_kep3d_npz(directory, filename)
			except:
				print(filename)
			
			
			
			Ys.append(Ys_temp[n_pre_steps:])
			Zs.append(Zs_temp[n_pre_steps:])
			#if you want to change this code to show the progenitor later you'll need to get rid of the slicing on the above lines
			
			shell_arr.append(shellnum)
				
	
	Ys = np.asarray(Ys)#convert regular lists to np arrays
	Zs = np.asarray(Zs) 
	
	shell_arr = np.asarray(shell_arr,dtype=int)
	
	##sort so shells are plotted in order
	idxs = np.arange(shell_arr.size)
	idxs = [q for _, q in sorted(zip(shell_arr,idxs))]
	
	Ys = Ys[idxs]
	Zs = Zs[idxs]
	
	shell_arr = shell_arr[idxs]
	
	#plot
	fig, ax = plt.subplots()#set up figure
	
	ax.set_xlabel('Y')
	ax.set_ylabel('Z')
	ax.set_title('RIP Alderaan')
	
	ax.axis('equal')
	ax.set_xlim((-3,3))
	#ax.set_ylim((-3,3))
	
	
	circle1 = plt.Circle((0, 0), R_star, color='black')
	ax.add_patch(circle1) #primary
	
	points = np.ndarray(Ys.shape[0], dtype=object) #arbitrary length based on # of filenames
	
	for idx1 in range(Ys.shape[0]):
		point, = ax.plot(Ys[idx1][0], Zs[idx1][0]) #sets up 2Dline objects
		points[idx1] = point
	
	def animate(frame):
		'''animation function for an arbitrary number of particles'''
		returns = []
		
		for idx2, point in enumerate(points):#steps through list of particles
			point.set_xdata(Ys[idx2][frame:frame+1])
			point.set_ydata(Zs[idx2][frame:frame+1])
			
			if shell_arr[idx2] > 10 : #progenitor
				print("ASDFADSFA")
				point.set_color('grey')
				point.set_markersize(2)
			else:
				#print("WORK")
				point.set_color(colors[shell_arr[idx2]])
				point.set_markersize(1)
			
			point.set_marker('o')
			
			returns.append(point)
		
		#print(returns)
		return returns
	
	ani = animation.FuncAnimation(fig, animate, frames = np.arange(1, Ys.shape[1]),interval = 10, blit=True, repeat = True)
	plt.show()
	
	if save:
		ani.save(save_file)

def animate_3d(directory = './data/', save = False, save_file = './3Danimation.gif', paramfile = 'params.csv'):
	'''sets up, runs, and saves a 3d animation
	
	Args:
		filepath: string, path to the directory with data to be read in
		save: boolean, if True, saves animation to .gif in pwd
		
	Returns:none
	'''
	colors = ['orange', 'green', 'red', 'purple', 'blue'] #outermost shell is colored by the last one
	writer = animation.PillowWriter()
	
	_,*_,R_star = read_param_file(paramfile)
	R_star = R_star.to(u.au).value

	Xs = [] #setting up orbit data, using regular lists because they append more nicely for our purposes
	Ys = []
	Zs = []
	shell_arr = []
	
	files = os.listdir(directory)
	#files.reverse()
	orbitflag = 0
	for filename in files:
		if not filename.startswith('.'):
			if 'orbit' in filename: # a file showing a single orbit of the progenitor
				orbitflag = 1
				orbitX, orbitY, orbitZ, *_ = read_kep3d_npz(directory, filename)
				
			else:
				Xs_temp, Ys_temp, Zs_temp, shellnum, *_ = read_kep3d_npz(directory, filename)
				Xs.append(Xs_temp)#adds to list
				Ys.append(Ys_temp)
				Zs.append(Zs_temp)
				
				shell_arr.append(shellnum)
	
	Xs = np.asarray(Xs)#convert regular lists to np arrays
	Ys = np.asarray(Ys)
	Zs = np.asarray(Zs)
	
	print(Xs.max())
	
	shell_arr = np.asarray(shell_arr,dtype=int)
	
	##sort so shells are plotted in order
	idxs = np.arange(shell_arr.size)
	idxs = [q for _, q in sorted(zip(shell_arr,idxs))]
	
	Xs = Xs[idxs]
	Ys = Ys[idxs]
	Zs = Zs[idxs]
	
	shell_arr = shell_arr[idxs]
	
	#plot
	fig, ax = plt.subplots(figsize=(8,7),subplot_kw={'projection': '3d'})#set up figure
	
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	ax.set_title('RIP Alderaan')
	
	h = 40
	ax.set_xlim(-h,h)
	ax.set_ylim(-h,h)
	ax.set_zlim(-h,h)
	
	ax.scatter(0., 0., 0., marker='o', color='black') #primary
	
	if orbitflag == 1:
		ax.plot(orbitX,orbitY, orbitZ, linestyle='-', alpha=.5, c = 'grey')
	
# 	for i in range(Xs.shape[0]):
# 		ax.plot(Xs[i],Ys[i], Zs[i], linestyle='--', alpha=.5) #if you want a reference orbit, that can be plotted
# 	ax.set_prop_cycle(None) #resets color cycler so orbits match up with particles
	
	points = np.ndarray(Xs.shape[0], dtype=object) #arbitrary length based on # of filenames
	
	for idx1 in range(Xs.shape[0]):
		point, = ax.plot(Xs[idx1][0], Ys[idx1][0], Zs[idx1][0]) #sets up 2Dline objects
		points[idx1] = point
	
	def animate(frame):
		'''animation function for an arbitrary number of particles'''
		returns = []
		
		for idx2, point in enumerate(points):#steps through list of particles
			point.set_xdata(Xs[idx2][frame:frame+1])
			point.set_ydata(Ys[idx2][frame:frame+1])
			point.set_3d_properties(Zs[idx2][frame:frame+1])  # update the data for each
			
			if shell_arr[idx2] == 1000: #progenitor
				point.set_color('grey')
				point.set_markersize(2)
			else:
				point.set_color(colors[shell_arr[idx2]])
				point.set_markersize(1)
			
			point.set_marker('o')
			
			returns.append(point)
		
		#print(returns)
		return returns
	
	ani = animation.FuncAnimation(fig, animate, frames = np.arange(1, Xs.shape[1]), interval = 1, blit=True, repeat = True)
	plt.show()
	
	if save:
		ani.save(save_file)

def animate_3d_subset(foi_list = 'foi.csv', directory = './data/', save = False, save_file = './3Danimation.gif', paramfile = 'params.csv'):
	'''sets up, runs, and saves a 3d animation
	
	Args:
		filepath: string, path to the directory with data to be read in
		save: boolean, if True, saves animation to .gif in pwd
		
	Returns:none
	'''
	colors = ['orange', 'green', 'red', 'purple'] #outermost shell is colored by the last one
	writer = animation.PillowWriter()
	
	_,*_,R_star = read_param_file(paramfile)
	R_star = R_star.to(u.au).value

	Xs = [] #setting up orbit data, using regular lists because they append more nicely for our purposes
	Ys = []
	Zs = []
	vel = []
	shell_arr = []
	
	files = np.loadtxt('foi.csv', dtype = 'str', delimiter=',')
	#files.reverse()
	orbitflag = 0
	for filename in files:
		#for making animations which show the original progenitor and/or orbit
		if not filename.startswith('.'):
			if 'orbit' in filename: # a file showing a single orbit of the progenitor
				orbitflag = 1
				orbitX, orbitY, orbitZ, *_ = read_kep3d_npz(directory, filename)
				
			else:
				Xs_temp, Ys_temp, Zs_temp, shellnum, *_ = read_kep3d_npz(directory, filename)
				Xs.append(Xs_temp)#adds to list
				Ys.append(Ys_temp)
				Zs.append(Zs_temp)
				
				shell_arr.append(shellnum)
	
	Xs = np.asarray(Xs)#convert regular lists to np arrays
	Ys = np.asarray(Ys)
	Zs = np.asarray(Zs)
	
	shell_arr = np.asarray(shell_arr,dtype=int)
	
	##sort so shells are plotted in order
	idxs = np.arange(shell_arr.size)
	idxs = [q for _, q in sorted(zip(shell_arr,idxs))]
	
	Xs = Xs[idxs]
	Ys = Ys[idxs]
	Zs = Zs[idxs]
	
	shell_arr = shell_arr[idxs]
	
	#plot
	fig, ax = plt.subplots(figsize=(8,7),subplot_kw={'projection': '3d'})#set up figure
	
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	ax.set_title('RIP Alderaan')
	
	h = 50
	ax.set_xlim(-h,h)
	ax.set_ylim(-h,h)
	ax.set_zlim(h,-h)
	
	ax.scatter(0., 0., 0., marker='o', color='black') #primary
	
	if orbitflag == 1:
		ax.plot(orbitX,orbitY, orbitZ, linestyle='-', alpha=.5, c = 'grey')
	
# 	for i in range(Xs.shape[0]):
# 		ax.plot(Xs[i],Ys[i], Zs[i], linestyle='--', alpha=.5) #if you want a reference orbit, that can be plotted
# 	ax.set_prop_cycle(None) #resets color cycler so orbits match up with particles
	
	points = np.ndarray(Xs.shape[0], dtype=object) #arbitrary length based on # of filenames
	
	for idx1 in range(Xs.shape[0]):
		point, = ax.plot(Xs[idx1][0], Ys[idx1][0], Zs[idx1][0]) #sets up 2Dline objects
		points[idx1] = point
	
	def animate(frame):
		'''animation function for an arbitrary number of particles'''
		returns = []
		
		for idx2, point in enumerate(points):#steps through list of particles
			point.set_xdata(Xs[idx2][frame:frame+1])
			point.set_ydata(Ys[idx2][frame:frame+1])
			point.set_3d_properties(Zs[idx2][frame:frame+1])  # update the data for each
			
			if shell_arr[idx2] == 1000: #progenitor
				point.set_color('grey')
				point.set_markersize(2)
			else:
				point.set_color(colors[shell_arr[idx2]])
				point.set_markersize(1)
			
			point.set_marker('o')
			
			returns.append(point)
		
		#print(returns)
		return returns
	
	ani = animation.FuncAnimation(fig, animate, frames = np.arange(1, Xs.shape[1]), interval = 1, blit=True, repeat = True)
	plt.show()
	
	if save:
		ani.save(save_file)

def staticdistribution(foi_list = 'foi.csv', save = False):
#TODO: Change to work with .npz
	
	files = np.loadtxt('foi.csv', dtype='str', delimiter = ',')
	
	thetas_flat = np.zeros(files.size)
	phis_flat = np.zeros(files.size)
	shellnums_flat = np.zeros(files.size)

	for idx in tqdm(range(files.size)):
		currdata = np.loadtxt('./data/'+files[idx], delimiter = ',')
		
		shellnums_flat[idx] = currdata[0, -3]
		thetas_flat[idx] = currdata[0, -2]
		phis_flat[idx] = currdata[0, -1]
	
	num_shells = np.unique(shellnums_flat).size
# 	print(num_shells)
# 	print(thetas_flat.size)

	thetas = np.zeros([num_shells, thetas_flat.size])
	phis = np.zeros([num_shells, phis_flat.size])

	shellnums = np.ones_like(thetas, dtype = int)

	for idx in np.unique(shellnums_flat):
		idx = int(idx)
		thetas[idx, 0:np.count_nonzero(shellnums_flat == idx)] = thetas_flat[shellnums_flat == idx]
		phis[idx, 0:np.count_nonzero(shellnums_flat == idx)] = phis_flat[shellnums_flat == idx]
			
		shellnums[idx] = shellnums[idx]*idx
	
	thetas[thetas == 0.] = np.nan
	phis[phis == 0.] = np.nan
	
	#print(np.unique(shellnums))
	
	draw_spherical(shellnums+1, thetas, phis, colors = ['orange', 'green', 'red', 'purple', 'blue'], plt_labels = np.unique(shellnums), save = save)

#depreciated

# def animate_los(directory = './data/', save = False):
# 	'''sets up, runs, and saves the 2d animation from earth's line of sight
# 	
# 	Args:
# 		directory: string, path to the directory with data to be read in
# 		save: boolean, if True, saves animation to .gif in pwd
# 		
# 	Returns:none
# 	'''
# 	writer = animation.PillowWriter()
# 	
# 	files = os.listdir(directory)
# 	
# 	Y = np.array(len(files), dtype = object)
# 	Z = np.array(len(files), dtype = object)
# 	print(Y)
# 	
# 	#print(len(files))
# 	
# 	for idx in tqdm(range(len(files))):
# 		file = files[idx]
# 		if 'orbit' in file:
# 			#print('FUCK')
# 			continue
# 			
# 		df = read_kep3d_pd(directory, file)
# 		
# 		try:
# 			shell_num = int(df['shellnum'][0])
# 		except:
# 			print(file)
# 			return
# 		
# 		Y[idx] = df['Ys']*u.au
# 		Z[idx] = df['Zs']*u.au
# 		#X = np.asarray(df['Xs']) * u.au
# 		
# 	
# 	print(Y.shape)
# 	return
# 	
# 	fig, ax = plt.subplots()#set up figure
# 	
# 	ax.axis('equal')
# 	offset = .1*u.AU
# 	ax.set_xlim(((Y.min()-offset).value,(Y.max()+offset).value)) #axis shape
# 	ax.set_ylim(((Z.min()-offset).value, (Z.max()+offset).value))
# 	
# 	ax.set_xlabel('Y')
# 	ax.set_ylabel('Z')
# 	ax.set_title('Alderaan System from Line of Sight')
# 	
# 	##"background" imagery
# 	circle1 = plt.Circle((0, 0), 0.1, color='grey')
# 	ax.add_patch(circle1) #primary
# 	
# # 	for i in range(Xs.shape[0]):
# # 		ax.plot(Xs[i],Ys[i], linestyle='--', alpha=.5) #if you want a reference orbit, that can be plotted
# # 	ax.set_prop_cycle(None)
# 	
# 	points = np.ndarray(Y.shape[0], dtype=object) #arbitrary length based on # of filenames
# 	
# # 	for idx1 in range(Y.shape[0]):
# # 		point, = ax.plot(Y[idx1][0], Z[idx1][0]) #sets up 2Dline objects
# # 		points[idx1] = point
# 		
# 	#print(points)
# 
# 	
# 	def animate(frame):
# 		'''animation function for an arbitrary number of particles'''
# 		returns = []
# 		
# 		for idx2, point in enumerate(points):#steps through list of particles
# 			point.set_xdata(Y[idx2][frame:frame+1])
# 			point.set_ydata(X[idx2][frame:frame+1])  # update the data for each
# 		
# 			point.set_marker('o')
# 			
# 			returns.append(point)
# 		
# 		#print(returns)
# 		return returns
# 	
# 	ani = animation.FuncAnimation(fig, animate, frames = np.arange(1, Y.shape[1]), interval=30, blit=True, repeat = True)
# 	plt.show()
# 	
# 	
# 	if save:
# 		ani.save('./2Danimation.gif')

