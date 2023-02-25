import numpy as np
from fileio import *
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

#Functions for animating the outputs of kep3d, position and velocity data can be saved from kep3d outputs
#using file write functions in fileio.py. Each steps through all the inidvidual files in a given directory,
#reading each and adding the corresponding point to the animation

#TODO: generalize animate_3d so it can read in either a file from save_kep3d or save_kep3d_short, Pandas would probably be easiest for this

#animate_2d is depreciated
# def animate_2d(directory = './data/', save = False):
# 	'''sets up, runs, and saves the 2d animation in the sky projection plane
# 	
# 	Args:
# 		directory: string, path to the directory with data to be read in
# 		save: boolean, if True, saves animation to .gif in pwd
# 		
# 	Returns:none
# 	'''
# 	writer = animation.PillowWriter()
# 
# 	Xs = [] #setting up orbit data, using regular lists because they append more nicely for our purposes
# 	Ys = []
# 	
# 	files = os.listdir(directory)
# 	files.reverse()
# 	for filename in files:
# 		#print (filename)
# 		if not filename.startswith('.'):
# 			Xs_temp, Ys_temp = read_kep3d(directory, filename)[2:4] #reads in relevant columns
# 			
# 			Xs.append(Xs_temp)#adds to list
# 			Ys.append(Ys_temp)
# 	
# 	Xs = np.asarray(Xs)#convert regular lists to np arrays
# 	Ys = np.asarray(Ys)
# 	
# 	#print(Ys.shape)
# 	
# 	fig, ax = plt.subplots()#set up figure
# 	
# 	ax.axis('equal')
# 	offset = .1
# 	ax.set_xlim((Xs.min()-offset,Xs.max()+offset)) #axis shape
# 	ax.set_ylim((Ys.min()-offset, Ys.max()+offset))
# 	
# 	ax.set_xlabel('X')
# 	ax.set_ylabel('Y')
# 	ax.set_title('Title')
# 	
# 	##"background" imagery
# 	ax.scatter(0., 0., marker='o', color='grey') #primary
# 	
# 	for i in range(Xs.shape[0]):
# 		ax.plot(Xs[i],Ys[i], linestyle='--', alpha=.5) #if you want a reference orbit, that can be plotted
# 	ax.set_prop_cycle(None)
# 	
# 	points = np.ndarray(Xs.shape[0], dtype=object) #arbitrary length based on # of filenames
# 	
# 	for idx1 in range(Xs.shape[0]):
# 		point, = ax.plot(Xs[idx1][0], Ys[idx1][0]) #sets up 2Dline objects
# 		points[idx1] = point
# 		
# 	#print(points)
# 
# 	
# 	def animate(frame):
# 		'''animation function for an arbitrary number of particles'''
# 		returns = []
# 		
# 		for idx2, point in enumerate(points):#steps through list of particles
# 			point.set_xdata(Xs[idx2][frame:frame+1])
# 			point.set_ydata(Ys[idx2][frame:frame+1])  # update the data for each
# 		
# 			point.set_marker('o')
# 			
# 			returns.append(point)
# 		
# 		#print(returns)
# 		return returns
# 	
# 	ani = animation.FuncAnimation(fig, animate, frames = np.arange(1, Xs.shape[1]), interval=30, blit=True, repeat = True)
# 	plt.show()
# 	
# 	
# 	if save:
# 		ani.save('./2Danimation.gif')
		
def animate_3d(directory = './data/', save = False, save_file = './3Danimation.gif', velocitycolor = True):
	'''sets up, runs, and saves a 3d animation
	
	Args:
		filepath: string, path to the directory with data to be read in
		save: boolean, if True, saves animation to .gif in pwd
		
	Returns:none
	'''
	colors = ['orange', 'green', 'red', 'purple'] #outermost shell is colored by the last one
	writer = animation.PillowWriter()

	Xs = [] #setting up orbit data, using regular lists because they append more nicely for our purposes
	Ys = []
	Zs = []
	vel = []
	shell_arr = []
	
	files = os.listdir(directory)
	#files.reverse()
	orbitflag = 0
	for filename in files:
		#for making animations which show the original progenitor and/or orbit
		if not filename.startswith('.'):
			df = read_kep3d_pd(directory, filename)
			if 'orbit' in filename: # a file showing a single orbit of the progenitor
				orbitflag = 1
				orbitX, orbitY, orbitZ = np.asarray(df['Xs']), np.asarray(df['Ys']), np.asarray(df['Zs'])
			else:
				Xs_temp, Ys_temp, Zs_temp = np.asarray(df['Xs']), np.asarray(df['Ys']), np.asarray(df['Zs'])
				Xs.append(Xs_temp)#adds to list
				Ys.append(Ys_temp)
				Zs.append(Zs_temp)
				
				if velocitycolor:
					Xv, Yv, Zv = np.asarray(df['Xs']), np.asarray(df['Ys']), np.asarray(df['Zs'])
					vel_temp = np.sqrt(Xv[0]**2 + Yv[0]**2 + Zv[0]**2)
					
					vel.append(vel_temp)
				if 'shellnum' in df.columns:
					shell_arr.append(np.asarray(df['shellnum'][-2:-1])[0])
				else:
					shell_arr.append(1000)
	
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

animate_3d(directory = './test_subset/', save=True)

