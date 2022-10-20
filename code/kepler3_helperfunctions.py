import numpy as np
import astropy.units as u
import astropy.constants as c
from os import path
import os

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


##Data handling:
def save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, overwrite = False, filepath = './data/', filename = 'data.csv'):
	'''saves the outputs of kep3d into a basic csv with one line header, added by Catherine Oct 12, 2022
	
	Args:
		X: same as kep3d return
		Y: same as kep3d return
		Xs: same as kep3d return
		Ys: same as kep3d return
		Zs: same as kep3d return
		Xv: same as kep3d return
		Yv: same as kep3d return
		Zv: same as kep3d return
		overwrite: boolean prevents previous data from being overwritten if False, optional
		filepath: string, path to where data should be saved, optional, defaults to pwd
		filename: string, file name with .csv extension, optional, defaults to "data"
		
	Returns: none
	'''
	fullname = filepath+filename
# 	print(fullname)
# 	print(path.exists(fullname))
	
	if overwrite:
		data = np.transpose([X, Y, Xs, Ys, Zs, Xv, Yv, Zv])
		np.savetxt(fullname, data, delimiter=',',header='X, Y, Xs, Ys, Zs, Xv, Yv, Zv')
		
	else:
		if path.exists(fullname):
			print('A file at '+fullname+' already exists, new data NOT saved')
		else:
			data = np.transpose([X, Y, Xs, Ys, Zs, Xv, Yv, Zv])
			np.savetxt(fullname, data, delimiter=',',header='X, Y, Xs, Ys, Zs, Xv, Yv, Zv')
	
	return

def read_kep3d(filepath = './data/', filename = 'data.csv'):
	'''reads in the kep3d data from a csv as saved by the function save_kep3d
	
	Args:
		filepath: string, path to where data should be saved, optional, defaults to pwd
		filename: string, file name with .csv extension, optional, defaults to "data"
		
	Returns:
		X: same as kep3d return
		Y: same as kep3d return
		Xs: same as kep3d return
		Ys: same as kep3d return
		Zs: same as kep3d return
		Xv: same as kep3d return
		Yv: same as kep3d return
		Zv: same as kep3d return
	'''
	
	fullname = filepath+filename
	#print(fullname)
	
	X, Y, Xs, Ys, Zs, Xv, Yv, Zv = np.loadtxt(fullname, comments = '#', delimiter = ',', unpack = True)
	
	return X, Y, Xs, Ys, Zs, Xv, Yv, Zv


#Plotting and animation:

def animate_2d(directory = './data/', save = False):
	'''sets up, runs, and saves the 2d animation in the sky projection plane
	
	Args:
		directory: string, path to the directory with data to be read in
		save: boolean, if True, saves animation to .gif in pwd
		
	Returns:none
	'''
	writer = animation.PillowWriter()

	Xs = [] #setting up orbit data, using regular lists because they append more nicely for our purposes
	Ys = []
	
	for filename in os.listdir(directory):
		#print (filename)
		if not filename.startswith('.'):
			Xs_temp, Ys_temp = read_kep3d(directory, filename)[2:4] #reads in relevant columns
			
			Xs.append(Xs_temp)#adds to list
			Ys.append(Ys_temp)
	
	Xs = np.asarray(Xs)#convert regular lists to np arrays
	Ys = np.asarray(Ys)
	
	print(Ys.shape)
	
	fig, ax = plt.subplots()#set up figure
	
	ax.axis('equal')
	ax.set_xlim((-0.05,0.075)) #axis shape
	ax.set_ylim((-0.025, 0.20))
	
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_title('Title')
	
	##"background" imagery
	ax.scatter(0., 0., marker='o', color='grey') #primary
	ax.plot(Xs[0],Ys[0], linestyle='--', color = 'darkgrey') #if you want a reference orbit, that can be plotted
	
	points = np.ndarray(Xs.shape[0], dtype=object) #arbitrary length based on # of filenames
	
	for idx1 in range(Xs.shape[0]):
		point, = ax.plot(Xs[idx1][0], Ys[idx1][0]) #sets up 2Dline objects
		points[idx1] = point
		
	#print(points)

	
	def animate(frame):
		'''animation function for an arbitrary number of particles'''
		returns = []
		
		for idx2, point in enumerate(points):#steps through list of particles
			point.set_xdata(Xs[idx2][frame:frame+1])
			point.set_ydata(Ys[idx2][frame:frame+1])  # update the data for each
		
			point.set_marker('o')
			
			returns.append(point)
		
		#print(returns)
		return returns
	
	ani = animation.FuncAnimation(fig, animate, frames = np.arange(1, Xs.shape[1]), interval=30, blit=True, repeat = True)
	plt.show()
	
	
	if save:
		ani.save('./2Danimation.gif')
		
def animate_3d(directory = './data/', save = False):
	'''sets up, runs, and saves a 3d animation
	
	Args:
		filepath: string, path to the directory with data to be read in
		save: boolean, if True, saves animation to .gif in pwd
		
	Returns:none
	'''
	writer = animation.PillowWriter()

	Xs = [] #setting up orbit data, using regular lists because they append more nicely for our purposes
	Ys = []
	Zs = []
	for filename in os.listdir(directory):
		if not filename.startswith('.'):
	
			Xs_temp, Ys_temp, Zs_temp = read_kep3d(directory, filename)[2:5] #reads in relevant columns
		
			Xs.append(Xs_temp)#adds to list
			Ys.append(Ys_temp)
			Zs.append(Zs_temp)
	
	Xs = np.asarray(Xs)#convert regular lists to np arrays
	Ys = np.asarray(Ys)
	Zs = np.asarray(Zs)
	
	fig, ax = plt.subplots(figsize=(8,7),subplot_kw={'projection': '3d'})#set up figure
	
	h = 0.2
	#ax.axis('equal')
	ax.set_xlim(-h,h)
	ax.set_ylim(-h,h)
	ax.set_zlim(h,-h)
	
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	ax.set_title('Title')
	
	ax.scatter(0., 0., 0., marker='o', color='grey') #primary
	ax.plot(Xs[0],Ys[0],Zs[0], linestyle='--', color = 'darkgrey') #if you want a reference orbit, that can be plotted
	
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
		
			point.set_marker('o')
			
			returns.append(point)
		
		#print(returns)
		return returns
	
	ani = animation.FuncAnimation(fig, animate, frames = np.arange(1, Xs.shape[1]), interval=30, blit=True, repeat = True)
	plt.show()
	
	if save:
		ani.save('./3Danimation.gif')
	