import numpy as np
import astropy.units as u
import astropy.constants as c
from os import path
import os

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

from kepler3 import *

##create point spheres:
@u.quantity_input
def random_sphere(n_layers, n_particles, max_vel:u.m/u.s):
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

@u.quantity_input
def fibonacci_sphere(n_layers, n_particles, max_vel:u.m/u.s):
	'''makes use of the fibonacci sphere algorithm to create several circular layers of points. 
	The radius output by this function will be used as our velocities in our simulation
	Sources:
		http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
		https://arxiv.org/pdf/0912.4540.pdf
		
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
	theta = np.zeros([n_layers, n_particles], dtype = float)
	phi = np.zeros([n_layers, n_particles], dtype = float)
	
	gr = (1 + np.sqrt(5))/2
	parts = np.arange(0,n_particles)
	
	vels = np.linspace(max_vel, 0, n_layers, endpoint=False)[::-1] #avoids a layer with 0 velocity
	
	for idx, vel in enumerate(vels):
		
		r[idx] = np.ones(n_particles)*vel
		phi[idx] = 2*np.pi * parts / gr
		theta[idx] = np.arccos(1 - 2*(parts+0.5)/n_particles)
		#theta and phi are flipped from the given sources because the conventions are flipped in math and physics
	
	#units
	r=r*max_vel.unit
	theta = theta*u.rad
	phi = phi*u.rad
	
	return r, theta, phi
	
def latlong_sphere(n_layers, n_particles, max_vel:u.m/u.s):
	'''ummmm'''
	
	r = np.zeros([n_layers, n_particles], dtype = float)
	theta = np.zeros([n_layers, n_particles], dtype = float)
	phi = np.zeros([n_layers, n_particles], dtype = float)
	pass

def sph2cart(r, theta, phi):
	'''eqs 4.4.1-3 in 
	https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Book%3A_Electromagnetics_I_(Ellingson)/04%3A_Vector_Analysis/4.04%3A_Spherical_Coordinates'''
	
	x = r * np.cos(phi) * np.sin(theta)
	y = r * np.sin(phi) * np.sin(theta)
	z = r * np.cos(theta)
	
	return x, y, z
	
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

##New Kep3d:
def kep3d_trueanom(epoch:u.year, P:u.year, true_anom, a, e, inc:u.deg, omega:u.deg, anode:u.deg):
    """
    Calculate the position and velocity of an orbiting body

    Given the Kepler elements for the secondary about the primary
    and in the coordinate frame where the primary is at the origin
    
    The same as the kep3d function in kepler3.py, but now takes in true anomaly instead of tperi

    Args:
        epoch (np.array):  epochs to evaluate (u.time)
        P (np.array): orbital period (u.time) 
        true_anom (float): the true anomaly of the orbit
        a (float): semi-major axis of the orbit
        e (float): eccentricity of the orbit
        inc (float): inclination of the orbit  (u.angle)
        omega (float): longitude of periastron (u.angle)
        anode (float): PA of the ascending node (u.angle)

    Returns:
       X,Y, Xs,Ys,Zs, Xsv,Ysv,Zsv

    Output frame has X,Y in computer plotting coordinates
    i.e. X is to the right, increasing (due West)

    Primary body is fixed at the origin.

    X,Y (float): 2D coordinates of in plane orbit with periapse
                 towards the +ve X axis.

    Xs,Ys,Zs (float): The 3D coordinates of the secondary body
        in the Position/velocity coords frame.

    Xsv, Ysv, Zsv (float): The 3D velocity of the secondary body
        in the Position/velocity coords frame.

    Coordinate frames are shown below.

    The 3D axes are NOT the usual right-handed coordinate frame. The
    Observer is located far away on the NEGATIVE Z axis. This is done
    so that +ve Zsv gives positive velocities consistent with most
    astronomers idea of redshift being positive velocity values.


    Sky coords         Computer coords   Position/velocity coords

      North                   Y                +Y    +Z
        ^                     ^                 ^   ^
        |                     |                 |  /
        |                     |                 | /
        |                     |                 |/
        +-------> West        +-------> X       +-------> +X
                                               /
                                              /
                                             /
                                           -Z

    +Y is North, +X is West and +Z is away from the Earth
    so that velocity away from the Earth is positive

    NOTE: Right Ascension is INCREASING to the left, but the
    (X,Y,Z) outputs have RA increasing to the right, as seen
    in the Computer coords. This is done to make plotting easier
    and to remind the user that for the sky coords they need
    to plot (-X,Y) and then relabel the RA/X axis to reflect 
    what is seen in the sky.

    Taken from ``Lecture Notes on Basic Celestial Mechanics''
    by Sergei A. Klioner (2011) page 22
    http://astro.geo.tu-dresden.de/~klioner/celmech.pdf

    Note that mean motion is defined on p.17 and that
    (P^2/a^3) = 4pi^2/kappa^2 and that the document is missing
    the ^2 on the pi.

    Written:
        Matthew Kenworthy, 2017

    """

    # mean motion n
    n = 2 * np.pi / P

    # calc eccentric anomaly E
    v = true_anom
    E = np.arctan2(np.sqrt(1 - e**2) * np.sin(v), e + np.cos(v))
    

    # calculate position and velocity in the orbital plane
    cE = np.cos(E)
    sE = np.sin(E)
    surde = np.sqrt(1 - (e*e))

    X  = a * (cE - e)
    Y  = a * surde * sE

    Xv = -(a * n * sE) / (1 - e * cE)
    Yv =  (a * n * surde * cE) / (1 - e * cE)

    # calculate Euler rotation matrix to get from orbital plane
    # to the sky plane

    mat = euler(anode, omega, inc)

    # rotate the coordinates from the orbital plane
    # to the sky projected coordinates

    # TODO we lose the dimensionality of X and Y and it
    # needs to be put back artificially
    # problem with building the np.array below and putting the Quantity through
    # the np.dot() routine
    
    (Xe, Ye, Ze)    = np.dot(mat, np.array([X.value,Y.value,0],dtype=object))
    #blog = np.array([X,Y,np.zeros(X.size)]) * X.unit
    #(Xe, Ye, Ze)    = np.dot(mat, blog)
    (Xev, Yev, Zev) = np.dot(mat, np.array([Xv.value,Yv.value,0],dtype=object))
    
    ##CATHERINE added dtype=object to these two lines to avoid the depreciation warning

    Xs = -Ye * X.unit
    Ys =  Xe * X.unit
    Zs =  Ze * X.unit
    Xsv = -Yev * Xv.unit
    Ysv =  Xev * Xv.unit
    Zsv =  Zev * Xv.unit

    # To transform from the Kep3D code coordinate system to the
    # celestial sphere, where:
    # Y is North, X is West and Z is away from the Earth
    # we have
    #   X,Y,Z,Xdot,Ydot,Zdot = -Ys, Xs, Zs, -Yv, Xv, Zv

    return(X,Y,Xs,Ys,Zs,Xsv,Ysv,Zsv)

def tperi_to_Tanom(tperi, epoch, P, e, derror=1e-6):
	# at epoch = tperi, mean anomoly M is 0
    # 
    # Y = time since epoch periastron
    Y = epoch - tperi

    # mean anomaly M varies smoothly from 0 to 2pi every orbital period
    # convert Y to angle in radians between 0 and 2PI
    Mt = Y / P
    M  = (2 * np.pi * (Mt - np.floor(Mt)))*u.radian
    
    (E,v) = kepler_solve(e, M, derror)
    
    #print(v)
    
    return E, v.value

##Data handling:
def save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, overwrite = False, filepath = './earthsuntest/', filename = 'data.csv'):
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
		data = np.transpose([X.value, Y.value, Xs.value, Ys.value, Zs.value, Xv.value, Yv.value, Zv.value])
		np.savetxt(fullname, data, delimiter=',',header='X, Y, Xs, Ys, Zs, Xv, Yv, Zv')
		
	else:
		if path.exists(fullname):
			print('A file at '+fullname+' already exists, new data NOT saved')
		else:
			data = np.transpose([X.value, Y.value, Xs.value, Ys.value, Zs.value, Xv.value, Yv.value, Zv.value])
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
	
	files = os.listdir(directory)
	files.reverse()
	for filename in files:
		#print (filename)
		if not filename.startswith('.'):
			Xs_temp, Ys_temp = read_kep3d(directory, filename)[2:4] #reads in relevant columns
			
			Xs.append(Xs_temp)#adds to list
			Ys.append(Ys_temp)
	
	Xs = np.asarray(Xs)#convert regular lists to np arrays
	Ys = np.asarray(Ys)
	
	#print(Ys.shape)
	
	fig, ax = plt.subplots()#set up figure
	
	ax.axis('equal')
	offset = .1
	ax.set_xlim((Xs.min()-offset,Xs.max()+offset)) #axis shape
	ax.set_ylim((Ys.min()-offset, Ys.max()+offset))
	
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_title('Title')
	
	##"background" imagery
	ax.scatter(0., 0., marker='o', color='grey') #primary
	
	for i in range(Xs.shape[0]):
		ax.plot(Xs[i],Ys[i], linestyle='--', alpha=.5) #if you want a reference orbit, that can be plotted
	ax.set_prop_cycle(None)
	
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
	
	files = os.listdir(directory)
	files.reverse()
	for filename in files:
		if not filename.startswith('.'):
	
			Xs_temp, Ys_temp, Zs_temp = read_kep3d(directory, filename)[2:5] #reads in relevant columns
		
			Xs.append(Xs_temp)#adds to list
			Ys.append(Ys_temp)
			Zs.append(Zs_temp)
	
	Xs = np.asarray(Xs)#convert regular lists to np arrays
	Ys = np.asarray(Ys)
	Zs = np.asarray(Zs)
	
	fig, ax = plt.subplots(figsize=(8,7),subplot_kw={'projection': '3d'})#set up figure
	
	offset = .1
	ax.set_xlim((Xs.min()-offset,Xs.max()+offset)) #axis shape
	ax.set_ylim((Ys.min()-offset, Ys.max()+offset))
	ax.set_zlim((Zs.min()-offset, Zs.max()+offset))
	
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	ax.set_title('Title')
	
	ax.scatter(0., 0., 0., marker='o', color='grey') #primary
	for i in range(Xs.shape[0]):
		ax.plot(Xs[i],Ys[i], Zs[i], linestyle='--', alpha=.5) #if you want a reference orbit, that can be plotted
	ax.set_prop_cycle(None) #resets color cycler so orbits match up with particles
	
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
	