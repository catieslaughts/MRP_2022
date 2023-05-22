import numpy as np
import pandas as pd
import os
import os.path as path
import glob
import astropy.units as u

def save_kep3d_npz(Xs, Ys, Zs, shellnum = np.nan, theta = np.nan, phi=np.nan, overwrite = False, filepath = './data/', filename = 'data'):
	
	fullname = filepath+filename
	
	Xs = Xs.to(u.au).value
	Ys = Ys.to(u.au).value
	Zs = Zs.to(u.au).value
	
	if overwrite:
		np.savez(fullname, Xs=Xs, Ys=Ys, Zs=Zs, shellnum=shellnum, theta=theta, phi=phi)
		
	else:
		if path.exists(fullname):
			print('A file at '+fullname+' already exists, new data NOT saved')
		else:
			np.savez(fullname, Xs=Xs, Ys=Ys, Zs=Zs, shellnum=shellnum, theta=theta, phi=phi)
	
	return

def read_kep3d_npz(filepath = './data/', filename = 'data'):
	
	fullname = filepath+filename
	
	npzfile = np.load(fullname)
	
	Xs = npzfile['Xs']
	Ys = npzfile['Ys']
	Zs = npzfile['Zs']
	
	Xs = Xs *u.au
	Ys = Ys *u.au
	Zs = Zs *u.au
	
	shellnum = int(npzfile['shellnum'])
	
	theta = npzfile['theta']
	phi = npzfile['phi']
	
	theta = theta *u.rad
	phi = phi *u.rad
	
	npzfile.close()
	
	return Xs, Ys, Zs, shellnum, theta, phi

#BE VERY CAREFUL CALLING THIS FUNCTION, FOR GOD'S SAKE
def empty_data_directory(directory):
	fileList = glob.glob(directory+'*')
	for filePath in fileList:
		os.remove(filePath)

#BE VERY CAREFUL CALLING THIS FUNCTION, FOR GOD'S SAKE
def empty_extra_data(directory):
	fileList = glob.glob(directory+'cloud_extra*')
	print(str(len(fileList))+' files removed from '+directory)
	for filePath in fileList:
		os.remove(filePath)
		
@u.quantity_input
def write_param_file(kick_vel:u.m/u.s = 1.*u.km/u.s, v_i:u.km/u.s = 6.*u.km/u.s, 
		t_firstpass:u.day = 750 * u.day, t_start:u.year = 2000*u.year, 
		delta_t:u.year = 50*u.year, e_i = 0., omega_i:u.rad = 0*u.rad, 
		anode_i:u.rad = 0*u.rad, M:u.kg = 1.0*u.Msun, m:u.kg = 1.0*u.Mearth, 
		num_per_shell = 400, n_steps = 500, n_pre_steps = 0, n_shells = 3, R_star:u.au = 1*u.Rsun, 
		paramfile = 'params.csv'):
		
	paramdata = np.asarray([kick_vel.to(u.km/u.s).value, v_i.to(u.km/u.s).value, t_firstpass.to(u.year).value, t_start.to(u.year).value, delta_t.to(u.year).value, e_i, omega_i.to(u.rad).value, anode_i.to(u.rad).value, M.to(u.kg).value, m.to(u.kg).value, num_per_shell, n_steps, n_pre_steps, n_shells, R_star.to(u.au).value])
	paramheader = 'kick_vel,v_i,t_firstpass,t_start,delta_t,e_i,omega_i,anode_i,M,m,num_per_shell,n_steps,n_pre_steps,n_shells,R_star\n(u.km/u.s),(u.km/u.s),(u.year),(u.year),(u.year), ,(u.rad),(u.rad),(u.kg),(u.kg), , , , ,(u.au)'
	
	np.savetxt(paramfile, [paramdata], delimiter=',', header = paramheader)
	
def read_param_file(paramfile = 'params.csv'):
	kick_vel,v_i,t_firstpass,t_start,delta_t,e_i,omega_i,anode_i,M,m,num_per_shell,n_steps,n_pre_steps,n_shells,R_star = np.loadtxt(paramfile, delimiter = ',')
	
	kick_vel, v_i, t_firstpass, omega_i, anode_i, M, m, t_start, delta_t, R_star = kick_vel*(u.km/u.s), v_i*(u.km/u.s), t_firstpass*u.year, omega_i*u.rad, anode_i*u.rad, M*u.kg, m*u.kg, t_start*u.year, delta_t*u.year, R_star*u.au
	
	n_shells, num_per_shell, n_steps, n_pre_steps = int(n_shells), int(num_per_shell), int(n_steps), int(n_pre_steps)
	
	return kick_vel,v_i,t_firstpass,t_start,delta_t,e_i,omega_i,anode_i,M,m,num_per_shell,n_steps,n_pre_steps,n_shells,R_star


# Depreciated:

# def save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, shellnums = np.asarray([], dtype=int), overwrite = False, filepath = './data/', filename = 'data.csv', pre_steps = 0, post_steps = 0):
# 	'''saves the outputs of kep3d into a basic csv with one line header, added by Catherine Oct 12, 2022
# 	
# 	Args:
# 		X: same as kep3d return
# 		Y: same as kep3d return
# 		Xs: same as kep3d return
# 		Ys: same as kep3d return
# 		Zs: same as kep3d return
# 		Xv: same as kep3d return
# 		Yv: same as kep3d return
# 		Zv: same as kep3d return
# 		overwrite: boolean prevents previous data from being overwritten if False, optional
# 		filepath: string, path to where data should be saved, optional, defaults to pwd
# 		filename: string, file name with .csv extension, optional, defaults to "data"
# 		
# 	Returns: none
# 	'''
# 	fullname = filepath+filename
# 	
# 	if shellnums.size == 0:
# 		data = np.transpose([X.to(u.au).value, Y.to(u.au).value, Xs.to(u.au).value, Ys.to(u.au).value, Zs.to(u.au).value, Xv.to(u.km/u.s).value, Yv.to(u.km/u.s).value, Zv.to(u.km/u.s).value])
# 		header = 'X,Y,Xs,Ys,Zs,Xv,Yv,Zv'
# 	else:
# 		data = np.transpose([X.to(u.au).value, Y.to(u.au).value, Xs.to(u.au).value, Ys.to(u.au).value, Zs.to(u.au).value, Xv.to(u.km/u.s).value, Yv.to(u.km/u.s).value, Zv.to(u.km/u.s).value, shellnums.value])
# 		header = 'X,Y,Xs,Ys,Zs,Xv,Yv,Zv,shellnum,theta,phi'
# 	
# 	if overwrite:
# 		if pre_steps > 0:
# 			empty_data = np.full((pre_steps,data.shape[1]), np.nan)
# 			np.savetxt(fullname, empty_data, delimiter=',',header=header)
# 			with open(fullname, 'ab') as f:
# 				np.savetxt(f, data, delimiter=',')
# 		else:
# 			np.savetxt(fullname, data, delimiter=',',header=header)
# 		
# 		if post_steps > 0:
# 			empty_data = np.full((post_steps,data.shape[1]), np.nan)
# 			with open(fullname, 'ab') as f:
# 				np.savetxt(f, empty_data, delimiter=',')
# 		
# 	else:
# 		if path.exists(fullname):
# 			print('A file at '+fullname+' already exists, new data NOT saved')
# 		else:
# 			if pre_steps > 0:
# 				empty_data = np.full((pre_steps,3), np.nan)
# 				np.savetxt(fullname, empty_data, delimiter=',',header='X, Y, Xs, Ys, Zs, Xv, Yv, Zv')
# 				with open(fullname, 'ab') as f:
# 					np.savetxt(f, data, delimiter=',')
# 			else:
# 				np.savetxt(fullname, data, delimiter=',',header='X, Y, Xs, Ys, Zs, Xv, Yv, Zv')
# 			
# 			if post_steps > 0:
# 				empty_data = np.full((post_steps,3), np.nan)
# 				with open(fullname, 'ab') as f:
# 					np.savetxt(f, empty_data, delimiter=',')
# 	
# 	return
 
# def read_kep3d(filepath = './data/', filename = 'data.csv'):
# 	'''reads in the kep3d data from a csv as saved by the function save_kep3d
# 	
# 	Args:
# 		filepath: string, path to where data should be saved, optional, defaults to pwd
# 		filename: string, file name with .csv extension, optional, defaults to "data"
# 		
# 	Returns:
# 		X: same as kep3d return
# 		Y: same as kep3d return
# 		Xs: same as kep3d return
# 		Ys: same as kep3d return
# 		Zs: same as kep3d return
# 		Xv: same as kep3d return
# 		Yv: same as kep3d return
# 		Zv: same as kep3d return
# 	'''
# 	
# 	fullname = filepath+filename
# 	#print(fullname)
# 	
# 	X, Y, Xs, Ys, Zs, Xv, Yv, Zv = np.loadtxt(fullname, delimiter = ',')
# 	
# 	return X, Y, Xs, Ys, Zs, Xv, Yv, Zv

# def save_kep3d_short(Xs, Ys, Zs, overwrite = False, filepath = './data/', filename = 'data.csv', pre_steps = 0, post_steps = 0):
# 	'''saves the outputs of kep3d into a basic csv with one line header, added by Catherine Oct 12, 2022
# 	
# 	Args:
# 		Xs: same as kep3d return
# 		Ys: same as kep3d return
# 		Zs: same as kep3d return
# 		overwrite: boolean prevents previous data from being overwritten if False, optional
# 		filepath: string, path to where data should be saved, optional, defaults to pwd
# 		filename: string, file name with .csv extension, optional, defaults to "data"
# 		pre_steps: empty steps to put before data, for visualization purposes
# 		post_steps: empty steps to put after data, for visualization purposes
# 		
# 	Returns: none
# 	'''
# 	
# 	fullname = filepath+filename
# # 	print(fullname)
# # 	print(path.exists(fullname))
# 	
# 	if overwrite:
# 		data = np.transpose([Xs.to(u.au).value, Ys.to(u.au).value, Zs.to(u.au).value])
# 		if pre_steps > 0:
# 			empty_data = np.full((pre_steps,3), np.nan)
# 			np.savetxt(fullname, empty_data, delimiter=',',header=' Xs,Ys,Zs')
# 			with open(fullname, 'ab') as f:
# 				np.savetxt(f, data, delimiter=',')
# 		else:
# 			np.savetxt(fullname, data, delimiter=',',header=' Xs,Ys,Zs')
# 		
# 		if post_steps > 0:
# 			empty_data = np.full((post_steps,3), np.nan)
# 			with open(fullname, 'ab') as f:
# 				np.savetxt(f, empty_data, delimiter=',')
# 		
# 	else:
# 		if path.exists(fullname):
# 			print('A file at '+fullname+' already exists, new data NOT saved')
# 		else:
# 			data = np.transpose([Xs.to(u.au).value, Ys.to(u.au).value, Zs.to(u.au).value])
# 			if pre_steps > 0:
# 				empty_data = np.full((pre_steps,3), np.nan)
# 				np.savetxt(fullname, empty_data, delimiter=',',header='Xs,Ys,Zs')
# 				with open(fullname, 'ab') as f:
# 					np.savetxt(f, data, delimiter=',')
# 			else:
# 				np.savetxt(fullname, data, delimiter=',',header='Xs,Ys,Zs')
# 			
# 			if post_steps > 0:
# 				empty_data = np.full((post_steps,3), np.nan)
# 				with open(fullname, 'ab') as f:
# 					np.savetxt(f, empty_data, delimiter=',')
# 	
# 	return

# def read_kep3d_short(filepath = './data/', filename = 'data.csv'):
# 	'''reads in the kep3d data from a csv as saved by the function save_kep3d
# 	
# 	Args:
# 		filepath: string, path to where data should be saved, optional, defaults to pwd
# 		filename: string, file name with .csv extension, optional, defaults to "data"
# 		
# 	Returns:
# 		Xs: same as kep3d return
# 		Ys: same as kep3d return
# 		Zs: same as kep3d return
# 	'''
# 	
# 	fullname = filepath+filename
# 	#print(fullname)
# 	
# 	Xs, Ys, Zs = np.loadtxt(fullname, comments = '#', delimiter = ',', unpack = True)
# 	
# 	return Xs, Ys, Zs

# def read_kep3d_pd(filepath = './data/', filename = 'data.csv'):
# 	'''reads in output of save_kep3d, but in a generalizable pandas dataframe'''
# 	
# 	fullname = filepath+filename
# 	
# 	df = pd.read_csv(fullname,skipinitialspace = True)
# 	
# 	firstcol = df.columns[0][3:]
# 	df = df.rename(columns={ df.columns[0]: firstcol})
# 	
# 	#print()
# 	
# 	return df
