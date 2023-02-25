import numpy as np
import pandas as pd
import os
import os.path as path
import glob
import astropy.units as u

def save_kep3d(X, Y, Xs, Ys, Zs, Xv, Yv, Zv, shellnums = np.asarray([], dtype=int), overwrite = False, filepath = './data/', filename = 'data.csv', pre_steps = 0, post_steps = 0):
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
	
	#print(shellnums)
	
	if shellnums.size == 0:
		data = np.transpose([X.to(u.au).value, Y.to(u.au).value, Xs.to(u.au).value, Ys.to(u.au).value, Zs.to(u.au).value, Xv.to(u.km/u.s).value, Yv.to(u.km/u.s).value, Zv.to(u.km/u.s).value])
		header = 'X,Y,Xs,Ys,Zs,Xv,Yv,Zv'
	else:
		data = np.transpose([X.to(u.au).value, Y.to(u.au).value, Xs.to(u.au).value, Ys.to(u.au).value, Zs.to(u.au).value, Xv.to(u.km/u.s).value, Yv.to(u.km/u.s).value, Zv.to(u.km/u.s).value, shellnums.value])
		header = 'X,Y,Xs,Ys,Zs,Xv,Yv,Zv,shellnum'
	
	if overwrite:
		if pre_steps > 0:
			empty_data = np.full((pre_steps,data.shape[1]), np.nan)
			np.savetxt(fullname, empty_data, delimiter=',',header=header)
			with open(fullname, 'ab') as f:
				np.savetxt(f, data, delimiter=',')
		else:
			np.savetxt(fullname, data, delimiter=',',header=header)
		
		if post_steps > 0:
			empty_data = np.full((post_steps,data.shape[1]), np.nan)
			with open(fullname, 'ab') as f:
				np.savetxt(f, empty_data, delimiter=',')
		
	else:
		if path.exists(fullname):
			print('A file at '+fullname+' already exists, new data NOT saved')
		else:
			if pre_steps > 0:
				empty_data = np.full((pre_steps,3), np.nan)
				np.savetxt(fullname, empty_data, delimiter=',',header='X, Y, Xs, Ys, Zs, Xv, Yv, Zv')
				with open(fullname, 'ab') as f:
					np.savetxt(f, data, delimiter=',')
			else:
				np.savetxt(fullname, data, delimiter=',',header='X, Y, Xs, Ys, Zs, Xv, Yv, Zv')
			
			if post_steps > 0:
				empty_data = np.full((post_steps,3), np.nan)
				with open(fullname, 'ab') as f:
					np.savetxt(f, empty_data, delimiter=',')
	
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
	
def save_kep3d_short(Xs, Ys, Zs, overwrite = False, filepath = './data/', filename = 'data.csv', pre_steps = 0, post_steps = 0):
	'''saves the outputs of kep3d into a basic csv with one line header, added by Catherine Oct 12, 2022
	
	Args:
		Xs: same as kep3d return
		Ys: same as kep3d return
		Zs: same as kep3d return
		overwrite: boolean prevents previous data from being overwritten if False, optional
		filepath: string, path to where data should be saved, optional, defaults to pwd
		filename: string, file name with .csv extension, optional, defaults to "data"
		pre_steps: empty steps to put before data, for visualization purposes
		post_steps: empty steps to put after data, for visualization purposes
		
	Returns: none
	'''
	
	fullname = filepath+filename
# 	print(fullname)
# 	print(path.exists(fullname))
	
	if overwrite:
		data = np.transpose([Xs.to(u.au).value, Ys.to(u.au).value, Zs.to(u.au).value])
		if pre_steps > 0:
			empty_data = np.full((pre_steps,3), np.nan)
			np.savetxt(fullname, empty_data, delimiter=',',header=' Xs,Ys,Zs')
			with open(fullname, 'ab') as f:
				np.savetxt(f, data, delimiter=',')
		else:
			np.savetxt(fullname, data, delimiter=',',header=' Xs,Ys,Zs')
		
		if post_steps > 0:
			empty_data = np.full((post_steps,3), np.nan)
			with open(fullname, 'ab') as f:
				np.savetxt(f, empty_data, delimiter=',')
		
	else:
		if path.exists(fullname):
			print('A file at '+fullname+' already exists, new data NOT saved')
		else:
			data = np.transpose([Xs.to(u.au).value, Ys.to(u.au).value, Zs.to(u.au).value])
			if pre_steps > 0:
				empty_data = np.full((pre_steps,3), np.nan)
				np.savetxt(fullname, empty_data, delimiter=',',header='Xs,Ys,Zs')
				with open(fullname, 'ab') as f:
					np.savetxt(f, data, delimiter=',')
			else:
				np.savetxt(fullname, data, delimiter=',',header='Xs,Ys,Zs')
			
			if post_steps > 0:
				empty_data = np.full((post_steps,3), np.nan)
				with open(fullname, 'ab') as f:
					np.savetxt(f, empty_data, delimiter=',')
	
	return

def read_kep3d_short(filepath = './data/', filename = 'data.csv'):
	'''reads in the kep3d data from a csv as saved by the function save_kep3d
	
	Args:
		filepath: string, path to where data should be saved, optional, defaults to pwd
		filename: string, file name with .csv extension, optional, defaults to "data"
		
	Returns:
		Xs: same as kep3d return
		Ys: same as kep3d return
		Zs: same as kep3d return
	'''
	
	fullname = filepath+filename
	#print(fullname)
	
	Xs, Ys, Zs = np.loadtxt(fullname, comments = '#', delimiter = ',', unpack = True)
	
	return Xs, Ys, Zs

def read_kep3d_pd(filepath = './data/', filename = 'data.csv'):
	'''reads in output of save_kep3d, but in a generalizable pandas dataframe'''
	
	fullname = filepath+filename
	
	df = pd.read_csv(fullname,skipinitialspace = True)
	
	firstcol = df.columns[0][3:]
	df = df.rename(columns={ df.columns[0]: firstcol})
	
	#print()
	
	return df

#BE VERY CAREFUL CALLING THIS FUNCTION, FOR GOD'S SAKE
def empty_data_directory(directory):
	fileList = glob.glob(directory+'*')
	for filePath in fileList:
		os.remove(filePath)
		
#read_kep3d_pd('./data/', 'cloud0.csv')



