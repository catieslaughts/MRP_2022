import numpy as np
import astropy.units as u
import astropy.constants as c

def findstartpos(vorb:u.km/u.s, vkick:u.km/u.s, M:u.kg, m:u.kg, t_start:u.yr):
	''' Determines the location of the explosion given some input parameters
		
		vorb: orbital velocity of the secondary pre-explosion
		vkick: maximum value of velocity kick
		M: Primary mass
		m: secondary mass
		deltat: time between explosion and first transit
	'''
	#before explosion:
	a_i = (c.G*(M+m)/(vorb**2)).to(u.au)
	
	P_i=np.sqrt((a_i**3)*4*np.pi*np.pi/(c.G*(M+m))).to(u.yr)
	
	
	#after explosion, 1st possilbe particle to cross
	v_f = vorb+vkick
	a_f = (c.G*(M+m)/(v_f**2)).to(u.au)
	P_f=np.sqrt((a_f**3)*4*np.pi*np.pi/(c.G*(M+m))).to(u.yr)
	
	l_ex = (v_f*t_start).decompose()
	
	theta_ex = (l_ex/a_f).decompose()*u.rad
	
	#calculate fraction of circle:
	circ_frac = theta_ex/(2.*np.pi*u.rad)
	
	#convert to time before crossing explosion occurs:
	deltat = (P_f * circ_frac).to(u.day)
	
	
	print(deltat)
	
	return