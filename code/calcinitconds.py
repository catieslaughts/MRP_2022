import numpy as np
import astropy.units as u

def findstartpos(vorb:u.km/u.s, vkick:u.km/u.s, M:u.kg, m:u.kg, deltat:u.yr):
	''' Determines the location of the explosion given some input parameters
		
		vorb: orbital velocity of the secondary pre-explosion
		vkick: maximum value of velocity kick
		M: Primary mass
		m: secondary mass
		deltat: time between explosion and first transit
	'''
	
	return