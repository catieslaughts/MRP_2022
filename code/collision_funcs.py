import numpy as np
import astropy.units as u
import astropy.constants as c

@u.quantity_input
def kick_kep_elements(delta_v:u.m/u.s, theta_i, phi_i, curr_epoch:u.year, P:u.year, tperi_i:u.year, a_i:u.m, e_i, inc_i:u.deg, omega_bar_i:u.deg, anode_i:u.deg, M:u.kg, m:u.kg):
	'''
	calculates the resulting orbital elements of a secondary given a velocity impulse is applied. 
	Adapted from section 2 and appendix A of Jackson et al. 2014
	
	Args:
        delta_v (float): magnitude of velocity kick
        theta_i, phi_i: theta and phi angles of velocity kick relative to the original orbit location IN RADIANS
        curr_epoch (float):  the time at which the explosion occurs, should just be the first in the list of epochs to evaluate
        P (np.array): orbital period (u.time)  
        tperi_i (float): initial epoch of periastron (u.time)
        a_i (float): semi-major axis of the initial orbit (u.distance)
        e_i (float): eccentricity of the initial orbit
        inc_i (float): inclination of the initial orbit  (u.angle)
        omega_bar_i (float): initial longitude of periastron (u.angle)
        anode_i (float): initial longitude of the ascending node (u.angle)
        M : Mass of the primary
        m : mass of the secondary

    Returns:
        tperi_prime (float): resulting epoch of periastron (u.time)
        a_prime (float): semi-major axis of the resulting orbit
        e_prime (float): eccentricity of the resulting orbit
        inc_prime (float): inclination of the resulting orbit  (u.angle)
        omega_prime (float): resulting longitude of periastron (u.angle)
        anode_prime (float): PA of the resulting ascending node (u.angle)
	'''
	
	#convert necessary values:
	f_i = tperi_to_f(tperi_i, curr_epoch, P, e_i)
	
	omega_i = omega_bar_i - anode_i #source: https://en.wikipedia.org/wiki/Longitude_of_the_periapsis
	
	theta1 = theta_i #redundancy in case I need to change the defs. of theta and phi
	phi1 = phi_i
	
	
	omega_i = omega_i.to_value(u.rad)
	inc_i = inc_i.to_value(u.rad)
	anode_i = anode_i.to_value(u.rad)
	
	#paper-defined values:
	
	#theta1 and phi1 trig
	Ct1 = np.cos(theta1)
	St1 = np.sin(theta1)
	Cp1 = np.cos(phi1)
	Sp1 = np.sin(phi1)
	
	#alpha and beta plus trig
	alpha = omega_i + f_i
	beta = phi1 - anode_i
	Cbeta = np.cos(beta)
	Sbeta = np.sin(beta)
	Calpha = np.cos(alpha)
	Salpha = np.sin(alpha)
	
	#input parameter trig
	Ci = np.cos(inc_i)
	Si = np.sin(inc_i)
	Comega = np.cos(omega_i)
	Somega = np.sin(omega_i)
	Cf = np.cos(f_i)
	Sf = np.sin(f_i)
	Sanode = np.sin(anode_i)
	Canode = np.cos(anode_i)
	
	#calculate theta and phi
	theta = np.arccos(Ct1*Ci - St1*Si*Sbeta) #Eqs A1
	phi = (St1*(Sbeta*Ci*Comega - Cbeta*Somega) + Ct1*Si*Comega) / (St1*(Sbeta*Ci*Somega + Cbeta*Comega) + Ct1*Si*Somega)
	
	#print('THETA: '+str(theta))
	
	#theta and phi trig
	St = np.sin(theta)
	Sp = np.sin(phi)
	Ct = np.cos(theta)
	Cp = np.cos(phi)
	
	Spf = np.sin(phi - f_i)
	
	#Other useful values
	vk = np.sqrt(c.G *(M+m)/a_i) #this line is why a needs to be given as unit of linear distance
	Q = vk/((1-e_i**2)**.5)
	
	#paper-defined vectors: 
	R_vec = (a_i*(1-e_i**2)/(1+e_i*Cf))*np.asarray([Cf,Sf,0]) #distance vector
	vel_i = Q*np.asarray([-Sf, (e_i+Cf), 0]) #initial velocity vector
	delta_v_vec = delta_v*np.asarray([St*Cp,St*Sp,Ct]) #impulse as a vector
	
	vel_prime = vel_i + delta_v_vec #velocities are additive
	
	h_i_vec = np.cross(R_vec, vel_i) #used two eqs here to test math was working out right, it is
	h_i = a_i*vk*(1-e_i**2)**.5
	
	h_prime_vec = np.cross(R_vec, vel_prime)
	h_prime = np.linalg.norm(h_prime_vec)
	
	#eq 4:
	frac = 1 - (delta_v/vk)**2 - (2/(np.sqrt(1-e_i**2)))*(delta_v/vk)*St*(Spf + e_i*Sp) #a/a'
	
	a_prime = a_i / frac
	
	#Using Kepler 3 to recalculate the period:
	P_prime = P*frac**(-3/2)
	
	#eq 6:
	e_prime = np.sqrt(1 - ((1 - e_i**2)*((h_prime/h_i)**2)*(a_i/a_prime)))
	
	#eq A2:
	factor = (np.sqrt(1 - e_i**2)/(1+e_i*Cf))*(delta_v/vk) #defined by me b/c it's used several times
	inc_prime = np.arccos((Ci + factor*St1*(Calpha*Sbeta - Salpha*Cbeta*Ci))*(h_prime**2/h_i**2)**-.5)
	
	#eq A3:
	top = Sanode*Si + factor* (Ct1*(Sanode*Calpha + Canode*Salpha*Ci) - St1*Cp1*Si*Salpha)
	bottom = Canode*Si + factor* (Ct1*(Sanode*Calpha - Canode*Salpha*Ci) - St1*Cp1*Si*Salpha)
	
	anode_prime = np.arctan(top/bottom)
	
	#eq 13:
	inner = (1/e_prime)*(((h_prime**2/h_i**2)*(1+e_i*Cf))-1)
	
	if inner < -1: #a rounding error occurs at multiples of pi when delta_v = 0, it's a rare case
		inner = -1 * inner.unit
	if inner > 1:
		inner = 1 * inner.unit
	
	f_prime = np.arccos(inner)
	
	#eq A4
	omega_prime = np.arcsin(Salpha*Si/np.sin(inc_prime)) - f_prime
	
	#convert back to tperi:
	tperi_prime = f_to_tperi(f_prime, e_prime, P_prime, curr_epoch)
	
	#tperi_prime = f_prime
	
	return P_prime, tperi_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime
	
	
def tperi_to_f(tperi, t, per, e):
	ma = ((2*np.pi/per) * (t-tperi)).value #mean anomaly, from Murray and Dermott
	
	f = ma + (2*e-.25*e**3)*np.sin(ma) + (5/4)*(e**2)*np.sin(2*ma) + (13/12)*(e**3)*np.sin(3*ma)#from https://en.wikipedia.org/wiki/True_anomaly#From_the_mean_anomaly
	
	return f

def f_to_tperi(f, e, P, t):
	
	f=f.value
	
	first = -np.sqrt(1-e**2)*np.sin(f)
	second = -e - np.cos(f)
	
	#print()
	
	ma = np.arctan2(first, second).value + np.pi - e*(-first/(1+e*np.cos(f))) #https://en.wikipedia.org/wiki/Mean_anomaly
	
	tperi = t - (ma*P/(2*np.pi)) #inverse of Murray and Dermott
	
	
	return tperi
	