import numpy as np
import astropy.units as u
import astropy.constants as c


@u.quantity_input
def kick_kep_elements(delta_v:u.m/u.s, theta_i:u.rad, phi_i:u.rad, P:u.year, f_i:u.rad, a_i:u.m, e_i, inc_i:u.deg, omega_i:u.deg, anode_i:u.deg, M:u.kg, m:u.kg):
	'''
	calculates the resulting orbital elements of a secondary given a velocity impulse is applied. 
	Adapted from seCtion 2 and appendix A of Jackson et al. 2014
	
	Args:
        delta_v (float): magnitude of velocity kick
        theta_i, phi_i: theta and phi angles of velocity kick relative to the original orbit location IN RADIANS
        P (np.array): orbital period (u.time)  
        f_i (float): initial true anomaly (u.angle)
        a_i (float): semi-major axis of the initial orbit (u.diStance)
        e_i (float): eccentricity of the initial orbit
        inc_i (float): inclination of the initial orbit  (u.angle)
        omega_bar_i (float): initial longitude of periaStron (u.angle)
        anode_i (float): initial longitude of the ascending node (u.angle)
        M : Mass of the primary
        m : mass of the secondary

    Returns:
        f_prime (float): resulting epoch of periaStron (u.time)
        a_prime (float): semi-major axis of the resulting orbit
        e_prime (float): eccentricity of the resulting orbit
        inc_prime (float): inclination of the resulting orbit  (u.angle)
        omega_prime (float): resulting longitude of periaStron (u.angle)
        anode_prime (float): PA of the resulting ascending node (u.angle)
	'''
	
	#convert necessary values:
	theta1 = theta_i 
	phi1 = phi_i
	
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
	theta = np.arccos(Ct1*Ci - St1*Si*Sbeta)#Eqs A1
	phi = np.arctan2((St1*(Sbeta*Ci*Comega - Cbeta*Somega) + Ct1*Si*Comega), (St1*(Sbeta*Ci*Somega + Cbeta*Comega) + Ct1*Si*Somega))
	
	#theta and phi trig
	St = np.sin(theta)
	Sp = np.sin(phi)
	Ct = np.cos(theta)
	Cp = np.cos(phi)
	
	Spf = np.sin(phi - f_i)
	Cpf = np.cos(phi - f_i)
	
	#Other useful values
	dvvk = delta_v / np.sqrt((c.G*(M+m))/a_i)
	hprimh2 = 1 + 2*np.sqrt(1.0-e_i*e_i)*dvvk*St*Spf/(1.+e_i*Cf) + (1-e_i*e_i)*dvvk*dvvk*(Ct*Ct + St*St*Spf*Spf)/((1.+ e_i*Cf)**2)
	
	faCtor = ((1-e_i**2)**.5/(1+e_i*Cf))*dvvk
	
	#eq 4: calculates a
	frac = 1 - dvvk**2 - (2/(np.sqrt(1-e_i**2)))*dvvk*St*(Spf + e_i*Sp) #a/a'
	a_prime = a_i / frac
	
	#Using Kepler 3 to recalculate the period:
	P_prime = ((P*P)*(a_prime*a_prime*a_prime)/(a_i*a_i*a_i))**(1./2)
	
	#eq 6: calculates e
	e_prime = np.sqrt(1 - ((1 - e_i**2)*hprimh2*(a_i/a_prime)))
	
	#eq A2: calculates inc 
	inc_prime = np.arccos((Ci + faCtor*St1*(Calpha*Sbeta - Salpha*Cbeta*Ci))*(hprimh2)**-.5)
	
	#eq A3: calculates anode
	top = Sanode*Si + faCtor* (Ct1*(Sanode*Calpha + Canode*Salpha*Ci) - St1*Cp1*Si*Salpha)
	bottom = Canode*Si + faCtor* (Ct1*(Canode*Calpha - Sanode*Salpha*Ci) - St1*Cp1*Si*Salpha)
	
	anode_prime = np.arctan2(top,bottom)
	
	#eq 14: calculates f
	coSf = (1/e_prime)*((hprimh2*(1+e_i*Cf))-1)
	sinf = (1/e_prime)*(hprimh2**.5)*(e_i*Sf + (1-e_i**2)**.5 * (dvvk)*St*Cpf)
	f_prime = np.arctan2(sinf, coSf)
	
	if f_prime < 0.0:
		f_prime = f_prime + 2*np.pi*u.rad

	#eq A4: calculates omega
	if inc_prime == 0:
		print('Warning: inc_prime == 0 creates a divide by zero error when calculating omega, setting inc_prime = 1e-5')
		inc_prime = 1e-5 * inc_prime.unit
	
	top = Salpha*Si/np.sin(inc_prime)
	bottom = (1/np.cos(anode_prime)) * (Canode*Calpha - Sanode*Salpha*Ci + np.sin(anode_prime)*Salpha*(Si*np.cos(inc_prime)/np.sin(inc_prime)))
	omega_prime = np.arctan2(top,bottom) - f_prime
			
	
	return P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime
	
@u.quantity_input
def kick_kep_elements_2d(delta_v:u.m/u.s, theta_i:u.rad, phi_i:u.rad, P_i:u.year, f_i:u.rad, a_i:u.m, e_i, omega_i:u.deg, anode_i:u.deg, M:u.kg, m:u.kg):
	'''
	calculates the resulting orbital elements of a secondary given a velocity impulse is applied. 
	Adapted from seCtion 2 and appendix A of Jackson et al. 2014
	
	Args:
        delta_v (float): magnitude of velocity kick
        theta_i, phi_i: theta and phi angles of velocity kick relative to the original orbit location IN RADIANS
        P (np.array): orbital period (u.time)  
        f_i (float): initial true anomaly (u.angle)
        a_i (float): semi-major axis of the initial orbit (u.diStance)
        e_i (float): eccentricity of the initial orbit
        inc_i (float): inclination of the initial orbit  (u.angle)
        omega_i (float): initial longitude of periaStron (u.angle)
        anode_i (float): initial longitude of the ascending node (u.angle)
        M : Mass of the central star
        m : mass of the planet pre-explosion
        pm: pass of the particle

    Returns:
        f_prime (float): resulting epoch of periaStron (u.time)
        a_prime (float): semi-major axis of the resulting orbit
        e_prime (float): eccentricity of the resulting orbit
        inc_prime (float): inclination of the resulting orbit  (u.angle)
        omega_prime (float): resulting longitude of periaStron (u.angle)
        anode_prime (float): PA of the resulting ascending node (u.angle)
	'''
	
	#trig values
	Spf = np.sin(phi_i - f_i, dtype = np.float128)
	Cpf = np.cos(phi_i - f_i, dtype = np.float128)
	Sf = np.sin(f_i, dtype = np.float128)
	Cf = np.cos(f_i, dtype = np.float128)
	
	St = np.sin(theta_i, dtype = np.float128)
	Ct = np.cos(theta_i, dtype = np.float128)
	
	Sp = np.sin(phi_i, dtype = np.float128)
	Cp = np.cos(phi_i, dtype = np.float128)
	
	#dimensionless initial velocity
	vorb=np.sqrt(c.G * (M+m)/ a_i, dtype = np.float128)
	
	dvvk=(delta_v/vorb)
	#print(dvvk.value.type)
	
	#calculate output a
	a_prime = a_i/(1. - dvvk*(dvvk + (2./np.sqrt(1.-e_i*e_i))*St*(Spf + e_i*Sp)))
	#print(a_prime)
	
	#useful factor:
	hprimh2 = np.float128(1 + 2*np.sqrt(1.0-e_i*e_i)*dvvk*St*Spf/(1.+e_i*Cf) + (1-e_i*e_i)*dvvk*dvvk*(Ct*Ct + St*St*Spf*Spf)/((1.+ e_i*Cf)**2))
	
	#calculate output eccentricity
	e_prime = np.sqrt(1.-(1-e_i*e_i)*hprimh2*a_i/a_prime)
	
	#calculate output inclination
	inside = (1 + np.sqrt(1-e_i*e_i)*dvvk*St*Spf/(1+e_i*Cf)) / np.sqrt(hprimh2)
	
	# if inside > 1:
# 		inside -= 2
	inc_prime = np.arccos(inside)
	
	#print((1 + np.sqrt(1-e_i*e_i)*dvvk*St*Spf/(1+e_i*Cf))/np.sqrt(hprimh2))
	
	#calculate output true anomaly (f):
	sinf_prime = (1./e_prime)*np.sqrt(hprimh2)*(e_i*Sf + np.sqrt(1-e_i*e_i)*St*Cpf*dvvk)
	
	cosf_prime = (1./e_prime)*(hprimh2*(1+e_i*Cf)-1)
	
	f_prime = np.arctan2(sinf_prime, cosf_prime)
	
	#AdjuSt for correCt range:
	if f_prime < 0:
		f_prime += 2*np.pi*u.rad
	
	#calculate output longitude of ascending node (anode):
	if Ct >= 0:
		anode_prime = f_i
	else:
		anode_prime = f_i+np.pi*u.rad
		
	#calculate output arg. of periaStron (omega)
	if Ct >= 0:
		omega_prime = np.pi*u.rad - f_prime
		if omega_prime < 0:
			omega_prime+=np.pi*2*u.rad
	else:
		omega_prime = 2*np.pi*u.rad - f_prime
		
	#calculate output period:
	P_prime = ((P_i*P_i)*(a_prime*a_prime*a_prime)/(a_i*a_i*a_i))**(1./2)
	
	return P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime
