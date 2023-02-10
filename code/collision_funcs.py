import numpy as np
import astropy.units as u
import astropy.constants as c

@u.quantity_input
def kick_kep_elements(delta_v:u.m/u.s, theta_i:u.rad, phi_i:u.rad, P:u.year, f_i:u.rad, a_i:u.m, e_i, inc_i:u.deg, omega_bar_i:u.deg, anode_i:u.deg, M:u.kg, m:u.kg):
	'''
	calculates the resulting orbital elements of a secondary given a velocity impulse is applied. 
	Adapted from section 2 and appendix A of Jackson et al. 2014
	
	Args:
        delta_v (float): magnitude of velocity kick
        theta_i, phi_i: theta and phi angles of velocity kick relative to the original orbit location IN RADIANS
        P (np.array): orbital period (u.time)  
        f_i (float): initial true anomaly (u.angle)
        a_i (float): semi-major axis of the initial orbit (u.distance)
        e_i (float): eccentricity of the initial orbit
        inc_i (float): inclination of the initial orbit  (u.angle)
        omega_bar_i (float): initial longitude of periastron (u.angle)
        anode_i (float): initial longitude of the ascending node (u.angle)
        M : Mass of the primary
        m : mass of the secondary

    Returns:
        f_prime (float): resulting epoch of periastron (u.time)
        a_prime (float): semi-major axis of the resulting orbit
        e_prime (float): eccentricity of the resulting orbit
        inc_prime (float): inclination of the resulting orbit  (u.angle)
        omega_prime (float): resulting longitude of periastron (u.angle)
        anode_prime (float): PA of the resulting ascending node (u.angle)
	'''
	
	#convert necessary values:
	omega_i = omega_bar_i #- anode_i #source: https://en.wikipedia.org/wiki/Longitude_of_the_periapsis
	
	theta1 = theta_i #redundancy in case I need to change the defs. of theta and phi
	phi1 = phi_i
	
	#print(omega_i)
	
# 	omega_i = omega_i.to_value(u.rad)
# 	inc_i = inc_i.to_value(u.rad)
# 	anode_i = anode_i.to_value(u.rad)
	
	#paper-defined values:
	
	#theta1 and phi1 trig
	Ct1 = np.cos(theta1)
	St1 = np.sin(theta1)
	Cp1 = np.cos(phi1)
	Sp1 = np.sin(phi1)
	
	#alpha and beta plus trig
	#print(omega_i)
	#print(f_i)
	
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
# 	theta = theta_i
# 	phi = phi_i
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
	vk = np.sqrt(c.G *(M+m)/a_i) #this line is why a needs to be given as unit of linear distance
	Q = vk/((1-e_i**2)**.5)
	
	#paper-defined vectors: 
	R_vec = (a_i*(1-e_i**2)/(1+e_i*Cf))*np.asarray([Cf,Sf,0]) #distance vector
	vel_i = Q*np.asarray([-Sf, (e_i+Cf), 0]) #initial velocity vector
	delta_v_vec = delta_v*np.asarray([St*Cp,St*Sp,Ct]) #impulse as a vector
	
	vel_prime = vel_i + delta_v_vec #velocities are additive
	
#	h_i_vec = np.cross(R_vec, vel_i) #used two eqs here to test math was working out right, it is
	h_i = a_i*vk*(1-e_i**2)**.5
	
	h_prime_vec = np.cross(R_vec, vel_prime)
	h_prime = np.linalg.norm(h_prime_vec)
	
	factor = ((1-e_i**2)**.5/(1+e_i*Cf))*(delta_v/vk)
	
	h_frac_sq = 1 + 2 * factor * St*Spf + factor**2 * (Ct**2 + (St**2)*(Spf**2))#h'^2 / h^2, eq 8 
	
	#eq 4: calculates a
	frac = 1 - (delta_v/vk)**2 - (2/(np.sqrt(1-e_i**2)))*(delta_v/vk)*St*(Spf + e_i*Sp) #a/a'
	#print(P*frac)
	
	a_prime = a_i / frac
	
	#Using Kepler 3 to recalculate the period:
	P_prime = P*frac**(-3/2)
	
	#eq 6: calculates e
	e_prime = np.sqrt(1 - ((1 - e_i**2)*h_frac_sq*(a_i/a_prime)))
	
	#eq A2: calculates inc 
	inc_prime = np.arccos((Ci + factor*St1*(Calpha*Sbeta - Salpha*Cbeta*Ci))*(h_frac_sq)**-.5)
	
# 	#eq 10: calculates inc
# 	inc_prime = np.arccos((1 + factor*St*Spf) * h_frac_sq**-.5)
	
	#if theta
	
	#eq A3: calculates anode
	top = Sanode*Si + factor* (Ct1*(Sanode*Calpha + Canode*Salpha*Ci) - St1*Cp1*Si*Salpha)
	bottom = Canode*Si + factor* (Ct1*(Canode*Calpha - Sanode*Salpha*Ci) - St1*Cp1*Si*Salpha)
	
# 	print(top)
# 	print(bottom)
	
	anode_prime = np.arctan2(top,bottom)
	
# 	print(anode_prime)
	
	#eq 11: calcuates anode
# 	if theta.value <= np.pi/2:
# 		anode_prime = f_i
# 	else:
# 		anode_prime = f_i + np.pi * u.rad

	#anode_prime = f_i + np.pi * u.rad
	
	#eq 14: calculates f
	cosf = (1/e_prime)*((h_frac_sq*(1+e_i*Cf))-1)
	sinf = (1/e_prime)*(h_frac_sq**.5)*(e_i*Sf + (1-e_i**2)**.5 * (delta_v/vk)*St*Cpf)
	f_prime = np.arctan2(sinf, cosf)
	
	if f_prime < 0.0:
		f_prime = f_prime + 2*np.pi

	#eq A4: calculates omega
	if inc_prime == 0:
		print('Warning: inc_prime == 0 creates a divide by zero error when calculating omega, setting inc_prime = 1e-5')
		inc_prime = 1e-5 * inc_prime.unit
	
	top = Salpha*Si/np.sin(inc_prime)
	bottom = (1/np.cos(anode_prime)) * (Canode*Calpha - Sanode*Salpha*Ci + np.sin(anode_prime)*Salpha*(Si*np.cos(inc_prime)/np.sin(inc_prime)))
	omega_prime = np.arctan2(top,bottom) - f_prime
	
	#eq 12: calculates omega:
# 	if theta.value <= np.pi/2:
# 		omega_prime = -f_prime
# 	else:
# 		omega_prime = np.pi*u.rad - f_prime
		
	omega_bar_prime =  omega_prime #+anode_prime
	
# 	anode_prime -= (np.pi)*u.rad
			
	
	return P_prime, f_prime, a_prime, e_prime, inc_prime, omega_bar_prime, anode_prime
	
	
	