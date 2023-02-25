import numpy as np
import astropy.units as u
import astropy.constants as c
from os import path
import os

from kepler3 import *

def kick_kep_elements(delta_v:u.km/u.s, theta_i:u.rad, phi_i:u.rad, 
    P_i:u.year, f_i:u.rad, a_i:u.au, e_i, omega_i:u.deg, anode_i:u.deg,
     M:u.kg, m:u.kg):
    '''
    calculates the resulting orbital elements of a secondary given a velocity impulse is applied. 
    Adapted from section 2 and appendix A of Jackson et al. 2014
    
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

    Returns:
    	P_prime (float): resulting Period (u.time)
        f_prime (float): resulting epoch of periastron (u.time)
        a_prime (float): semi-major axis of the resulting orbit
        e_prime (float): eccentricity of the resulting orbit
        inc_prime (float): inclination of the resulting orbit  (u.angle)
        omega_prime (float): resulting longitude of periaStron (u.angle)
        anode_prime (float): PA of the resulting ascending node (u.angle)
    '''
    
    #np.seterr(all='raise')
    
    #trig values
    Spf = np.sin(phi_i - f_i)
    Cpf = np.cos(phi_i - f_i)
    Sf = np.sin(f_i)
    Cf = np.cos(f_i)
    
    St = np.sin(theta_i)
    Ct = np.cos(theta_i)
    
    Sp = np.sin(phi_i)
    Cp = np.cos(phi_i)
    
    #dimensionless initial velocity (line after equation 3)
    vorb=np.sqrt(c.G * (M+m)/ a_i).to(u.km/u.s)
    
    dvvk=(delta_v/vorb)
    
    #calculate output a
    a_prime = a_i/(1. - dvvk*(dvvk + (2./np.sqrt(1.-e_i*e_i))*St*(Spf + e_i*Sp)))
    
    #a_prime = 1 *u.au
    
    #useful factor: equation 8
    hprimh2 = np.float128(1 + 2*np.sqrt(1.0-e_i*e_i)*dvvk*St*Spf/(1.+e_i*Cf) + \
        (1-e_i*e_i)*dvvk*dvvk*(Ct*Ct + St*St*Spf*Spf)/((1.+ e_i*Cf)**2))
    
    #calculate output eccentricity: equation 6
    e_prime = np.sqrt(1.-(1-e_i*e_i)*hprimh2*a_i/a_prime)
    # for idx, e in enumerate(e_prime):
#     	if e > 1:
#     		#print('FUCK')
#     		e_prime[idx] = np.nan
    
    #calculate output inclination: equation 10
    inside = (1 + np.sqrt(1-e_i*e_i)*dvvk*St*Spf/(1+e_i*Cf)) / np.sqrt(hprimh2)
    
    inc_prime = np.arccos(inside)
    
#     try:
#     	inc_prime = np.arccos(inside)
#     except:
#     	inc_prime = np.nan
#     	print('Invalid value encountered in arccosine')
    
    #calculate output true anomaly (f): equation 14
    sinf_prime = (1./e_prime)*np.sqrt(hprimh2)*(e_i*Sf + np.sqrt(1-e_i*e_i)*St*Cpf*dvvk)
    
    cosf_prime = (1./e_prime)*(hprimh2*(1+e_i*Cf)-1)
    
    f_prime = np.arctan2(sinf_prime, cosf_prime) # is this true???? yes, see below

    f_prime = np.mod(f_prime, 2*np.pi*u.rad)
    
    #calculate output longitude of ascending node (anode): equation 11
    tltpi2 = (Ct >= 0) # mask for theta < pi/2

    # equation 11
    anode_prime = f_i*np.ones_like(Ct)+np.pi*u.rad
    anode_prime[tltpi2] = f_i

    # equation 12
    omega_prime = np.pi*u.rad -f_prime # f_prime is an output array and has full size
    omega_prime[tltpi2] = -f_prime[tltpi2]


    omega_prime = np.mod(omega_prime, 2*np.pi*u.rad)

    P_prime = P_i*np.power(a_prime/a_i,3./2.)

    return P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime


##New Kep3d:
# def kep3d_trueanom(epoch:u.year, P:u.year, true_anom, a, e, inc:u.deg, omega:u.deg, anode:u.deg):
#     """
#     Calculate the position and velocity of an orbiting body
# 
#     Given the Kepler elements for the secondary about the primary
#     and in the coordinate frame where the primary is at the origin
#     
#     The same as the kep3d function in kepler3.py, but now takes in true anomaly instead of tperi
# 
#     Args:
#         epoch (np.array):  epochs to evaluate (u.time)
#         P (np.array): orbital period (u.time) 
#         true_anom (float): the true anomaly of the orbit
#         a (float): semi-major axis of the orbit
#         e (float): eccentricity of the orbit
#         inc (float): inclination of the orbit  (u.angle)
#         omega (float): longitude of periastron (u.angle)
#         anode (float): PA of the ascending node (u.angle)
# 
#     Returns:
#        X,Y, Xs,Ys,Zs, Xsv,Ysv,Zsv
# 
#     Output frame has X,Y in computer plotting coordinates
#     i.e. X is to the right, increasing (due West)
# 
#     Primary body is fixed at the origin.
# 
#     X,Y (float): 2D coordinates of in plane orbit with periapse
#                  towards the +ve X axis.
# 
#     Xs,Ys,Zs (float): The 3D coordinates of the secondary body
#         in the Position/velocity coords frame.
# 
#     Xsv, Ysv, Zsv (float): The 3D velocity of the secondary body
#         in the Position/velocity coords frame.
# 
#     Coordinate frames are shown below.
# 
#     The 3D axes are NOT the usual right-handed coordinate frame. The
#     Observer is located far away on the NEGATIVE Z axis. This is done
#     so that +ve Zsv gives positive velocities consistent with most
#     astronomers idea of redshift being positive velocity values.
# 
# 
#     Sky coords         Computer coords   Position/velocity coords
# 
#       North                   Y                +Y    +Z
#         ^                     ^                 ^   ^
#         |                     |                 |  /
#         |                     |                 | /
#         |                     |                 |/
#         +-------> West        +-------> X       +-------> +X
#                                                /
#                                               /
#                                              /
#                                            -Z
# 
#     +Y is North, +X is West and +Z is away from the Earth
#     so that velocity away from the Earth is positive
# 
#     NOTE: Right Ascension is INCREASING to the left, but the
#     (X,Y,Z) outputs have RA increasing to the right, as seen
#     in the Computer coords. This is done to make plotting easier
#     and to remind the user that for the sky coords they need
#     to plot (-X,Y) and then relabel the RA/X axis to reflect 
#     what is seen in the sky.
# 
#     Taken from ``Lecture Notes on Basic Celestial Mechanics''
#     by Sergei A. Klioner (2011) page 22
#     http://astro.geo.tu-dresden.de/~klioner/celmech.pdf
# 
#     Note that mean motion is defined on p.17 and that
#     (P^2/a^3) = 4pi^2/kappa^2 and that the document is missing
#     the ^2 on the pi.
# 
#     Written:
#         Matthew Kenworthy, 2017
# 
#     """
# 
#     # mean motion n
#     n = 2 * np.pi / P
# 
#     # calc eccentric anomaly E
#     v = true_anom
#     E = np.arctan2(np.sqrt(1 - e**2) * np.sin(v), e + np.cos(v))
#     
# 
#     # calculate position and velocity in the orbital plane
#     cE = np.cos(E)
#     sE = np.sin(E)
#     surde = np.sqrt(1 - (e*e))
# 
#     X  = a * (cE - e)
#     Y  = a * surde * sE
# 
#     Xv = -(a * n * sE) / (1 - e * cE)
#     Yv =  (a * n * surde * cE) / (1 - e * cE)
# 
#     # calculate Euler rotation matrix to get from orbital plane
#     # to the sky plane
# 
#     mat = euler(anode, omega, inc)
# 
#     # rotate the coordinates from the orbital plane
#     # to the sky projected coordinates
# 
#     # TODO we lose the dimensionality of X and Y and it
#     # needs to be put back artificially
#     # problem with building the np.array below and putting the Quantity through
#     # the np.dot() routine
#     
#     (Xe, Ye, Ze)    = np.dot(mat, np.array([X.value,Y.value,0],dtype=object))
#     #blog = np.array([X,Y,np.zeros(X.size)]) * X.unit
#     #(Xe, Ye, Ze)    = np.dot(mat, blog)
#     (Xev, Yev, Zev) = np.dot(mat, np.array([Xv.value,Yv.value,0],dtype=object))
#     
#     ##CATHERINE added dtype=object to these two lines to avoid the depreciation warning
# 
#     Xs = -Ye * X.unit
#     Ys =  Xe * X.unit
#     Zs =  Ze * X.unit
#     Xsv = -Yev * Xv.unit
#     Ysv =  Xev * Xv.unit
#     Zsv =  Zev * Xv.unit
# 
#     # To transform from the Kep3D code coordinate system to the
#     # celestial sphere, where:
#     # Y is North, X is West and Z is away from the Earth
#     # we have
#     #   X,Y,Z,Xdot,Ydot,Zdot = -Ys, Xs, Zs, -Yv, Xv, Zv
# 
#     return(X,Y,Xs,Ys,Zs,Xsv,Ysv,Zsv)

def tperi_to_Tanom(tperi, epoch, P, e, derror=1e-6):
	# at epoch = tperi, mean anomoly M is 0
    # 
    # Y = time since epoch periastron
    Y = epoch - tperi

    # mean anomaly M varies smoothly from 0 to 2pi every orbital period
    # convert Y to angle in radians between 0 and 2PI
    Mt = Y / P
    M  = (2 * np.pi * (Mt - np.floor(Mt)))*u.radian
    
    (E,v) = kepler_solve(e, M, derror) #returned in rad units
    
    #print(E)
    
    return E, v 
    
def Tanom_to_tperi(f, e, P):
	#calc E
	top = np.tan(f/2) * np.sqrt(1-e)
	bottom = np.sqrt(1+e)
	E = 2*np.arctan2(top, bottom) #new-orbital.pdf pg 6
	
	#calc M
	M = E - e*np.sin(E)*u.rad #kepler's equation
	
	#calc tperi:
	tperi = -M*P/(2*np.pi*u.rad)
	#print(tperi)
	
	return tperi
