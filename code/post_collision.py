# import numpy as np
# import astropy.units as u
# import astropy.constants as c
# import matplotlib.pyplot as plt
# from astropy.visualization import quantity_support
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.collections import PolyCollection
# from kepler3_catherine import *
# 
# from kepler3_addon import *
# 
# @u.quantity_input
# 
# P     = 12.0 * u.year # years
# eperi = 2000 * u.year
# M1 = 1.0*u.Msun
# M2 = 1.0*u.Mjup 
# 
# 
# 
# e     = 0.5
# i     = 0 * u.deg
# anode = 0 * u.deg
# w     = 0 * u.deg
# 
# empty_data_directory()
# 
# start = 2000*u.year
# steps_per_year = 4/u.year
# 
# pre_yrs = 5*u.year
# post_yrs = 30*u.year
# 
# total_steps = ((pre_yrs+post_yrs)*steps_per_year).astype(int)
# 
# pre_steps = (pre_yrs*steps_per_year).astype(int)
# 
# post_extra_lines = total_steps - pre_steps
# 
# # make complete orbit for progenitor
# t_orbit_all = np.linspace(start,start+P,((start+P)*steps_per_year).astype(int))
# _, _, Xa, Ya, Za, _, _, _ = kep3d(t_orbit_all,P,eperi,a,e,i,w,anode)
# 
# poststeps = ((start+P)*steps_per_year).astype(int)
# save_kep3d_short(Xa, Ya, Za, filename = 'ogorbit.csv', overwrite=True, post_steps = poststeps)
# 
# # orbit up to point of explosion
# t_orbit = np.linspace(start,start+pre_yrs,pre_steps)
# X, Y, Xs, Ys, Zs, Xv, Yv, Zv = kep3d(t_orbit,P,eperi,a,e,i,w,anode)
# 
# save_kep3d_short(Xs, Ys, Zs, filename = 'progenitor.csv', overwrite=True, post_steps = post_extra_lines)
# 
# # make a ball of particles
# (r,theta,phi) = fibonacci_sphere(1, 1500, 0.2)
# (xf,yf,zf) = sph2cart(r.flatten(),theta.flatten(),phi.flatten())
# 
# 
# # texplosion
# t_explode = start+pre_yrs
# print(f'time of explosion is at {t_explode:9.3f}')
# 
# # calculate mean anomaly for t_explode
# 
# M_explode = u.rad * 2*np.pi*(t_explode-eperi)/P
# 
# # calculate the true anomaly for this mean anomaly
# 
# (E_explode, f_explode) = kepler_solve(e,M_explode)
# 
# v_expand = 3 * u.km/u.s
# 
# (P_prime, f_prime, a_prime, e_prime, inc_prime, omega_prime, anode_prime) = \
#     kick_kep_elements_2d(v_expand, theta.flatten(), phi.flatten(), \
#     P, f_explode, a, e, w, anode, 1*u.Msun, 1*u.Mjup)
# 
# 
# 
# # now we want to get the periastron times for all the exploded particle orbits...
# 
# # calculate eccentric anomaly from true anomaly
# E_prime = np.arctan2(np.sqrt(1 - e_prime**2) * np.sin(f_prime), e_prime + np.cos(f_prime))
# 
# # calculate the mean anomaly fom the eccentric anomaly
# M_prime = E_prime - e_prime*np.sin(E_prime)*u.rad
# 
# T_peri_prime = t_explode - P_prime*M_prime/(2*np.pi*u.rad)
# 
# # time post explosion
# t_after = post_yrs
# t_post = t_explode + t_after
# post_steps = (post_yrs * steps_per_year).astype(int)
# t_post = np.linspace(t_explode, t_post, post_steps)
# 
# # calculate positions of all exploded particles
# 
# # Xcloud = np.zeros_like(a_prime)
# # Ycloud = np.zeros_like(a_prime)
# # Zcloud = np.zeros_like(a_prime)
# 
# extra_lines = total_steps - post_steps
# 
# for u in np.arange(P_prime.size): # may God forgive me for using this for loop
#     _, _, Xcloud, Ycloud, Zcloud, _, _, _ = \
#         kep3d(t_post,P_prime[u],T_peri_prime[u],a_prime[u],e_prime[u],inc_prime[u],omega_prime[u],anode_prime[u])
#     
#     filename = 'cloud'+str(u)+'.csv'
#     save_kep3d_short(Xcloud, Ycloud, Zcloud, filename = filename, overwrite=True, empty_steps = extra_lines)
# 
# #animate_3d()
# # fig7 = plt.figure(figsize=(8,8))
# # ax5 = fig7.add_subplot(111, projection='3d')
# # 
# # 
# # 
# # 
# # ax5.plot(Xa, Ya, Za, marker=None, color='red', alpha=0.3)
# # 
# # ax5.plot(Xs, Ys, Zs, marker=None, color='red')
# # 
# # _, _, X2, Y2, Z2, _, _, _ = kep3d(t_post,P,eperi,a,e,i,w,anode)
# # ax5.scatter(X2, Y2, Z2, marker=None, color='red')
# # 
# # 
# # # 
# # ax5.scatter(xf, yf, zf, marker=None, color='green')
# # 
# # ax5.scatter(Xcloud, Ycloud, Zcloud, marker=None, color='blue')
# # 
# # ax5.scatter(0.0, 0.0, 0.0, marker='o', color='blue')
# # ax5.set_xlabel('X')
# # ax5.set_ylabel('Y')
# # ax5.set_zlabel('Z')
# # ax5.set_title(f'Cloud {t_after:5.3f} with expansion velocity {v_expand:4.1f}')
# # 
# # h = 10
# # ax5.set_xlim(-h,h)
# # ax5.set_ylim(-h,h)
# # ax5.set_zlim(h,-h)
# # 
# # plt.show()
