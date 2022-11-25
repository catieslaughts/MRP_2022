import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from kepler3 import *
from kepler3_addon import *
		

n_layers = 3
n_particles = 1000
max_vel = 5 *u.km/u.s

vel, theta, phi = random_sphere(n_layers, n_particles, max_vel)
vel, theta, phi = fibonacci_sphere(n_layers, n_particles, max_vel)

draw_spherical(vel,theta,phi)



