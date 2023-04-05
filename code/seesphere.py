import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from spheres import *
		

n_layers = 2
n_particles = 100
max_vel = 5 *u.km/u.s

# vel, theta, phi = random_sphere(n_layers, n_particles, max_vel)
vel, theta, phi = fibonacci_sphere(n_layers, n_particles, max_vel)

print(theta)

draw_spherical(vel,theta,phi)



