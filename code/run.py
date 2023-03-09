from simulation_funcs import *
from animate_kep3d import *

simulate_spherical(n_shells = 1, e_i = .5, num_per_shell = 100, tperi_i = 2025*u.year, 
	delta_t= 20*u.year, n_steps = 100, show_pre = True, overwrite = True)
animate_3d(save=False)
