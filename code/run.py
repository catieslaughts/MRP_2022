from simulation_funcs import *
from animate_kep3d import *

from lightcurves import *

from spheres import *

steprate = simulate_spherical_velinput(n_shells = 5, num_per_shell = 300,
	v_i = 5.5*u.km/u.s, kick_vel = 1.75 * u.km/u.s, delta_t= 5*u.year, n_steps = 100, overwrite = True)
# # # animate_3d(save=False)
# # 
# # 
# # 
# steprate = 200./u.year
# # 
save_lightcurve_data(R = .3 * u.AU, num_shells = 5)
# # plot_from_saved(steprate = steprate, shell_weights = [1, 1, 1, 1, 1], save_plt = False)
# 
# #animate_3d(save=False)
# 
# #staticdistribution(save = True)
# 
animate_los_subset()

#fibonacci_sphere(1, 10000000, [1]*u.km/u.s)