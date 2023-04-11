from simulation_funcs import *
from animate_kep3d import *
from lightcurves import *
from spheres import *
from fileio import *
from huge_sphere import *

from tqdm import tqdm

# write_param_file(n_shells = 5, num_per_shell = 50, v_i = 5.5*u.km/u.s, 
# 		kick_vel = 1.75 * u.km/u.s, delta_t= 5*u.year, n_steps = 100)
		
# simulate_spherical_docinput(overwrite = True)
# 
# save_lightcurve_data(cut_run = True)
# animate_los_subset()
#
make_huge_sphere(num_parts = 1e4)
simulate_postprune(overwrite = True)
# 
# save_lightcurve_data()
# animate_los_subset()
# 
# plot_from_saved()


# make_huge_sphere(num_parts = 1e8)