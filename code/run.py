from simulation_funcs import *
from animate_kep3d import *
from lightcurves import *
from spheres import *
from fileio import *
from huge_sphere import *
from data_fitting import *


# make_huge_sphere(num_parts = 1e6)
## MAKE INITIAL RUN
# write_param_file(n_shells = 100, num_per_shell = 1000, v_i = 5.5*u.km/u.s, 
# 		kick_vel = 1.75 * u.km/u.s, delta_t= 5*u.year, n_steps = 500)
# 		
# steprate = simulate_spherical_docinput(overwrite = True)

##FILTER OUT IRRELEVANT STARS
# save_lightcurve_data_midsaves(cut_run = True)
# steprate = get_steprate()
# plot_from_saved(steprate=steprate, legend=False)
# 
# # animate_los_subset()

##SIMULATE EXTENDED DATA
# steprate = simulate_postprune(save_dir = './data_extended/', overwrite = True)

##CREATE MODEL LIGHT CURVES
# save_lightcurve_data_midsaves(directory = './data_extended/', foi_file = 'foi_extended.csv', lc_file = 'lc_data_extended')
# steprate = get_steprate()
# plot_from_saved(readfile = 'lc_data_extended.npz', steprate=steprate, legend=False)
#
# animate_los_subset(directory = './extra_data/', foi_list = 'foi_extra.csv')

##FIT REAL DATA:
trim_models(readfile = 'lc_data_extended.npz', writefile = 'trimmed_models.npz', buffer = 100)


### TESTING ZONE:
# time, flux, flux_err = read_real_data()
# timestep = get_steprate()
# 
# binned_flux, time_bins = bin_real_data (time, flux, timestep = timestep)

# plot_lc_basic(binned_flux)
# plot_lc_basic(flux)

fit_data(modelfile = 'trimmed_models.npz', datafile = './asassn_data/21qj_fromgithub.csv', filter_name = 'g', paramfile = 'params.csv')
