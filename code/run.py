from simulation_funcs import *
from animate_kep3d import *
from lightcurves import *
from spheres import *
from fileio import *
from huge_sphere import *


# make_huge_sphere(num_parts = 1e6)
# write_param_file(n_shells = 100, num_per_shell = 1000, v_i = 5.5*u.km/u.s, 
# 		kick_vel = 1.75 * u.km/u.s, delta_t= 5*u.year, n_steps = 500)
# 		
# steprate = simulate_spherical_docinput(overwrite = True)
# 
# save_lightcurve_data_midsaves(cut_run = True)
# steprate = get_steprate()
# plot_from_saved(steprate=steprate, legend=False)
# 
# # animate_los_subset()
# 
# steprate = simulate_postprune(save_dir = './data_extended/', overwrite = True)

# save_lightcurve_data_midsaves(directory = './data_extended/', foi_file = 'foi_extended.csv', lc_file = 'lc_data_extended')
steprate = get_steprate()
plot_from_saved(readfile = 'lc_data_extended.npz', steprate=steprate, legend=False)

shell_weights = np.linspace(0,1,100)

plot_from_saved(readfile = 'lc_data_extended.npz', shell_weights = shell_weights, plt_subs = True, steprate=steprate, legend=False)

# animate_los_subset(directory = './extra_data/', foi_list = 'foi_extra.csv')


# make_huge_sphere(num_parts = 1e6)
# 
# for num in [100, 1000]:
# 	write_param_file(n_shells = 5, num_per_shell = num, v_i = 5.5*u.km/u.s, 
# 		kick_vel = 1.75 * u.km/u.s, delta_t= 5*u.year, n_steps = 500)
# # 
# # 	simulate_spherical_docinput(save_dir = './data_'+str(int(num))+'/', overwrite = True)
# # 	save_lightcurve_data(directory = './data_'+str(int(num))+'/', foi_file = 'foi_'+str(int(num))+'.csv', cut_run = True)
# # 	
# # 	simulate_postprune(read_dir = './data_'+str(int(num))+'/', save_dir = './data_extended_'+str(int(num))+'/', foi_file = 'foi_'+str(int(num))+'.csv', overwrite = True)
# 	save_lightcurve_data(directory = './data_extended_'+str(int(num))+'/', foi_file = 'foi_'+str(int(num))+'.csv', lc_file = 'lc_'+str(int(num)))
# 	plot_from_saved(readfile = 'lc_'+str(int(num))+'.npz', save_plt = True, savefile = './lightcurve_'+str(int(num))+'.pdf')
# # 	staticdistribution(from_list = False, directory = './data_extended_'+str(int(num))+'/' , save = True, savefile = 'allpoints_'+str(int(num))+'.png')

