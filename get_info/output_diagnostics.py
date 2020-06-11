'''
a code to check for NaN values in a simulation
to determine whether or not to keep going
'''

# import libraries
import numpy as np
import os
import pandas as pd

# load project to run diagnostics
proj_path = os.path.join('maptest','maptest_beaufo_2000_aug31_192cube')
scratch = os.path.join(os.sep,'glade','scratch','jfogart')
total_path = os.path.join(scratch,proj_path)
print(f'\n  Checking on project: {total_path}')
nt = np.loadtxt(os.path.join(total_path,'total_time.dat'))
print(f'\n  For this project, nt = {nt}')

# relevant parameters for CFL
dt_dim = 0.03
Nx = 192
z_i = 1000.0
u_g = 8.0
p_count = 1000
total_time = nt*dt_dim

Ny = Nx; Nz = Nx; Lz = z_i
dx = 2.0*np.pi/(Nx); dy = 2.0*np.pi/(Ny)
dz = (Lz/z_i)/(Nz)
asp_ratio = dx/dz
Lx = Lz*asp_ratio; Ly = Lz*asp_ratio
dx_dim = dx*Lx; dy_dim = dy*Ly; dz_dim = dz

#### FIRST: check for NaN values using u,v,w, and theta ####

# check certain variables
var_to_check = ['u', 'v', 'w', 'theta']

for v in var_to_check:

    # load matrix
    print(f'\n  Checking variable: {v}')
    path = os.path.join(total_path,'output','aver_'+v+'.out')
    loaded_mat = np.loadtxt(path)

    #check for nan values, and where they are
    nan_loc = np.argwhere(np.isnan(loaded_mat))

    # if no NaN values...
    if nan_loc.size == 0:
        print('\n    There are no NaN values.')
        print('    Enjoy your day! :)')

    # if NaN values...
    else:
        print('\n    NaN values detected!')
        print(f'    The first NaN value occurs at row {nan_loc[0][0]}')

    #### THEN get the shape of these output arrays ####
    #### in case save_texts.py has errors ####
    print(f"    The shape of all {v} is {np.shape(loaded_mat)}")
    print(f"    Nz is given by {np.shape(loaded_mat)[0]/(np.shape(loaded_mat)[1]-1):.2f}")


#### THEN  estimate max CFL ####

# get max velocity from arrays above
#all_velocity = np.dstack((u_all,v_all,w_all))
##max_vel = u_g*np.max(all_velocity)
#
## get min dx
#min_dx = min(dx_dim, dy_dim, dz_dim)
#
## calculate and print CFL
#CFL_est = (max_vel*dt_dim)/min_dx
#print(f'\n  An estimate for CFL is {CFL_est:.3f}\n')

#### give output for time slices ####

# get information
time_slices = nt/p_count
index = np.arange(time_slices+1)
tt_index_s = p_count*dt_dim*index
tt_index_h = tt_index_s/3600

# convert to pd dataframe
output_info_list = pd.DataFrame(
    {'output_number': index,
     'total_time_s': tt_index_s,
     'total_time_h': tt_index_h
    })

# print pandas datafram to csv
output_info_list.to_csv('timing_output_info.csv',index=False,float_format='%.2f')
print("\n  CSV printed for output timing information.")
print("  For more information, 'vi timing_output_info.csv'")


print('\n  Finished, exiting program\n')
