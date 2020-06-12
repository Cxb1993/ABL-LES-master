"""
LES Averaging - A simple script to obtain output from LES and save these
as .txt files (arrays and vectors) to be post-processed how you wish
"""

# Import needed libraries
import numpy as np
import os, sys

# import constants
from constants import cnst

########## Preliminaries - fill in from PARAM.F90 ##########

# IMPORTANT! Be sure to change these parameters from param.f90
project = cnst.project
save_text_path = cnst.save_text_path
scratch_path = cnst.scratch_path
N_x = cnst.N_x
N_z = cnst.N_z 
z_i = cnst.z_i
dt_dim = cnst.dt_dim
total_time = cnst.total_time
p_count = cnst.p_count
u_g = cnst.u_g
theta_s = cnst.theta_s
t_slices = cnst.t_slices
nt = cnst.nt

# start analysis
while True:
    print(f'\n  Conducting analysis for {project}')
    print(f'  Loading from {scratch_path}')
    print(f'  Saving information to {save_text_path}')
    cont = input("\n  With this information, continue? y/n: ")
    while cont.lower() not in ("y","n"):
        cont = input("  Please choose, y or n: ")
    if cont == "y":
        break
    if cont == "n":
        sys.exit()


print(f"\n  For {project}, nt = {nt} with {t_slices} total time slices")

# First and last timesteps for the averaging and animation
i0 = int(input("    First time slice you want to see: ")) # first timestep slice you want to see
ie = int(input("    Last time slice you want to see:  ")) # last timestep slice you want to see

# the final string that this codes goes to get data
# make sure this is right:
data_loc_string = os.path.join(scratch_path,'output')

############### Retriving data ##################

# velocity components, stress
var_list = ['u', 'v', 'w', 'tke', 'txz', 'tyz']
# potential temp, heat flux, subgrid heat flux
scalar_list = ['theta', 'wt', 'sgs_t3']

tot_list = var_list + scalar_list

for var in tot_list:

    # load the file, check if this variable exists
    filename = 'aver_' + var + '.out'
    filepath = os.path.join(data_loc_string,filename)
    print(f"\n  Now retrieving data for {var} from {i0} to {ie} time slices.")
    print(f"    Loading from {filepath}")
    loaded_mat = np.loadtxt(filepath)

    # Check for NaN values
    if np.isnan(loaded_mat).any() == True:
        print("\n  WARNING!")
        slice_row = np.argwhere(np.isnan(loaded_mat))[0,0] # first row with NaN
        print(f"    There are NaN values in the array at row {slice_row}. Slicing these off...")
        print("    Stopping Analysis, Breaking.")
        exit()

    # Reshape to a 3D numpy array, each slice is in time
    mat_3d = np.reshape(loaded_mat[:,1:], (-1, N_z, N_x))

    # Print some statistics about the data
    # print(f"\n  General Statistics for {var}:")
    # print(f"    The shape of the new data array is {np.shape(mat_3d)}")
    # print(f"    There are {np.shape(mat_3d)[0]} total timeslices")
    # print(f"    This simulation ran for {total_time_dim/3600.0:.2f} hours")
    # print(f"    The lifetime of one large eddy is {large_eddy_lifetime/3600.0:.2f} hours")
    # print(f"    ~{total_time_dim/large_eddy_lifetime/20.0:.2f} large eddies have been simulated")

    # Average the results to 2D
    if ie >= int(np.shape(mat_3d)[0]):
        print("\n  Warning: Last time slice exceeds final time slice in LES run")
        print("    Exiting...")
        exit()

    # Save the time-averaged 2d view of the data from i0 to ie
    time_avg = np.mean(mat_3d[i0:ie], axis=0)
    np.savetxt(os.path.join(save_text_path,'txts',f'{project}_{var}_nd_yt_{i0}_{ie}_tsteps.txt'),
               time_avg,fmt='%.10f',delimiter=' ')
    print(f"\n  2D average saved for {var}")

    # Save the time and x averaged vector
    time_x_avg = np.mean(time_avg,axis=1)
    np.savetxt(os.path.join(save_text_path,'txts',f'{project}_{var}_nd_yt_{i0}_{ie}_tsteps.txt'),
               time_x_avg,fmt='%.10f',delimiter=' ')
    print(f"  1D average saved for {var}")

print(f"\n  LES averaging complete, exiting...\n")
