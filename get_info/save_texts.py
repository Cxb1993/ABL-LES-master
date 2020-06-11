"""
LES Averaging - A simple script to obtain output from LES and save these
as .txt files (arrays and vectors) to be post-processed how you wish
"""

# Import needed libraries
import numpy as np
import os

########## Preliminaries - fill in from PARAM.F90 ##########

# project to analyze - corresponds to output folder in /glade/scratch/jfogart
#project = str(input("\n  Project name: "))

# IMPORTANT! Be sure to change these parameters from param.f90
N_z = 192
dt_dim = 0.03
N_x = N_z
z_i = 1000
vonk = 0.4
p_count = 1000

# project name and file path for the solutions
project = 'maptest_beaufo_2000_aug31_192cube'
filepath = os.path.join(os.sep,'glade','scratch','jfogart','maptest',project)

# nt is taken from total_time.dat
tt_file = open(os.path.join(filepath,'total_time.dat'), 'r')
nt = int(float(tt_file.readlines()[0]))
tt_file.close()


t_slices = int(nt/p_count)
print(f"\n  For {project}, nt = {nt} with {t_slices} total time slices")

# First and last timesteps for the averaging and animation
i0 = int(input("    First time slice you want to see: ")) # first timestep slice you want to see
ie = int(input("    Last time slice you want to see:  ")) # last timestep slice you want to see

# Non-dimensional parameters
u_g = 2.0
u_star = u_g
rho = 1.0
theta_s = 300

# Large eddies
total_time_dim = nt*dt_dim # in seconds
large_eddy_lifetime = z_i/(u_star) # rough estimate, in seconds

# also add a subfolder here for retirving data, if need be
data_loc_string = os.path.join(filepath,'output')
#data_loc_string = f"/glade/scratch/jfogart/striptest/{project}/output/"

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
    np.savetxt(f"LES-results/{project}/txts/{project}_{var}_nd_yt_{i0}_{ie}_tsteps.txt",
               time_avg,fmt='%.10f',delimiter=' ')
    print(f"\n  2D average saved for {var}")

    # Save the time and x averaged vector
    time_x_avg = np.mean(time_avg,axis=1)
    np.savetxt(f"LES-results/{project}/txts/{project}_{var}_nd_xyt_{i0}_{ie}_tsteps.txt",
               time_x_avg,fmt='%.10f',delimiter=' ')
    print(f"  1D average saved for {var}")

print(f"\n  LES averaging complete, exiting...\n")
