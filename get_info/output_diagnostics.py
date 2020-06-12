'''
a code to check for NaN values in a simulation
to determine whether or not to keep going
'''

# import libraries
import numpy as np
import os
import pandas as pd
import sys

# import constants
from constants import cnst

# import constants
project = cnst.project
od_path = cnst.od_path
scratch_path = cnst.scratch_path
N_x = cnst.N_x
N_z = cnst.N_z 
z_i = cnst.z_i
dt_dim = cnst.dt_dim
total_time = cnst.total_time
p_count = cnst.p_count
u_g = cnst.u_g

# start analysis
while True:
    print(f'\n  Conducting analysis for {project}')
    print(f'  Loading from {scratch_path}')
    print(f'  Saving information to {od_path}')
    cont = input("\n  With this information, continue? y/n: ")
    while cont.lower() not in ("y","n"):
        cont = input("  Please choose, y or n: ")
    if cont == "y":
        break
    if cont == "n":
        sys.exit()


# load project to run diagnostics
nt = int(np.loadtxt(os.path.join(scratch_path,'total_time.dat')))
print(f'\n  For this project, nt = {nt}')

#### FIRST: check for NaN values using u,v,w, and theta ####

# check certain variables
var_to_check = ['u', 'v', 'w', 'theta']

for v in var_to_check:

    # load matrix
    print(f'\n  Checking variable: {v}')
    path = os.path.join(scratch_path,'output','aver_'+v+'.out')
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
