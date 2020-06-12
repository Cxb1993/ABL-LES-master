"""
shape of array
"""

# import needed libraries
import numpy as np
import os

# import constants
from constants import cnst

# input folder path
input_folder = os.path.join(cnst.root_home,'ABL-LES-master','get_info','input')

print("\n  Getting information for files in input folder.")

# for each file in the input folder
for filename in os.listdir(input_folder):
    
    # print filename
    print(f"\n For {filename},")

    # load the matrix
    loaded_mat = np.loadtxt(os.path.join(input_folder,filename))
    
    # get info on shape and unique values
    shape = np.shape(loaded_mat)
    unique = np.unique(loaded_mat)
    
    # print info to screen
    print(f"\n    The shape of {filename} is {shape}",)
    print(f"    The values in {filename} are {unique}",)
    

# end
print('\n  All files have been accounted for. Exiting....\n')
    



#### for when I need to input constants by going up a folder:
#os.chdir(os.path.join("D:",os.sep,"surface-heterogeneity-analysis"))
#from constants import cnst
#os.chdir(os.path.join("create_surfaces"))