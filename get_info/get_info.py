"""
shape of array
"""

import numpy as np
import os

print("\n  Getting information for files in input folder.")

root = os.path.join(os.sep,'glade','u','home','jfogart','get-info-on-arrays')

for filename in os.listdir(os.path.join(root,'input')):
    
    # print filename
    print(f"\n For {filename},")

    # load the matrix
    loaded_mat = np.loadtxt(os.path.join(root,'input',filename))
    
    # get info
    shape = np.shape(loaded_mat)
    unique = np.unique(loaded_mat)
    
    # print info
    print(f"\n    The shape of {filename} is {shape}",)
    print(f"    The values in {filename} are {unique}",)
    
    

# print done
print('\n  All files have been accounted for. Exiting....\n')
    
