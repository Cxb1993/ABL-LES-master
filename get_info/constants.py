""""
A class of cosntants that will be used for LES output analysis
Please use this to change load path, project name, etc

Ideally this would read straight from param.f90, but for now,
manual input is required (nbd)

"""

import numpy as np
import os

### constants class ###

class cnst(object):
    
    # root paths (don't touch these)
    root_scratch = os.path.join(os.sep,'glade','scratch','jfogart')
    root_home = os.path.join(os.sep,'glade','u','home','jfogart')
    
    # project-specific path
    project = 'beaufo_2000_aug31'
    scratch_path = os.path.join(root_scratch,'maptest',project)
    save_path = os.path.join(root_home,'LES-results','maptest_beaufo_2000_aug31')
    
    # spatial mesh info
    N_z = 192
    N_x = N_z
    z_i = 1000.0


    # temporal mesh info
    dt_dim = 0.03
    
    # output mesh info
    c_count = 100
    p_count = 1000
    
    # other constants and ND parameters
    vonk = 0.4
    u_g = 8.0
    u_star = u_g/15.0
    rho = 1.0
    theta_s = 300

total_time = nt*dt_dim