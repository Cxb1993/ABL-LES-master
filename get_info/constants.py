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
    project = 'maptest_beaufo_2000_aug31'
    scratch_path = os.path.join(root_scratch,'maptest',project)
    save_text_path = os.path.join(root_home,'LES-results',project)
    od_path = os.path.join(root_home,'get_info','timing_output_info.csv')
    
    # spatial mesh info
    N_z = 192
    N_x = N_z
    N_y = N_x
    z_i = 1000.0
    Lz = z_i
    dx = 2.0*np.pi/(N_x)
    dy = 2.0*np.pi/(N_y)
    dz = (Lz/z_i)/(N_z)
    asp_ratio = dx/dz
    Lx = Lz*asp_ratio
    Ly = Lz*asp_ratio
    dx_dim = dx*Lx
    dy_dim = dy*Ly
    dz_dim = dz

    # temporal info
    dt_dim = 0.03
    nt = 1000000
    #nt = int(np.loadtxt(os.path.join(scratch_path,'total_time.dat')))
    total_time = nt*dt_dim
    t_slices = int(nt/p_count)
    
    # output mesh info
    c_count = 100
    p_count = 1000
    
    # other constants and ND parameters
    vonk = 0.4
    u_g = 8.0
    u_star = u_g/15.0
    rho = 1.0
    theta_s = 300

