from pycchem import *
import numpy as np

infile_path = "/Users/mohan/Desktop/Research/paolucci/data/freq/pd2/square_planar/1_nh3_1_h2o_2_c2h4_c2h4180apart.log"

temp = np.linspace(100, 300, 100)

s_trans, s_vib, s_rot = calc_entropy(infile = infile_path, temperature = temp, use_low_freq = True)