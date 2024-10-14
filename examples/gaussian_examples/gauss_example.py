import sys
sys.path.append("/Users/mohan/pycchem")

from pycchem.gaussian import *
import numpy as np
from pathlib import Path

script_dir = Path(__file__).parent # find parent folder of current file, `gaussian_example.py`

infile_path = script_dir / "1_nh3_3_h2o.log"

temp = np.linspace(100, 300, 100)

s_trans, s_vib, s_rot = calc_entropy(infile = infile_path, temperature = temp, use_low_freq = True)