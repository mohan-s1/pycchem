import sys
sys.path.append("/Users/mohan/pycchem")

import pycchem.gaussian.entropy as pge
import pycchem.vasp.vasp_utils as pvv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


script_dir = Path(__file__).parent # find parent folder of current file, `gaussian_example.py`

infile_path = script_dir / "1_nh3_3_h2o.log"

temp = np.linspace(200, 500, 100)

frequencies = pvv.vasp_parse_low_frequencies("infile")

s_vib = pge.vibrational_entropy(frequencies_cm1 = frequencies, temperatures = temp)
