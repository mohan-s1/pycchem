import sys
sys.path.append("/Users/mohan/pycchem")

import pycchem.gaussian.entropy as pge
import pycchem.gaussian.gauss_utils as pgg
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


script_dir = Path(__file__).parent # find parent folder of current file, `gaussian_example.py`

infile_path = script_dir / "1_nh3_3_h2o.log"

temp = np.linspace(200, 500, 100)

s_trans, s_vib, s_rot = pge.calc_entropy(infile = infile_path, temperature = temp, use_low_freq = True) # all entropies in terms of J/mol*K

print(f"The energy of the optimized structure is {pgg.last_gauss_energy(infile_path)} Hartree")

plt.plot(temp, s_trans, label = 'Translation')
plt.plot(temp, s_vib, label = 'Vibration')
plt.plot(temp, s_rot, label = 'Rotation')
plt.xlabel("Temperature [K]", fontsize = 14, fontweight = "bold")
plt.ylabel(r"Entropy $ \mathbf{ \left[ \frac{J}{mol \cdot K} \right] } $", fontsize = 14, fontweight = "bold")
plt.xticks(fontsize=12, weight = 'bold')  # Adjust the fontsize as needed (e.g., 12, 14)
plt.yticks(fontsize=12, weight = 'bold')  # Adjust the fontsize as needed (e.g., 12, 14)
plt.legend(frameon = False)
plt.show();