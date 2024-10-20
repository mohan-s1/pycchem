import sys
sys.path.append("/Users/mohan/pycchem")

import pycchem.gaussian as pg
import pycchem.vasp as pv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


script_dir = Path(__file__).parent # find parent folder of current file, `gaussian_example.py`

infile_path = script_dir / "OUTCAR"

temp = np.linspace(200, 500, 100)

frequencies = pv.vasp_low_freq(infile_path)

print(f"The energy of the optimized structure is {pv.last_vasp_energy(infile_path)} eV")

s_vib = pg.vib_entropy(frequencies_cm1 = frequencies, temperatures = temp) # function from `entropy.py` in gaussian directory 

plt.plot(temp, s_vib)
plt.xlabel("Temperature [K]", fontsize = 14, fontweight = "bold")
plt.ylabel(r"Entropy $ \mathbf{ \left[ \frac{J}{mol \cdot K} \right] } $", fontsize = 14, fontweight = "bold")
plt.xticks(fontsize=12, weight = 'bold')  # Adjust the fontsize as needed (e.g., 12, 14)
plt.yticks(fontsize=12, weight = 'bold')  # Adjust the fontsize as needed (e.g., 12, 14)
plt.show();