# import sys
# sys.path.append("/Users/mohan/pycchem")

import pycchem.misc_utils as pm
import numpy as np
import matplotlib.pyplot as plt

nh3_url = "https://janaf.nist.gov/tables/H-083.html"

nh3_results = pm.nist_collector(nh3_url)

nist_temp = nh3_results[:, 0] # Kelvin
nist_entropy = nh3_results[:, 1] # kJ/mol*K
nist_enthalpy = nh3_results[:, 2] # kJ/mol

zero_k_nh3_entropy = nh3_results[0][1] # S at 0 K
zero_k_nh3_enthalpy = nh3_results[0][2] # H - H0 at 0 K

temp_step = 0.05

interp_temp = np.arange(200, 4000 + temp_step, temp_step) # interpolation temperature range in Kelvin

entropy_nh3 = pm.cs_interp(nh3_results[1:-1, 0], nh3_results[1:-1, 1], interp_temp) # [1:-1, 0] Ignore first and last element of first column (0, 298.15 K)

enthalpy_nh3 = pm.cs_interp(nh3_results[1:-1, 0], nh3_results[1:-1, 2], interp_temp)

fig, ax1 = plt.subplots(figsize=(8, 8))
ax2 = ax1.twinx() # initiate twin axes 

plt.setp(ax1.get_yticklabels(), fontweight='bold') # set y1 label weight to bold 

ax1.scatter(nist_temp, nist_entropy, color = 'royalblue') # plot entropy vs. temperature 

ax1.plot(interp_temp, entropy_nh3, color = 'k') # plot interpolated entropy
ax1.set_ylabel(r"Entropy $ \mathbf{ \left[ \frac{kJ}{mol \cdot K} \right] } $", fontsize = 14, fontweight = "bold", color = 'royalblue') # set left axis y-label
ax1.tick_params(axis="y", labelcolor='royalblue', color='royalblue', labelsize=12) # set label color and size

plt.setp(ax2.get_yticklabels(), fontweight='bold') # set y2 label weight to bold 

ax2.scatter(nist_temp, nist_enthalpy, color = 'red') # plot enthalpy vs. temperature 

ax2.plot(interp_temp, enthalpy_nh3, color = 'k') # plot interpolated enthalpy
ax2.set_ylabel(r"Enthalpy $ \mathbf{ \left[ \frac{kJ}{mol} \right] } $", fontsize = 14, fontweight = "bold", color = 'red') # set left axis y-label
ax2.tick_params(axis="y", labelcolor='red', color='red', labelsize=12) # set label color and size
ax2.spines['right'].set_color('red')  # Set y-axis (right) spine color
ax2.spines['left'].set_color('royalblue')  # Set y-axis (left) spine color; ax2 is instantiated after ax1 so we have to set the color here

plt.setp(ax1.get_xticklabels(), fontweight='bold')

ax1.set_xlabel("Temperature [K]", fontsize = 14, fontweight = "bold") # set x-axis aesthetics 
ax1.tick_params(axis="x", labelsize=12)

plt.xticks(fontsize=12, weight = 'bold')  # Adjust the fontsize as needed (e.g., 12, 14)
plt.show();