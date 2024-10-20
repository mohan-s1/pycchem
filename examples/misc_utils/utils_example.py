import sys
sys.path.append("/Users/mohan/pycchem")

import pycchem.misc_utils as pm
import matplotlib.pyplot as plt

nh3_url = "https://janaf.nist.gov/tables/H-083.html"

nh3_results = pm.nist_collector(nh3_url)

nist_temp = nh3_results[:, 0]
nist_entropy = nh3_results[:, 1]
nist_enthalpy = nh3_results[:, 2]

fig, ax1 = plt.subplots(figsize=(8, 8))
ax2 = ax1.twinx()

ax1.scatter(nist_temp, nist_entropy, color = 'royalblue') # plot entropy vs. temperature 

# ax1.set_yticks(fontsize=12, weight = 'bold')  # Adjust the fontsize as needed (e.g., 12, 14)
ax1.set_ylabel(r"Entropy $ \mathbf{ \left[ \frac{kJ}{mol \cdot K} \right] } $", fontsize = 14, fontweight = "bold", color = 'royalblue')
ax1.tick_params(axis="y", color = 'royalblue')

ax2.scatter(nist_temp, nist_enthalpy, color = 'red') # plot entropy vs. temperature 
ax2.set_ylabel(r"Enthalpy $ \mathbf{ \left[ \frac{kJ}{mol} \right] } $", fontsize = 14, fontweight = "bold", color = 'red')

ax1.set_xlabel("Temperature [K]", fontsize = 14, fontweight = "bold")
plt.xticks(fontsize=12, weight = 'bold')  # Adjust the fontsize as needed (e.g., 12, 14)
plt.show();