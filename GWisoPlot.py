#################################################################################################################################
'''
Created on: 3rd December 2018 (Harsh Beria) 
Last updated on 
What it does?

Used for not so important plots

SourceFiles used:
OutputFiles/GW_conceptual/	
	
OutputFiles processed: NONE
Figures made:
OutputFiles/GW_conceptual/Figures/GW_iso_0.3_rain_0.6_snow.jpeg
	
'''
#################################################################################################################################
# %% Imports

from __future__ import division
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

#################################################################################################################################
# %% Custom functions


#################################################################################################################################
# %% Main variables

RAIN_EFF, SNOW_EFF = 0.3, 0.6
PATH = "OutputFiles/GW_conceptual/"

#################################################################################################################################
# %% Plotting isotopic ratio in groundwater for the entire simulation period

filename = PATH + "RAIN_" + str(RAIN_EFF) + "_SNOW_" + str(SNOW_EFF) + ".csv"
df = pd.read_csv(filename)
GW_iso = df["Storage isotopic ratio"].values

fig, ax = plt.subplots(1)
ax.plot(GW_iso)

# To display yticks in multiple of 10 years
ax.set_xticks(np.arange(0, 365*100+1, 365*10))
ax.set_xticklabels([year for year in range(0, 101, 10)])

# Setting the labels
ax.set_ylabel("Groundwater " + r'$\delta^{2}$H' + u' \u2030 ', fontsize=14)
ax.set_xlabel("Year number", fontsize=14)

plt.tight_layout()
plt.grid(linestyle='dotted')
figOutput = PATH + "Figures/GW_iso_" + str(RAIN_EFF) + "_rain_" + str(SNOW_EFF) + "_snow.jpeg"
plt.savefig(figOutput, dpi=300)
plt.close()
print (figOutput)







#################################################################################################################################
