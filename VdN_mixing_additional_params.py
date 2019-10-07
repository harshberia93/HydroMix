#################################################################################################################################
'''
Created on: 12th November 2018 (Harsh Beria) 
Last updated on 
What it does?

Computes the proportion of groundwater that is made of snow vs rain
Introduced additional model parameters (isotopic lapse rate)

SourceFiles used:
SourceFiles/Rain_SP_GW_VdN.xlsx
SourceFiles/Hypsometric_curve_data.xlsx	
	
OutputFiles processed:
OutputFiles/VdN/Lapse_rate/H2_results.csv
OutputFiles/VdN/Lapse_rate/O18_results.csv
	
Figures made:
OutputFiles/VdN/Lapse_rate/lapse_rate_posterior.jpeg
OutputFiles/VdN/Lapse_rate/posterior.jpeg

'''
#################################################################################################################################
# %% Imports

from __future__ import division
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt

#################################################################################################################################
# %% Custom functions

def catchment_avg_isotope(isotope_ratio, elev_sample_point, slope_lapse_rate, hypsometric_dic):
	'''
	Computes the catchment averaged isotope ratio from isotope_ratio at a point using the lapse rate and hypsometric curve data
	
	isotope_ratio => Ratio of isotope at a given point
	elev_sample_point => Elevation at which the sample was taken and isotope ratio was obtained 
	slope_lapse_rate => slope of the lapse rate line
	hypsometric_dic => Key is elevation and value is the percent of catchment at that elevation

	Returns the catchment averaged isotope ratio
	'''

	isotope_dic = [0., 0.] # First is for sum and second is the normalizing factor

	for elevation in hypsometric_dic:
		estimated_ratio = isotope_ratio + slope_lapse_rate*1.*(elevation - elev_sample_point)
		isotope_dic[0] += (estimated_ratio*hypsometric_dic[elevation]*1.)
		isotope_dic[1] += hypsometric_dic[elevation]

	return (isotope_dic[0]/isotope_dic[1])


def HydroMix(source1, elev1, source2, elev2, mixture, LAMBDA_prior, ErrorSTD_prior, NUMBER_ITERATIONS, lapse_rate_param, HYPS_DIC):
	'''
	Contains the core of HydroMix
	
	source1 => Concentration data of source 1
	elev1 => Elevation at which source 1 data is collected
	source2 => Concentration data of source 2
	elev2 => Elevation at which source 2 data is collected
	mixture => Concentration data of the mixture made up of a linear combination of sources 1 and 2
	LAMBDA_prior => prior distribution of the proportion of source1 in the mixture
	ErrorSTD_prior => prior distribution of the error standard deviation
	NUMBER_ITERATIONS => Number of times HydroMix has be run
	lapse_rate_param => List of potential isotopic lapse rates
	HYPS_DIC => Key is elevation and value is the percent of catchment at that elevation

	Returns a tuple containing log likelihood values, lambda values and error standard deviations for all the model runs
	
	'''

	LIKELIHOOD = []
	for i in range(0, NUMBER_ITERATIONS, 1):
		# For displaying purposes
		if (i%100 == 0):
			print ("Iteration number:" + str(i+1))
		temp_source1, temp_source2 = [], [] # Values for source 1 and source 2 isotopic ratio after transformation
		
		# Incorporating lapse rate effect in the isotopic ratio of source1 (snow) and source2(rain)
		for index in range(len(source1)):
			temp_iso, temp_elev = source1[index], elev1[index]
			temp_source1.append(catchment_avg_isotope(temp_iso, temp_elev, lapse_rate_param[i], HYPS_DIC))
		for index in range(len(source2)):
			temp_iso, temp_elev = source2[index], elev2[index]
			temp_source2.append(catchment_avg_isotope(temp_iso, temp_elev, lapse_rate_param[i], HYPS_DIC))
			
		# Computing likelihood
		temp_likelihood = []
		for temp_mix in mixture:	
			temp_value = 0.
			temp_residual = []
			for temp_s1 in temp_source1:
				for temp_s2 in temp_source2:
					temp_estimated_mix = LAMBDA_prior[i]*1.*temp_s1 + (1-LAMBDA_prior[i])*temp_s2*1.
					
					# Likelihood computation
					temp_value += (-1 * (temp_estimated_mix - temp_mix)**2) / (2 * (ErrorSTD_prior[i]**2))
					temp_value -= (0.5* np.log(2 * np.pi * ErrorSTD_prior[i] * ErrorSTD_prior[i]))
					temp_residual.append(temp_estimated_mix-temp_mix)
			
			temp_likelihood.append(temp_value)
		LIKELIHOOD.append(np.sum(temp_likelihood))
		
	# Sorting the parameter sets by best runs (highest likelihood values)
	zipped = sorted(zip(LIKELIHOOD, LAMBDA_prior, ErrorSTD_prior, lapse_rate_param), reverse = True)
	LIKELIHOOD, LAMBDA_prior, ErrorSTD_prior, lapse_rate_param = zip(*zipped)
	
	return ((LIKELIHOOD, LAMBDA_prior, ErrorSTD_prior, lapse_rate_param))

# This is a csv writer
import csv
def csv_writer(data, path):
	with open(path, "wb") as csv_file:
		writer = csv.writer(csv_file, delimiter=',')
		for line in data:
			writer.writerow(line)

#################################################################################################################################
# %% Main variables

NUMBER_ITERATIONS = 5000

# Bayesian model parameters (in the given range (inclusive))
LAMBDA_RANGE = [0., 1.] # LAMBDA values imply the fraction of snow in groundwater
O18_LAPSE_SLOPE = [-0.0081, 0.0081] # Obtained from GNIP Swiss data (-0.0027)
H2_LAPSE_SLOPE = [-0.0582, 0.0582] # Obtained from GNIP Swiss data (-0.0194)

SWISS_LAPSE = {"H2": -0.0194, "O18": -0.0027} # Swiss lapse rate according to GNIP data
SWISS_LAPSE_lowBound = {"H2": -0.029, "O18": -0.0039} # Swiss lapse rate according to GNIP data
SWISS_LAPSE_highBound = {"H2": -0.0076, "O18": -0.0013} # Swiss lapse rate according to GNIP data

# Number of best simulations using which lambda is computed
BEST_SIM_PER = 5. # In percentage

HYPS_DIC = {} # Key is elevation and value is the percent of catchment at that elevation

OUTPUTFILEPATH = "OutputFiles/VdN/Lapse_rate/"


# %% Setting up a random seed
np.random.seed(1) # Setting up a common seed number

#################################################################################################################################
# %% Reading the hypsometric curve data

filename = "SourceFiles/Hypsometric_curve_data.xlsx"
df = pd.read_excel(filename, sheetname='Sheet1')
for index, row in df.iterrows():
	temp_elev_lis = re.findall(r'\d+', row['Elevation band']) # Identifying the numbers in the string and putting them in a list
	temp_elev_lis = [float(k) for k in temp_elev_lis] # Converting characters to floats
	elevation_avg = sum(temp_elev_lis)*1./len(temp_elev_lis) # Average elevation in a given elevation band
	HYPS_DIC[elevation_avg] = float(row['Percentage of grids'])


# %% Reading all the isotopic data

filename = "SourceFiles/Rain_SP_GW_VdN.xlsx"
df = pd.read_excel(filename, sheet_name='Sheet1')

# Separating the isotopic ratios in rain, snow and groundwater
rain_df = df[(df["Type"] == "Rain")]
snow_df = df[(df["Type"] == "Snow")]
gw_df = df[(df["Type"] == "Groundwater")]

#################################################################################################################################
# %% Initializing the model parameters

LAMBDA_params = np.random.uniform(LAMBDA_RANGE[0], LAMBDA_RANGE[1], NUMBER_ITERATIONS)
O18_lapse_param_lis = np.random.uniform(O18_LAPSE_SLOPE[0], O18_LAPSE_SLOPE[1], NUMBER_ITERATIONS)
H2_lapse_param_lis = np.random.uniform(H2_LAPSE_SLOPE[0], H2_LAPSE_SLOPE[1], NUMBER_ITERATIONS)

# Assuming constant error variance computed from data
LIKELIHOOD_std_params_H2 = np.full(NUMBER_ITERATIONS, np.std(gw_df["H2 isotope"].values, ddof=1))
LIKELIHOOD_std_params_O18 = np.full(NUMBER_ITERATIONS, np.std(gw_df["O18 isotope"].values, ddof=1))

# %% Running HydroMix for H2

rain, snow, gw = rain_df["H2 isotope"].values, snow_df["H2 isotope"].values, gw_df["H2 isotope"].values
rain_elev, snow_elev = rain_df["Elevation (m)"].values, snow_df["Elevation (m)"].values
LIKELIHOOD_H2, LAMBDA_H2, LIKELIHOOD_std_params_H2, H2_lapse_param_lis = HydroMix(snow, snow_elev, rain, rain_elev, gw, 
	   LAMBDA_params, LIKELIHOOD_std_params_H2, NUMBER_ITERATIONS, H2_lapse_param_lis, HYPS_DIC)

# %% Writing output in csv file

final_lis = [["Snow ratio", "Log likelihood", "Error std", "H2 lapse rate"]]
path = OUTPUTFILEPATH + "H2_results.csv"
for index in range(0, len(LIKELIHOOD_H2)):
	final_lis.append([ round(LAMBDA_H2[index], 2), round(LIKELIHOOD_H2[index], 2), round(LIKELIHOOD_std_params_H2[index], 2), 
				   round(H2_lapse_param_lis[index], 2) ])
csv_writer(final_lis, path)
print (path)

# %% Running HydroMix for O18

rain, snow, gw = rain_df["O18 isotope"].values, snow_df["O18 isotope"].values, gw_df["O18 isotope"].values
LIKELIHOOD_O18, LAMBDA_O18, LIKELIHOOD_std_params_O18, O18_lapse_param_lis = HydroMix(snow, snow_elev, rain, rain_elev, gw, 
	 LAMBDA_params, LIKELIHOOD_std_params_O18, NUMBER_ITERATIONS, O18_lapse_param_lis, HYPS_DIC)

# %% Writing output in csv file

final_lis = [["Snow ratio", "Log likelihood", "Error std", "O18 lapse rate"]]
path = OUTPUTFILEPATH + "O18_results.csv"
for index in range(0, len(LIKELIHOOD_O18)):
	final_lis.append([ round(LAMBDA_O18[index], 2), round(LIKELIHOOD_O18[index], 2), round(LIKELIHOOD_std_params_O18[index], 2),
				   round(O18_lapse_param_lis[index], 2) ])
csv_writer(final_lis, path)
print (path)

# %% Histogram plot showing snow ratio in groundwater using H2 and O18

plt.figure()
plt.hist(LAMBDA_H2[0:int(0.01 * BEST_SIM_PER * NUMBER_ITERATIONS)], color='blue', alpha=0.4, label=r'$\delta^{2}$H' + u' \u2030 (VSMOW)', normed='True')
plt.hist(LAMBDA_O18[0:int(0.01 * BEST_SIM_PER * NUMBER_ITERATIONS)], color='red', alpha=0.4, label=r'$\delta^{18}$O' + u' \u2030 (VSMOW)', normed='True')
plt.xlim(0., 1.)
plt.grid(linestyle='dotted')
plt.xlabel("Fraction of snow in groundwater", fontsize=14)
plt.ylabel("Normalised frequency", fontsize=14)
plt.legend()
plt.tight_layout()
path = OUTPUTFILEPATH + "posterior.jpeg"
plt.savefig(path, dpi=300)
plt.close()
print (path)

# %% Histogram plot showing posterior distributions of isotopic lapse rates in 2H and 18O

f, axarr = plt.subplots(2, figsize=(10, 10))
axarr[0].hist(H2_lapse_param_lis[0:int(0.01 * BEST_SIM_PER * NUMBER_ITERATIONS)], color='blue', alpha=0.4, normed='True')
axarr[0].set_xlim(H2_LAPSE_SLOPE[0], H2_LAPSE_SLOPE[1])
axarr[0].set_ylabel("Normalised frequency", fontsize=14)
axarr[0].set_xlabel("Lapse rate in" + r'$\ ^{2}$H', fontsize=14)
axarr[0].grid(linestyle='dotted')
#axarr[0].axvline(x=SWISS_LAPSE["H2"], color='black', label="Swiss lapse rate in" + r'$\ ^{2}$H')
axarr[0].axvspan(SWISS_LAPSE_lowBound["H2"], SWISS_LAPSE_highBound["H2"], facecolor='green', alpha=0.3, label="Swiss lapse rate range in" + r'$\ ^{2}$H')
axarr[0].legend(fontsize=14)
axarr[0].tick_params(labelsize=14)

axarr[1].hist(O18_lapse_param_lis[0:int(0.01 * BEST_SIM_PER * NUMBER_ITERATIONS)], color='red', alpha=0.4, normed='True')
axarr[1].set_xlim(O18_LAPSE_SLOPE[0], O18_LAPSE_SLOPE[1])
axarr[1].set_ylabel("Normalised frequency", fontsize=14)
axarr[1].set_xlabel("Lapse rate in" + r'$\ ^{18}$O', fontsize=14)
axarr[1].grid(linestyle='dotted')
#axarr[1].axvline(x=SWISS_LAPSE["O18"], color='black', label="Swiss lapse rate in" + r'$\ ^{18}$O')
axarr[1].axvspan(SWISS_LAPSE_lowBound["O18"], SWISS_LAPSE_highBound["O18"], facecolor='green', alpha=0.3, label="Swiss lapse rate range in" + r'$\ ^{18}$O')
axarr[1].legend(fontsize=14)
axarr[1].tick_params(labelsize=14)

plt.tight_layout()


path = OUTPUTFILEPATH + "lapse_rate_posterior.jpeg"
plt.savefig(path, dpi=300)
plt.close()
print (path)

#################################################################################################################################
