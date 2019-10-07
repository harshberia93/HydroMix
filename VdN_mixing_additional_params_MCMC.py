#################################################################################################################################
'''
Created on: 30th September 2019 (Harsh Beria) 
Last updated on 
What it does?

Computes the proportion of groundwater that is made of snow vs rain using a MCMC sampler
Introduced additional model parameters (isotopic lapse rate)

SourceFiles used:
SourceFiles/Rain_SP_GW_VdN.xlsx
SourceFiles/Hypsometric_curve_data.xlsx	
	
OutputFiles processed:
OutputFiles/VdN/Lapse_rate/H2_results_MCMC_1000.csv
OutputFiles/VdN/Lapse_rate/O18_results_MCMC_1000.csv
	
Figures made: NONE
'''
#################################################################################################################################
# %% Imports

from __future__ import division
import pandas as pd
import numpy as np
import re

#################################################################################################################################
# %% Custom functions

def Random_walk(initialParam, paramLimit, step):
	'''
	Generates the new parameter set

	initialParam => Initial parameter values (LIST)
	paramLimit => Lower and upper limits of the list of parameters ([[lambdaLowVal, lambdaHighVal], [StdLowVal, StdHighVal])
	step => Defines the maximum extent of jump from the current parameter state (in %)

	Returns the updated parameter state (LIST)
	'''
	
	stepChange = [(i[1]-i[0])*0.005*step for i in paramLimit] # List of step changes (equal jump on both sides 0.01/2 = 0.005)
	
	# Lower and upper limits for all parameters
	lowerParLimit = [max(initialParam[index]-stepChange[index], paramLimit[index][0]) for index in range(len(initialParam))]
	upperParLimit = [min(initialParam[index]+stepChange[index], paramLimit[index][1]) for index in range(len(initialParam))]
	
	# Generating new parameter set by sampling uniformly in the given range
	newParam = [np.random.uniform(lowerParLimit[index], upperParLimit[index]) for index in range(len(initialParam))]
	
	return (newParam) # Returns the list containing new parameters

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


def HydroMix(source1, elev1, source2, elev2, mixture, stdDev, param_init, param_limit, nb_iter, HYPS_DIC, RandomWalkStep=5):
	'''
	Contains the core of HydroMix
	
	source1 => Concentration data of source 1
	elev1 => Elevation at which source 1 data is collected
	source2 => Concentration data of source 2
	elev2 => Elevation at which source 2 data is collected
	mixture => Concentration data of the mixture made up of a linear combination of sources 1 and 2
	stdDev => Standard deviation of the mixture (specified apriori and not calibrated)
	param_init => Initial value of the list of parameters ([lambdaValue, stdValue])
	param_limit => Lower and upper limits of the list of parameters ([[lambdaLowVal, lambdaHighVal]])
	nb_iter => Number of MCMC runs
	RandomWalkStep => In percentage to define the maximum extent of jump from the current parameter state
	HYPS_DIC => Key is elevation and value is the percent of catchment at that elevation

	Returns a tuple containing [[LOG LIKELIHOOD VALUES], [PARAMETER VALUES], [RESIDUAL VALUES]]
	
	'''

	logLikelihoodLis, paramLis, residual_list = [], [param_init], []
	std_param = stdDev
	n_accepts = 0
	for i in range(nb_iter):
		
		# Compute new parameter set
		updatedParam = Random_walk(paramLis[-1], param_limit, RandomWalkStep)
		lambda_param, lapse_param,  = updatedParam[0], updatedParam[1]
	
		# Incorporating lapse rate effect in the isotopic ratio of source1 (snow) and source2(rain)
		temp_source1, temp_source2 = [], []
		for index in range(len(source1)):
			temp_iso, temp_elev = source1[index], elev1[index]
			temp_source1.append(catchment_avg_isotope(temp_iso, temp_elev, lapse_param, HYPS_DIC))
		for index in range(len(source2)):
			temp_iso, temp_elev = source2[index], elev2[index]
			temp_source2.append(catchment_avg_isotope(temp_iso, temp_elev, lapse_param, HYPS_DIC))
	
		# Log likelihood computation
		temp_loglikelihood, temp_residual = [], []
		for temp_mix in mixture:	
			temp_value = 0.
			
			for temp_s1 in temp_source1:
				for temp_s2 in temp_source2:
					temp_estimated_mix = lambda_param*1.*temp_s1 + (1-lambda_param)*temp_s2*1.
					
					# Log likelihood computation
					temp_value -= (0.5* np.log(2 * np.pi * std_param * std_param))
					temp_value -= (0.5*(temp_estimated_mix - temp_mix)*(temp_estimated_mix - temp_mix)) / (std_param*std_param)
					temp_residual.append(temp_estimated_mix-temp_mix)
			
			temp_loglikelihood.append(temp_value)
		residual_list.append(np.average(temp_residual))
		LLValue = np.sum(temp_loglikelihood)
	
		# Hastings test
		if (i == 0): # First iteration (accept it)
			logLikelihoodLis.append(LLValue)
			n_accepts += 1
		else:
			alpha = np.exp(LLValue - logLikelihoodLis[-1])
			if ( (alpha > 1) or (np.random.rand() > (1-alpha)) ): # Accept the new move
				paramLis.append(updatedParam)
				logLikelihoodLis.append(LLValue)
				n_accepts += 1
			else: # Keep the last move as it is
				paramLis.append(paramLis[-1])
				logLikelihoodLis.append(logLikelihoodLis[-1])
				
		# For displaying purposes
		if (i%100 == 0):
			print ("Iteration number:" + str(i+1) + ", Acceptance: " + str(n_accepts) + ", Acceptance %: " + str(n_accepts/(i+1)))
		
	return ((logLikelihoodLis, paramLis, residual_list))

# This is a csv writer
import csv
def csv_writer(data, path):
	with open(path, "wb") as csv_file:
		writer = csv.writer(csv_file, delimiter=',')
		for line in data:
			writer.writerow(line)

#################################################################################################################################
# %% Main variables

NUMBER_ITERATIONS = 10000

# Bayesian model parameters (in the given range (inclusive))
LAMBDA_RANGE = [0., 1.] # LAMBDA values imply the fraction of snow in groundwater
O18_LAPSE_SLOPE = [-0.0081, 0.0081] # Obtained from GNIP Swiss data (-0.0027)
H2_LAPSE_SLOPE = [-0.0582, 0.0582] # Obtained from GNIP Swiss data (-0.0194)

SWISS_LAPSE = {"H2": -0.0194, "O18": -0.0027} # Swiss lapse rate according to GNIP data
SWISS_LAPSE_lowBound = {"H2": -0.029, "O18": -0.0039} # Swiss lapse rate according to GNIP data
SWISS_LAPSE_highBound = {"H2": -0.0076, "O18": -0.0013} # Swiss lapse rate according to GNIP data

JUMP_PERCENTAGE = 5. # In percentage (JUMP_PERCENTAGE/2 in both directions)

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
# %% Running HydroMix for H2

rain, snow, gw = rain_df["H2 isotope"].values, snow_df["H2 isotope"].values, gw_df["H2 isotope"].values
rain_elev, snow_elev = rain_df["Elevation (m)"].values, snow_df["Elevation (m)"].values

H2_std = np.std(gw_df["H2 isotope"].values, ddof=1) # Standard deviation of H2 in groundwater

# List of initial parameter values 
initParam = [np.random.uniform(LAMBDA_RANGE[0], LAMBDA_RANGE[1]), np.random.uniform(H2_LAPSE_SLOPE[0], H2_LAPSE_SLOPE[1])]

# Lower and upper limits of the model parameters
paramLimit = [LAMBDA_RANGE, H2_LAPSE_SLOPE]

# Running the mixing model
LOGLIKELIHOOD_H2, PARAM_H2, RESIDUAL_H2 = HydroMix(snow, snow_elev, rain, rain_elev, gw, H2_std, 
												   initParam, paramLimit, NUMBER_ITERATIONS, HYPS_DIC, JUMP_PERCENTAGE)

snowRatioLis_H2 = [i[0] for i in PARAM_H2]
lapseLis_H2 = [i[1] for i in PARAM_H2]

# %% Writing output in csv file

final_lis = [["Snow ratio", "Log likelihood", "Error std", "H2 lapse rate", "Residual"]]
path = OUTPUTFILEPATH + "H2_results_MCMC_" + str(NUMBER_ITERATIONS) + ".csv"
for index in range(0, len(LOGLIKELIHOOD_H2)):
	final_lis.append([ round(snowRatioLis_H2[index], 4), round(LOGLIKELIHOOD_H2[index], 4), round(H2_std, 4), 
				   round(lapseLis_H2[index], 4), round(RESIDUAL_H2[index], 4) ])
csv_writer(final_lis, path)
print (path)

# %% Running HydroMix for O18

rain, snow, gw = rain_df["O18 isotope"].values, snow_df["O18 isotope"].values, gw_df["O18 isotope"].values
rain_elev, snow_elev = rain_df["Elevation (m)"].values, snow_df["Elevation (m)"].values

O18_std = np.std(gw_df["O18 isotope"].values, ddof=1) # Standard deviation of O18 in groundwater

# List of initial parameter values 
initParam = [np.random.uniform(LAMBDA_RANGE[0], LAMBDA_RANGE[1]), np.random.uniform(O18_LAPSE_SLOPE[0], O18_LAPSE_SLOPE[1])]

# Lower and upper limits of the model parameters
paramLimit = [LAMBDA_RANGE, O18_LAPSE_SLOPE]

# Running the mixing model
LOGLIKELIHOOD_O18, PARAM_O18, RESIDUAL_O18 = HydroMix(snow, snow_elev, rain, rain_elev, gw, O18_std, 
												   initParam, paramLimit, NUMBER_ITERATIONS, HYPS_DIC, JUMP_PERCENTAGE)

snowRatioLis_O18 = [i[0] for i in PARAM_O18]
lapseLis_O18 = [i[1] for i in PARAM_O18]

# %% Writing output in csv file

final_lis = [["Snow ratio", "Log likelihood", "Error std", "O18 lapse rate", "Residual"]]
path = OUTPUTFILEPATH + "O18_results_MCMC_" + str(NUMBER_ITERATIONS) + ".csv"
for index in range(0, len(LOGLIKELIHOOD_O18)):
	final_lis.append([ round(snowRatioLis_O18[index], 4), round(LOGLIKELIHOOD_O18[index], 4), round(O18_std, 4),
				   round(lapseLis_O18[index], 4), round(RESIDUAL_O18[index], 4) ])
csv_writer(final_lis, path)
print (path)

##################################################################################################################################
