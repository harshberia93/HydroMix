#################################################################################################################################
'''
Created on: 5th October 2018 (Harsh Beria) 
Last updated on: 22th November 2018 (Put the mixing part in another script GW_mixing_HydroMix.py) 
Last updated on: 21st December 2018 (Plots figure for snowfall vs snowmelt isotopic ratios) 
What it does?

Creates synthetic time series of rain and snow isotopic ratio, & simulates isotopic ratio in groundwater
Makes boxplots of isotopic ratio of snowfall vs snowmelt

SourceFiles used: NONE
OutputFiles processed:
OutputFiles/GW_conceptual/	
	
Figures made: 
OutputFiles/GW_conceptual/SF_SM_boxplot.jpeg

'''
#################################################################################################################################
# %% Imports

from __future__ import division
import numpy as np
import random
import matplotlib
import matplotlib.pyplot as plt

#################################################################################################################################
# %% Custom functions

def airTemp_Gen(meanTemp=4., amplTemp=8., years=1., offset=-np.pi/2.):
	'''
	To generate daily time series of air temperature
	Assumes a sine wave and adds a normally distributed error with variance = 0.2 * amplitude of temperature
	
	meanTemp => Mean annual temperature
	amplTemp => Amplitude of the sine wave of the temperature time series
	years => Number of years for which the time series is to be obtained
	offset => To take into account that January temperature is lowest and July temperature is highest
	
	Returns a numpy array containing the time series of airTemp
	'''
	
	numb_of_days = int(365 * years)
	dayNumb = np.linspace(1, numb_of_days, numb_of_days) # Corresponding to the number of days in a year
	airTemp = amplTemp * np.sin(2 * np.pi * dayNumb / 365. + offset) + meanTemp + np.random.normal(loc=0, scale=0.2 * amplTemp, size=numb_of_days)
	return (airTemp)
	
def PoissonPrcp(NumbEvent=30, meanPrcp=1500., years=1.):
	'''
	To generate daily time series of precipitation assuming the time between precipitation events comes from a poisson distribution
	& precipitation amount comes from an exponential distribution
	
	NumbEvent => Number of precipitation events in a year
	meanPrcp => Average annual precipitation
	years => Number of years for which time series is to be obtained
	
	Returns a numpy array with daily precipitation values
	'''
	intermittent_times = np.random.poisson(lam = 365. / NumbEvent, size = int(round(years * NumbEvent)))
	precip_days = np.cumsum(intermittent_times)
	precip_days = precip_days[(precip_days > 1) & (precip_days < int(round(years * 365)))] # Removing dates outside the possible range
	precip_amt = np.random.exponential(scale = meanPrcp/NumbEvent, size = len(precip_days))
	
	prcp = np.zeros(int(years * 365))
	for index in range(len(precip_days)):
		prcp[precip_days[index] - 1] = precip_amt[index]
	return(prcp)

def prcp_Iso(meanIso=-80., amplIso=40., years=1., offset=-np.pi/2.):
	'''
	To generate time series of precipitation isotopic ratio
	Assumes a sine wave and adds a normally distributed error with variance = 0.2 * amplitude of isotopic ratio
	
	meanIso => Mean isotopic ratio of precipitation
	amplIso => Amplitude of the sine wave of the precipitation isotopic time series
	years => Number of years for which time series is to be obtained
	offset => To take into account that January isotopic ratio is the most negative and July isotopic ratio is the least negative
	
	Returns a numpy array containing the time series of precipitation isotopic ratio
	'''
	
	numb_of_days = int(365 * years)
	dayNumb = np.linspace(1, numb_of_days, numb_of_days) # Corresponding to the number of days in a year
	prcpIso = amplIso * np.sin(2 * np.pi * dayNumb / 365. + offset) + meanIso + np.random.normal(loc=0, scale=0.2 * amplIso, size=numb_of_days)
	return (prcpIso)

def HydroMix(source1, source2, mixture, LAMBDA_prior, ErrorSTD_prior, NUMBER_ITERATIONS):
	'''
	Contains the core of HydroMix
	
	source1 => Concentration data of source 1
	source2 => Concentration data of source 2
	mixture => Concentration data of the mixture made up of a linear combination of sources 1 and 2
	LAMBDA_prior => prior distribution of the proportion of source1 in the mixture
	ErrorSTD_prior => prior distribution of the error standard deviation
	NUMBER_ITERATIONS => Number of times HydroMix has be run

	Returns a tuple containing log likelihood values, lambda values and error standard deviations for all the model runs
	
	'''

	LIKELIHOOD = []
	for i in range(0, NUMBER_ITERATIONS, 1):
		# For displaying purposes
		if (i%100 == 0):
			print ("Iteration number:" + str(i+1))
	
		# Computing likelihood
		temp_likelihood = []
		for temp_mix in mixture:	
			temp_value = 0.
			temp_residual = []
			for temp_s1 in source1:
				for temp_s2 in source2:
					temp_estimated_mix = LAMBDA_prior[i]*1.*temp_s1 + (1-LAMBDA_prior[i])*temp_s2*1.
					
					# Likelihood computation
					temp_value += (-1 * (temp_estimated_mix - temp_mix)**2) / (2 * (ErrorSTD_prior[i]**2))
					temp_value -= (0.5* np.log(2 * np.pi * ErrorSTD_prior[i] * ErrorSTD_prior[i]))
					temp_residual.append(temp_estimated_mix-temp_mix)
			
			temp_likelihood.append(temp_value)
		LIKELIHOOD.append(np.sum(temp_likelihood))
		
	# Sorting the parameter sets by best runs (highest likelihood values)
	zipped = sorted(zip(LIKELIHOOD, LAMBDA_prior, ErrorSTD_prior), reverse = True)
	LIKELIHOOD, LAMBDA_prior, ErrorSTD_prior = zip(*zipped)
	
	return ((LIKELIHOOD, LAMBDA_prior, ErrorSTD_prior))

# This is a csv writer
import csv
def csv_writer(data, path):
	with open(path, "wb") as csv_file:
		writer = csv.writer(csv_file, delimiter=',')
		for line in data:
			writer.writerow(line)

#################################################################################################################################
# %% Main variables

# Percentage of rain and snow recharging groundwater
RAIN_EFF, SNOW_EFF = 1., 1.
K_Q = 0.1 # For Q = KS
DEGDAYFACT = 2.5 # Degree day factor in mm/C/day
TMELT = 0. # Temperature at which snowpack starts melting
CONST_STORAGE, INIT_SNOW = 1000., 0. # Constant storage in groundwater which is not fed to streamflow and initial snow height in catchment
varLSnowTemp, varHSnowTemp = -1., 1. # Threshold temperature for bounding thresholding method for snowfall estimation

# Parameters for generation of air temperature, precipitation and precipitation isotopic ratio
YEARS = 100 # Number of years for which simulation is carried out
AirTempMean, AirTempAmpl = 4., 8.
PrecipEventNumb, PrecipMean = 30, 1000.
PrecipIsoMean, PrecipIsoAmpl = -80., 40.

# To save records into a csv file
FINAL_LIS = [["Day number", "Rainfall (mm)", "Snowfall (mm)", "Storage (mm)", "Snowmelt (mm)", "Snowheight (mm)", "Q from GW (mm)",
	  "Rain recharge (mm)", "Snow recharge (mm)", "Precip isotopic ratio", "Snowpack isotopic ratio", "Storage isotopic ratio"]]
snowMeltIso = [] # For storing snowmelt isotopic ratios
OUTPUTPATH = "OutputFiles/GW_conceptual/"

# Mixing model parameters
NUMBER_ITERATIONS = 1000
LAMBDA_RANGE = [0., 1.] # LAMBDA values imply the fraction of snow in groundwater
# Number of best simulations using which lambda is computed
BEST_SIM_PER = 5. # In percentage

LAST_YEARS = 1 # Number of years at the end of the timeseries from which isotopic data is sampled

#################################################################################################################################
# %%
np.random.seed(15544) # Setting up a common seed number
random.seed(55452) # Setting up random seed for the random function

# %% Time series of air temperature, precipitation and isotopic ratio of precipitation
airTemp = airTemp_Gen(meanTemp = AirTempMean, amplTemp = AirTempAmpl, years = YEARS)
prcp = PoissonPrcp(NumbEvent = PrecipEventNumb, meanPrcp = PrecipMean, years = YEARS)
prcpIso = prcp_Iso(meanIso = PrecipIsoMean, amplIso = PrecipIsoAmpl, years = YEARS)
dayNumb = np.linspace(1, 365*YEARS, 365*YEARS) # Day numbering

# %% Running the hydrologic model
storage = [CONST_STORAGE] # Catchment storage
snow_height = [INIT_SNOW] # Snow height
snow_melt = [0.] # Snow melt
storage_C = [0.] # Isotopic ratio in catchment storage (groundwater)
snow_C = [0.] # Isotopic ratio in snowpack
rain_lis, snow_lis = [0.], [0.] # Daily rainfall and snowfall amount
recharge_rain_lis, recharge_snow_lis = [0.], [0.] # Daily groundwater recharge values from rainfall and snowmelt
Q_flow = [0.] # Daily flow into the stream from groundwater

for index in range(len(dayNumb)):
	
	# Estimating amount of rain and snow
	if (airTemp[index] <= varLSnowTemp):
		rain_frac = 0.
	elif (airTemp[index] >= varHSnowTemp):
		rain_frac = 1.
	else:
		rain_frac = (airTemp[index] - varLSnowTemp) * 1. / (varHSnowTemp - varLSnowTemp)
	rain_mag = rain_frac * prcp[index] *1.
	snow_mag = (1-rain_frac) * prcp[index] *1.
	rain_lis.append(rain_mag)
	snow_lis.append(snow_mag)
	
	# Computing snow height and snowmelt volume
	if (airTemp[index] >= TMELT):
		melt = DEGDAYFACT * (airTemp[index] - TMELT)
		if (melt > snow_height[-1]): # Applying limit of available storage
			melt = snow_height[-1] + snow_lis[-1]
	else:
		melt = 0.
	snow_height.append(snow_height[-1] + snow_lis[-1] - melt)
	snow_melt.append(melt)
	
	# Isotopic ratio of snowpack
	if (snow_height[-1] == 0.): # No more snowpack left
		newC = 0.
	else:
		newC = (snow_height[-2]*snow_C[-1] + snow_lis[-1]*prcpIso[index] - melt*snow_C[-1]) * 1. / (snow_height[-2] + snow_lis[-1] - melt)
		snowMeltIso.append(newC)
	snow_C.append(newC)
	
	# Computing groundwater recharge from rain and snow
	recharge_rain, recharge_snow = rain_lis[-1] * RAIN_EFF * 1., snow_melt[-1] * SNOW_EFF * 1.
	recharge_rain_lis.append(recharge_rain)
	recharge_snow_lis.append(recharge_snow)
		
	recharge = recharge_rain + recharge_snow
	storage.append(storage[-1] + recharge)
	Q = K_Q * (storage[-1] - CONST_STORAGE) # Flow into streamflow from groundwater
	storage[-1] -= Q # Updating groundwater volume
	Q_flow.append(Q)
	
	# Isotopic ratio of groundwater storage
	if (round(storage[-1], 2) == 0.):
		newC = 0.
	else:
		newC = (storage[-2]*storage_C[-1] + rain_lis[-1]*RAIN_EFF*prcpIso[index] + snow_melt[-1]*SNOW_EFF*snow_C[-2] - Q*storage_C[-1]) * 1. / (storage[-2] + rain_lis[-1]*RAIN_EFF + snow_melt[-1]*SNOW_EFF - Q)
	storage_C.append(newC)
	
	FINAL_LIS.append([dayNumb[index], rain_lis[-1], snow_lis[-1], storage[-1], snow_melt[-1], snow_height[-1], Q_flow[-1],
				   recharge_rain_lis[-1], recharge_snow_lis[-1], prcpIso[index], snow_C[-1], storage_C[-1]])
	
# %% Saving the simulated data into a csv file

path = OUTPUTPATH + "RAIN_" + str(RAIN_EFF) + "_SNOW_" + str(SNOW_EFF) + ".csv"
csv_writer(FINAL_LIS, path)

print (path)
#################################################################################################################################
# %% Plotting boxplots of snowfall and snowmelt

snowFallIso = prcpIso[(prcp > 0.) & (airTemp < (0.5 * (varLSnowTemp + varHSnowTemp)))] # Snowfall isotopic ratios

matplotlib.rc('xtick', labelsize=14)
plt.figure(figsize=(4, 6))
plt.boxplot([snowFallIso, snowMeltIso], labels=["Snowfall", "Snowmelt"])
plt.ylabel(r'$\delta^{2}$' + 'H ' + u'\u2030' + ' (VSMOW)', fontsize=14)
plt.tight_layout()
path = OUTPUTPATH + "Figures/SF_SM_boxplot.jpeg"
plt.savefig(path, dpi=300)
plt.close()
print (path)