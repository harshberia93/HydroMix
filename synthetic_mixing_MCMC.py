#################################################################################################################################
'''
Created on: 27th September 2019 (Harsh Beria) 
What it does?

A mixing model to estimate source contribution in a two-source linear system using an average likelihood approach.
Error in estimation is assumed to be normally distributed with zero mean and constant standard deviation (which is computed from data)

Model parameters:
NUMBER_ITERATIONS => Number of model runs
LAMBDA_params => Fraction of source 1 in the mixture
LIKELIHOOD_std_params => Fixed value

SourceFiles used: NONE

OutputFiles processed:
OutputFiles/Synthetic_exp/
	
Figures made:
'''
#################################################################################################################################
# %% Imports

from __future__ import division

import numpy as np
from matplotlib import pyplot as plt

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


def HydroMix(source1, source2, mixture, stdDev, param_init, param_limit, nb_iter, RandomWalkStep=5):
	'''
	Contains the core of HydroMix
	
	source1 => Concentration data of source 1
	source2 => Concentration data of source 2
	mixture => Concentration data of the mixture made up of a linear combination of sources 1 and 2
	stdDev => Standard deviation of the mixture (specified apriori and not calibrated)
	param_init => Initial value of the list of parameters ([lambdaValue, stdValue])
	param_limit => Lower and upper limits of the list of parameters ([[lambdaLowVal, lambdaHighVal]])
	nb_iter => Number of MCMC runs
	RandomWalkStep => In percentage to define the maximum extent of jump from the current parameter state

	Returns a tuple containing [[LOG LIKELIHOOD VALUES], [PARAMETER VALUES], [RESIDUAL VALUES]]
	
	'''

	logLikelihoodLis, paramLis, residual_list = [], [param_init], []
	std_param = stdDev
	n_accepts = 0
	for i in range(nb_iter):
		
		# Compute new parameter set
		updatedParam = Random_walk(paramLis[-1], param_limit, RandomWalkStep)
	
		# Log likelihood computation
		lambda_param, temp_residual = updatedParam[0], []
		temp_loglikelihood = []
		for temp_mix in mixture:	
			temp_value = 0.
			for temp_s1 in source1:
				for temp_s2 in source2:
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
			print ( "Iteration number:" + str(i+1) + ", Acceptance: "+ str(n_accepts) + ", Acceptance %: " + 
																			  str(n_accepts/(i+1)) )
		
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

# HydroMix model parameters (in the given range (inclusive))
NUMBER_ITERATIONS = 1000
LAMBDA_RANGE = [0., 1.]
JUMP_PERCENTAGE = 5. # In percentage (JUMP_PERCENTAGE/2 in both directions)

# Parameters to generate synthetic time series
# Lower and upper bounds of rain and snow isotopic ratio (in H2)
N_SAMPLES = 100 # Number of samples
S1_MEAN, S1_STD, N_S1_SAMPLES = 10., .5, N_SAMPLES
S2_MEAN, S2_STD, N_S2_SAMPLES = 20., .5, N_SAMPLES
N_MIX_SAMPLES = N_SAMPLES

OUTPUTPATH = "OutputFiles/Synthetic_exp/Normal_distribution_0.5Std_" + str(NUMBER_ITERATIONS) + "iter_"
OUTPUTPATH += str(N_SAMPLES) + "sample/"

# %% Setting up a random seed
np.random.seed(1) # Setting up a common seed number

# %% Mixing for lambda ranging from 0.05 to 0.95
LAMBDA = 0.05 # Ratio of source 1 in the mixture
scatterplot_orig_lambda, scatterplot_sim_lambda = [], [] # For displaying in the scatterplot

while (LAMBDA <= 0.96):

	# Computing GW mean and variance
	MIX_MEAN = LAMBDA * S1_MEAN + (1-LAMBDA) * S2_MEAN
	MIX_STD = np.sqrt((LAMBDA*S1_STD)**2 + ((1-LAMBDA)*S2_STD)**2)
	
	# Generating synthetic values for the sources and the mixture
	S1_val = np.random.normal(S1_MEAN, S1_STD, N_S1_SAMPLES)
	S2_val = np.random.normal(S2_MEAN, S2_STD, N_S2_SAMPLES)
	MIX_val = np.random.normal(MIX_MEAN, MIX_STD, N_MIX_SAMPLES)
	
	# Saving it into a csv file
	path = OUTPUTPATH + "input_" + str(LAMBDA) + ".csv"
	output_lis = [["Source 1", "Source 2", "Mixture"]]
	for temp_index in range(len(S1_val)):
		output_lis.append([S1_val[temp_index], S2_val[temp_index], MIX_val[temp_index]])
	csv_writer(output_lis, path)
	print (path)
	
	# List of initial parameter values 
	initParam = [np.random.uniform(LAMBDA_RANGE[0], LAMBDA_RANGE[1])]
	paramLimit = [LAMBDA_RANGE] # Lower and upper limits of the model parameters
	
	# Running the mixing model
	LOGLIKELIHOOD, PARAM, RESIDUAL = HydroMix(S1_val, S2_val, MIX_val, MIX_STD, initParam, paramLimit, NUMBER_ITERATIONS, JUMP_PERCENTAGE)

	mixingRatioLis = [i[0] for i in PARAM]

	# Writing into a csv file
	path = OUTPUTPATH + "output_" + str(LAMBDA) + "_MCMC.csv"
	output_lis = [["Lambda value", "Log likelihood", "Error std", "Residual"]] 
	for index in range(0, len(LOGLIKELIHOOD)):
		output_lis.append([ round(mixingRatioLis[index], 4), round(LOGLIKELIHOOD[index], 4), round(MIX_STD, 4), 
					 round(RESIDUAL[index], 4) ])
	csv_writer(output_lis, path)
	print (path)
	
	
#	scatterplot_orig_lambda.append(LAMBDA)
#	scatterplot_sim_lambda.append(np.median(LAMBDA_params[0:int(BEST_SIM_PER * NUMBER_ITERATIONS / 100.)]))
#	print (LAMBDA, np.median(LAMBDA_params[0:int(BEST_SIM_PER * NUMBER_ITERATIONS / 100.)]))
	LAMBDA += 0.05
	
	
#################################################################################################################################
## %% Plotting scatterplot
#
#plt.figure()
#plt.scatter(scatterplot_orig_lambda, scatterplot_sim_lambda)
#plt.plot([0., 1.], [0., 1.], 'k--', color = 'r')
#plt.xlabel("Original lambda")
#plt.ylabel("Simulated labmda")
#plt.xlim(0., 1.)
#plt.ylim(0., 1.)
#plt.tight_layout()
#path = OUTPUTPATH + "Figures/scatterplot.jpeg"
#plt.savefig(path)
#
#plt.close()

#################################################################################################################################
