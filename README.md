# HydroMix

HydroMix is a Bayesian mixing model used to quantify the proportion of different sources in a given mixture. It is implemented in Python version 2.7.

HydroMix was originally implemented to solve mixing problems in hydrology, such as 2-component hydrograph separation, resolving sources of tree water uptake, role of snowmelt in groundwater recharge, etc. HydroMix formulates the linear mixing problem in a Bayesian inference framework, where the model error is parameterized instead of the source distributions. This makes the framework suitable for problems where the source sample size is small. Also, it is easy to introduce additional model parameters that can transform the source composition and infer them jointly. For problems with large sample size, HydroMix becomes computationally very expensive. In such cases, parameterizing the source distributions, and then inferring the probability density function of the mixing ratio using a probabilistic programming language (eg: Stan) is more effective.

A detailed manuscript describing the model can be seen at:
https://www.geosci-model-dev-discuss.net/gmd-2019-69/

Code for the following case studies have been included in this document:

1. Synthetic case study
	This has been explained in Section 3.1 of the manuscript

	A MCMC sampler has been implemented to solve the mixing problem with two normally distributed random variables

	Script: synthetic_mixing_MCMC.py
		Generates a given number of samples (specified by user) and generates samples for two random variables with distinct means and variances (specified by user). Generates sample from the mixture distribution for a given mixing ratio (specified by user) and then uses the mixture and source samples to infer back the mixing ratio
	Example runs for 100 samples and standard deviations of 0.5 and 5 are stored in 
	OutputFiles/Synthetic_exp/Normal_distribution_0.5Std_1000iter_100sample/
	OutputFiles/Synthetic_exp/Normal_distribution_5Std_1000iter_100sample/

2. Conceptual groundwater model case study
	This has been explained in Section 3.2 of the manuscript
	
	Model run with a given rain and snow recharge efficiency can be made with the script: GW_conceptual.py. This script was run 100 times by varying rain and snow efficiency between 0.1 to 1 in steps of 0.1 and the results are saved in OutputFiles/GW_conceptual/ folder

	An importance sampling and a MCMC sampler has been implemented to solve the mixing problem

	Script for importance sampling: GW_mixing_HydroMix.py
	The script allows using snowfall or snowmelt to do mixing. It also allows weighting samples.
	OutputFiles produced in folder:
		OutputFiles/GW_conceptual/Rainfall_Snowfall_mixing_last_2Yr/
		OutputFiles/GW_conceptual/Rainfall_Snowfall_mixing_last_2Yr_weighted/
		OutputFiles/GW_conceptual/Rainfall_Snowmelt_mixing_last_2Yr/
		OutputFiles/GW_conceptual/Rainfall_Snowmelt_mixing_last_2Yr_weighted/

	Script for MCMC sampling: GW_mixing_HydroMix_MCMC.py
	The script allows using snowfall or snowmelt to do mixing. It also allows weighting samples.
	OutputFiles produced in folder:
		OutputFiles/GW_conceptual/Rainfall_Snowfall_mixing_last_2Yr_MCMC/
		OutputFiles/GW_conceptual/Rainfall_Snowfall_mixing_last_2Yr_weighted_MCMC/
		OutputFiles/GW_conceptual/Rainfall_Snowmelt_mixing_last_2Yr_MCMC/
		OutputFiles/GW_conceptual/Rainfall_Snowmelt_mixing_last_2Yr_weighted_MCMC/

3. Real case study (Vallon de Nant):
	This has been explained in Section 3.4 of the manuscript
	Uses the following source files:
		SourceFiles/Rain_SP_GW_VdN.xlsx
			DESCRIPTION: The source file provides data for isotopic ratio in snowpack, rainfall and groundwater collected in Vallon de Nant (an alpine catchment in Western Switzerland). The timestamp is provided in the excel sheet.
	
	An importance sampling and a MCMC sampler has been implemented to solve the mixing problem

	Script for importance sampling: VdN_mixing.py
	OutputFiles produced:
		OutputFiles/VdN/H2_results.csv
		OutputFiles/VdN/O18_results.csv
		OutputFiles/VdN/posterior.jpeg

	Script for MCMC sampling: VdN_mixing_MCMC.py
	OutputFiles produced:
		OutputFiles/VdN/H2_results_MCMC_1000.csv
		OutputFiles/VdN/O18_results_MCMC_1000.csv

4. Real case study (with additional model parameters)
	This has been explained in Section 3.5 of the manuscript
	Uses the following source files:
		SourceFiles/Rain_SP_GW_VdN.xlsx (described in point 3)
		SourceFiles/Hypsometric_curve_data.xlsx
			DESCRIPTION: Hypsometric curve of Vallon de Nant discretized at 100 m resolution

		An importance sampling and a MCMC sampler has been implemented to solve the mixing problem using an additional model parameter, lapse rate

		Script for importance sampling: VdN_mixing_additional_params.py
		OutputFiles produced:
			OutputFiles/VdN/Lapse_rate/H2_results.csv
			OutputFiles/VdN/Lapse_rate/O18_results.csv
			OutputFiles/VdN/Lapse_rate/posterior.jpeg
			OutputFiles/VdN/Lapse_rate/lapse_rate_posterior.jpeg

		Script for MCMC sampling: VdN_mixing_additional_params_MCMC.py
		OutputFiles produced:
			OutputFiles/VdN/Lapse_rate/H2_results_MCMC_10000.csv
			OutputFiles/VdN/Lapse_rate/O18_results_MCMC_10000.csv
