# HydroMix
HydroMix is a Bayesian mixing model used to quantify the proportion of different sources in a given mixture. It is implemented in Python version 2.7.

HydroMix was originally implemented to solve mixing problems in hydrology, such as 2-component hydrograph separation, resolving sources of tree water uptake, role of snowmelt in groundwater recharge, etc. HydroMix formulates the linear mixing problem in a Bayesian inference framework, where the model error is parameterized instead of the source distributions. This makes the framework suitable for problems where the source sample size is small. Also, it is easy to introduce additional model parameters that can transform the source composition and infer them jointly. For problems with large sample size, HydroMix becomes computationally very expensive. In such cases, parameterizing the source distributions, and then inferring the probability density function of the mixing ratio using a probabilistic programming language (eg: Stan) is more effective.

A detailed manuscript describing the model is currently under review in Geoscientific Model Development.
