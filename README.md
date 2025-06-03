# songbird-for-R
Reimplementation of the songbird algorithm by Morton et. al. (2019) in R with additional functionalities.

Original publication:\
https://doi.org/10.1038/s41467-019-10656-5

You can install this package via: \
devtools::install_github("RobertasTupikas/songbird-for-R")

You will need to create a conda environment via Anaconda prompt. The environment has to contain the following: \
tensorflow 2.14.0 \
tensorflow-probability 0.23.0 \
numpy 1.26.4 

You will have to activate the environment in R via use_condaenv("environment", required = TRUE)

An example workflow is provided in: \
songtest.R

The data used is from the teeth brushing examples in the original publication by Morton et. al. (2019)

Here it's: \
physeq.RDS

Workflow:
1. Run the regression model to obtain beta coefficients. Input either a phyloseq object or OTU table + metadata, while specifying which covariate(s) to model in the formula. Also specify model training parameters. You may run the regression 2 times while changing the first OTU in the OTU table, as it is removed from the analysis. The beta coefficients do not change based on this reference OTU, but it is worth checking whether you have not removed a potentially differential OTU. 
2. Make a rank plot, which ranks OTUs based on their beta coefficients. The ones with the highest absolute beta coefficients (on the sides of the plot) are changing the most, while those with the lowest absolute beta coefficients are the most stable (middle of the plot). You may also highlight a OTU of interest to see where it stands in the ranking. Use the plot to identify potentially differential OTUs and a reference OTU 
3. Select a stable OTU as the reference in the log-ratio calculation based on its stability and prevalence in samples. 
4. Using BMDD, model the proportions of all OTUs, efectively imputing zeros. Use the proportions to calculate the log ratios of all OTUs to a specified reference and insert a covariate column that shows which group they belong to (discrete variable). 
5. Run a statistical test (T-test or Wilcoxon test) to compare log ratios of a particular OTU between conditions. 
