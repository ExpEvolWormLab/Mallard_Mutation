# Data and code for analysis for the MS "Variation in mutational (co)variances"
# Please contact François Mallard (mallard@bio.ens.psl.eu) in case you have problems with these files.

The data directory contains:
- the raw transition rates in  Final_merged_data_MA_Lines
- the A6140 G matrix released in the Mallard Stasis manuscript
- a table with the sample size per assay block

The Rcode directory contains all code to reproduce the analysis.

The RData directory contains intermediate files for the analysis, mainly the output from the MCMCglmm models and the null expectations that are computationally demanding.

The Output_files directory contains txt output files showing the different computed M matrices (txt directory) and the auto-correlation plots from the MCMCglmm models.

DETAILS OF THE RDATA FILES IN THE OUTPUT DIRECTORY (4 files):

In detail, the RData files contain the output from the MCMCglmm models that are long to compute. All remaining R code is rapid to run.
 
The MCMCglmm are saved as R lists in 'MCMC_glmm_output.RData' and 'MCMC_glmm_output_SCALED.RData'. Each list contains the following 6 elements:

- Population (chr): the population ID 
- N measurement (int): the number of measurements used to compute the model
- G1_mat, G2_mat, R_mat (num [1:6,1:6]): posterior mode of the 6x6 matrices computed by the MCMCglmm. G1 is the random effect for line ID (the G matrix), G2 is the random effect for block ID and R is the matrix residual variance.
- VCV_Mat (mcmc num [1:XXX,1:108]): the posterior distribution of the three matrices described above. The number of rows might vary depending on the model run.

The RData file "Mmatrix_Null.RData" contains the outputs of the Null distribution for the N2 and PB306 G matrices (N2_G1_null and PB_G1_null respectively). Each variable is a matrix of size [1:6000,1:6] containing 1000 G matrices of size 6x6 concatenated by row. Each of these 1000 "random" G matrix is computed by an independent random shuffling of the line ID label and block ID label using MCMCglmm. Only the posterior mode of the MCMCglmm are kept here.


The RData file "Tensor_processed_N2PB_withA6140.Rdata" contains a list of 1000 random G array as computed in Morrissey and Bonnet (2019).

DETAILS OF PLOTS OF THE MCMC OUTPUTS:

There are two types of plots. 1/ Outputs of the MCMCglmm plots from R (see ?plot.MCMCglmm) are saved as "M_for_Model_pdf_MCMC_XX.pdf" with XX being the genetic background. 2/Autocorrelation plots performed using the acf function in R are saved as "M_for_Model_pdf_MCMC_autocorr_XX.pdf" with XX being the genetic background. 
Plots are produced for the initial matrices computation (auto_corr_plots_MCMCglmm directory) and for the scaled ones (auto_corr_plots_MCMCglmm_scaled directory).


DETAILS OF THE DATA FILES:

- Final_merged_data_MA_Lines is a tab-separated file containing headers (1:261,1:24). Here is a summary description of the columns :
	- exper_data_id: unique identifier of the movie with block ID and date and time stamps
	- pop_label: line ID / concatenation of the background (N2anc and PBanc for the ancestors and N2250 and PB250 for MA lines, then followed by the line number (L0 for the ancestor).
	- data_group_name: block ID
	- data_id: unique identifier of the movie with the line ID (concatenation of lineID,env_label and blockID)
	- date_str: date of the movie (YYYYMMDD)
	- time_str: starting time of the movie (HHMMSS)
	- exper_duration_s: duration of the movie
	- trackL: mean track length
	- n_tracks: number of tracks (total)
	- mean_nb_tracks: number of tracks per unit of time
	- logD,temperature,rel_humidity: log number of tracks per unit of time, temperature and relative humidity at the time of measurement
	- location_label: location name, city where the measurements have been done (values: Lisbon, Paris, Paris2)
	- env_label: environment of culture growth (values: NGM for standard environmental conditions, NaCl for high salt) 
	- T12,T13,T21,T23,T31,T32: transition rates measurements. 1:still,2:forward,3:backward
	- mean_still, mean_forward, mean_backward: mean proportion of time spend in each state


- The "VCV_A6140.RData" file contains two elements:
	- final_merged, a data frame (1833 abs of 23 variables) containing the transition rate measurements. Column names match the data frame detailed above.
	- VCV_mat_A6140, a list containing the MCMCglmm output of the A6140 populations from the Mallard_Stasis manuscript with the same structure described above.

- sample size : table presenting the number of phenotypes lines per block (first column). Number of MA lines per ancestral background are shown in columns 2 and 3. Columns 4 and 5 are binary values (0/1) indicating whether the ancestor was included in the block.

 
DETAILS OF RCODE FILES:

- 01_Mutational_Bias
## This code produce the phenotypic analysis of the MA lines from the 2 founder genotype produce figure 1 of the manuscript

- 02A_Produce_M_matrices
## This code computes the two M-matrices using MCMCglmm - Figure 3, Figure S1 and S2

- 02B_Produce_M_matrices_NULL
## This code produce the Null expectation for the G matrices by shuffling the line ID and block ID columns

- 03_Produce_M_matrices_Scaled
## This code compute the scaled G and M matrices - plot for Figure 4

- 04_MMatrice_plots
## Code that produce the plots for Figure 2 and some exports for tables

- 05_Project_Ms
#Projections and angles between Mmax and Gmax vectors - Figure 3C and 5D

- 06_Tensor_Analysis_withA6140_Scaled
## Code that performs the Tensor analysis and the related plots. Updated from Aguirre et al. 2014 (Supplementary information) and Morrissey and Bonnet (2019).

- functions_tensor
## Functions for the tensor analysis as published by Aguirre et al. 2014 (Supplementary information)