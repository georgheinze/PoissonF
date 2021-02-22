# Run SAS simulation

## Caution

The SAS simulation to be executable requires a SAS installation and folder structure similar to the one used in the original simulation study. We used SAS 9.4TS1M7 with Cytel's PROC LOGXACT. The simulation study was carried out in portions and required several weeks to complete.

## 1. Set paths correctly in the following files:

* PoissonF\Simulation study\Follow-up time batch\FIRTHPoisson_simulatoR_FUT_Rev4.r
* PoissonF\Simulation study\Follow-up time batch\run_simu_EX_ncov2.sas
* PoissonF\Simulation study\Follow-up time batch\run_simu_EX_ncov5.sas
* PoissonF\Simulation study\Follow-up time batch\run_simu_EX_ncov10.sas

*Githubpath* is the path to where the code available in this repository is cloned to and to the place where results files will be stored (in csv format),
*Localpath* is the path to a folder where SAS results files will be stored (in sas7bdat format). *Don't push the resulting csv files back to github as the required space is excessive!*

## 2. Run R program to generate data

Run 'PoissonF\Simulation study\Follow-up time batch\FIRTHPoisson_simulatoR_FUT_Rev4.r'. This code references to the 'simdata' R package which is available at https://github.com/matherealize/simdata.

Each scenario is defined by an approximate event rate of 0.1,  number of covariates (2, 5 or 10), beta1 (9 equidistant values from -log(16) to log(16)), and sample size (such that events per variable ratios of approximately 3, 5, and 10 result). For each of these 81 scenarios, the R program will create 1 big data set in csv format containing the simulated data sets 1-10000. It will also create a file with the intercepts (and beta1, beta2) iteratively found to create an event rate of approximately 0.10.

The data sets are saved in folder 'Data'. 

## 3. Run SAS program to analyze simulated data sets

In batch mode (right-click in explorer),

* Run PoissonF\Simulation study\Follow-up time batch\run_simu_EX_ncov2.bat
* Run PoissonF\Simulation study\Follow-up time batch\run_simu_EX_ncov5.bat
* Run PoissonF\Simulation study\Follow-up time batch\run_simu_EX_ncov10.bat

These files will call SAS programs in batch mode and with Cytel PROCs configuration.

Upon completing, they will generate .sas7bdat and .csv files with the results from PROC GENMOD for Firth, FLAC and ML methods.

These files are saved in folder 'Results' of the Githubpath and Localpath folders.

A pooled (maximum likelihood) analysis of all data sets simultaneously serves to estimate asymptotic bias which could be caused by the restrictions
imposed in generating the data sets (minimum event rule/minimum frequency of each level of a binary covariate rule).
The pooled analysis is saved in POOLED_ML_FUT.csv and POOLED_ML_FUT.sas7bdat.

## 4. Run SAS file to generate and export summaries

Run 'PoissonF\Simulation study\Follow-up time batch\extract_results.sas' 
to generate summary statistics like 
bias, MSE, RMSE, RMSE*sqrt(N), lower coverage, upper coverage, power. 
These will be saved as .sas7bdat (on localpath) as well as .csv files (on Githubpath).

All summaries are saved in folder 'Summaries'. Based on the two summary files 'PoissonF_sep_nonconv.csv' and 'PoissonFResultsTable.csv', the R-markdown file 'SimulationReport.Rmd' generates the full summary of the simulation study.
