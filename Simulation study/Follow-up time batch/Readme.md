# Run SAS simulation

## 1. Set github path correctly in the following files:

PoissonF\Simulation study\Follow-up time batch\FIRTHPoisson_simulatoR_FUT_Rev3.r
PoissonF\Simulation study\Follow-up time batch\run_simu_FUT_N60.sas
PoissonF\Simulation study\Follow-up time batch\run_simu_FUT_N200.sas

## 2. Run R program to generate data

Run PoissonF\Simulation study\Follow-up time batch\FIRTHPoisson_simulatoR_FUT_Rev3.r. 

Per scenario defined by number of covariates, number of beta (1-5), and sample size it will create 1 big data set in csv format containing the simulated data sets 1-10000.

## 3. Run SAS program to analyze simulated data sets

In batch mode (right-click in explorer),
Run PoissonF\Simulation study\Follow-up time batch\run_simu_FUT_N60.sas
Run PoissonF\Simulation study\Follow-up time batch\run_simu_FUT_N200.sas

It will generate .sas7bdat files with the results from PROC GENMOD for Firth, FLAC and ML methods.

## 4. Run SAS file to generate and export summaries

Run PoissonF\Simulation study\Follow-up time batch\export parms FUT.sas 
to generate summary statistics like 
bias, MSE, RMSE, RMSE*sqrt(N), lower coverage, upper coverage, power. 
These will be saved as .sas7bdat as well as summary_<Method>_FUT.csv files.


