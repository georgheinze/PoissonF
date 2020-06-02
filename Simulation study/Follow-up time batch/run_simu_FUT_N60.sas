%let Githubpath = %str(E:\1Projects\2018-08 NegBinom FC\Github\);


%include "&githubpath.PoissonF\Simulation study\Sas macro\flacpoisson.sas";

libname poislrci "&githubpath.PoissonF\Simulation study\Follow-up time batch\Results";


%let ncov=2;
%let n=60;
%let nbeta=1;

%include "&githubpath.PoissonF\Simulation study\Follow-up time batch\run_simu_fut.sas";

%let ncov=2;
%let n=60;
%let nbeta=2;

%include "&githubpath.PoissonF\Simulation study\Follow-up time batch\run_simu_fut.sas";

%let ncov=2;
%let n=60;
%let nbeta=3;

%include "&githubpath.PoissonF\Simulation study\Follow-up time batch\run_simu_fut.sas";

%let ncov=2;
%let n=60;
%let nbeta=4;

%include "&githubpath.PoissonF\Simulation study\Follow-up time batch\run_simu_fut.sas";


%let ncov=2;
%let n=60;
%let nbeta=5;

%include "&githubpath.PoissonF\Simulation study\Follow-up time batch\run_simu_fut.sas";

