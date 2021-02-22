%include "D:\1Projects\2018-08 NegBinom FC\Github\PoissonF\Simulation study\Sas macro\flacpoisson.sas";
%include "D:\1Projects\2018-08 NegBinom FC\Github\PoissonF\Simulation study\Sas macro\dataugpoisson.sas";

PROC IMPORT OUT= WORK.dataset 
            DATAFILE= "D:\1Medizin\K
uchler Zahnklinik\Prediction of complications\dataSet.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data dataset_tmp;
set dataset;
length patno_n 8;
Comp_n=(compinf="yes") + (compbleed="yes") + (compnerve="yes") + (comphema="yes") + (Compfail="yes");
compinf_n=(compinf="yes");
comphema_n=(comphema="yes");
if facttime ne . then do;
	FactTime_O1v3 = (FactTime <=1);
	FactTime_O2v3 = (facttime = 2);
end;
else do;
	FactTime_O1v3 = .;
	FactTime_O2v3 = .;
end;
Patno_n=patno;
factage_2=(factage-55)**2;
factage_cent=factage-55;
factage_cent_dec=floor(factage_cent/10);
run;
data dataset_t1;
set dataset_tmp;
by patno_n;
if first.patno_n;
run;

proc means data=dataset_tmp noprint;
by patno_n;
var comp_n Compinf_n comphema_n factind_D1 factind_d2 factind_d3 factind_d4 factreg_d1 factreg_d2 factreg_d3 factreg_d4; 
output out=dataset_agg sum=comp_n compinf_n comphema_n  factind_D1 factind_d2 factind_d3 factind_d4 factreg_d1 factreg_d2 factreg_d3 factreg_d4 n=nimpl;
run;

data dataset_m;
merge dataset_tmp(drop=factind_d1 factind_d2 factind_d3 factind_d4 factreg_d1 factreg_d2 factreg_d3 factreg_d4) dataset_agg;
by patno_n;
if first.patno_n;
log_nimpl=log(nimpl);
log_unif =log(ranuni(538401)*2);
run;

proc univariate data=dataset_m;
var nimpl;
histogram;
run;

proc freq data=dataset_m;
tables factage_cent_dec;
tables factSmoke;
tables FactDm;
tables nimpl;
tables comphema_n;
run;

proc means min max data=dataset_m;
var factage;
class factage_cent_dec;
run;


proc genmod data=dataset_m;
model comphema_n = FactAge_cent_Dec FactSmoke_O2 FactSmoke_O3 FactAugm_D FactOsteo_D FactInd_D1 FactInd_D2 FactInd_D3 FactInd_d4 FactReg_D1 FactReg_D2 FactReg_D3 FactReg_D4 FactTime_O2v3 FactTime_O1v3 FactDM_D FactAnti_d
	/dist=poisson offset=log_nimpl;
run;

proc genmod data=dataset_m;
model comphema_n =  FactAge_cent_dec FactSmoke_O2  FactSmoke_O3 FactDM_D
	/dist=poisson offset=log_nimpl;
run;


proc genmod data=dataset_m;
model comphema_n =  FactAge_cent_dec FactSmoke_O2  FactSmoke_O3 FactDM_D
	/dist=poisson offset=log_nimpl lrci;
run;

*** data augmentation prior (prior interval for IRR e.g. 1/100, 100);


%dataugpoisson(varlist=a b c d, priorinterval=1.05 1.1 5 10, printprior=1);

%dataugpoisson(data=dataset_m, y=comphema_n, offset=log_nimpl, varlist=FactAge_cent_dec FactSmoke_O2  FactSmoke_O3 FactDM_D, priorinterval=100  1000, S=25);

%dataugpoisson(data=dataset_m, y=comphema_n, offset=log_nimpl, varlist=FactAge_cent_dec FactSmoke_O2  FactSmoke_O3 FactDM_D, priorinterval=100 1000);

%dataugpoisson(data=dataset_m, y=comphema_n, offset=log_nimpl, varlist=FactAge_cent_dec FactSmoke_O2  FactSmoke_O3 FactDM_D, priorinterval=5 50);

proc logxact data=dataset_m;
model comphema_n =  FactAge_cent_dec FactSmoke_O2 FactSmoke_O3  FactDM_D
	/link=poisson ;
es /ex intercept estimatefile=est_logxact  FactAge_cent_dec FactSmoke_O2 FactSmoke_O3  FactDM_D;
rate nimpl;
run;



proc logxact data=dataset_m;
model comphema_n =  FactAge_cent_dec FactSmoke_O2 FactSmoke_O3  FactDM_D
	/link=poisson ;
es /as estimatefile=est_lxpmle pmle plci   factSmoke_o2 factSmoke_o3 factage_cent_dec factDM_d ;
rate nimpl;
run;


%flacpoisson(data=dataset_m, y=comphema_n, 
varlist= FactAge_cent_dec FactSmoke_O2 FactSmoke_O3 FactDM_D,
offset=log_nimpl);
/*
proc univariate data=dataset_m;
var log_nimpl log_unif;
histogram;
run;
*/


* compare predicted values;

proc means data=_firthpredictions sum mean;
var comphema_n _firthpred;
run;

proc means data=_flacpredictions sum mean;
var comphema_n _flacpred;
run;
