*** run simulation for varying follow-up times;
*** needs xlist as the set of variables, e.g. xlist = x1 x2 x3 x4 x5;

%include "&githubpath.PoissonF\Simulation study\Sas macro\dataugpoisson.sas";

%include "&githubpath.PoissonF\Simulation study\Sas macro\datauginportion.sas";

%include "&githubpath.PoissonF\Simulation study\Sas macro\flacpoisson.sas";

%include "&githubpath.PoissonF\Simulation study\Sas macro\flacinportion.sas";

%include "&githubpath.PoissonF\Simulation study\Sas macro\mlexinportion.sas";

%include "&githubpath.PoissonF\Simulation study\Sas macro\mlLXinportion.sas";

%include "&githubpath.PoissonF\Simulation study\Sas macro\data2r.sas";


* read FUT data set as it also holds y simulated without FUT;
PROC IMPORT OUT= WORK.simdata10
	            DATAFILE= "&localpath.\PoissonF\Simulation study\Follow-up time batch\Data\Simdata_FUT_&ncov._&nbeta._&n..csv" 
	            DBMS=CSV REPLACE;
	     GETNAMES=YES;
	     DATAROW=2; 
	RUN;



data simdata10;
set simdata10;
logFUT=log(FUT);
if isim<=10000;
run;
* sanity check: analyze all data sets together;

ods select none;
proc genmod data=simdata10;
ods output parameterestimates=poislrci.POOLED_ML_FUT_&ncov._&nbeta._&n.;
model y_FUT=&xlist /dist=poisson offset=logFUT;
*by isim;
run;
ods select all;

%data2r(data=poislrci.POOLED_ML_FUT_&ncov._&nbeta._&n., dir=%str(&localpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=POOLED_ML_FUT_&ncov._&nbeta._&n..csv);


options notes;
%flacinportion(portions=10, nsim=10000, y_2go=y_FUT, FUT=1, keep=y mu FUT mu_FUT);
%datauginportion(portions=10, nsim=10000, y_2go=y_FUT, FUT=1, keep=mu FUT mu_FUT);

options notes;


%if &do_exact=1 %then %do; 
	%mlexinportion(portions=200, nsim=10000, y_2go=y_FUT, exact=&do_exact, lxpmle=&do_lxpmle, FUT=1);
%end;
%if &ml=1 & &do_exact=0 %then %do;
	%mlexinportion(portions=10, nsim=10000, y_2go=y_FUT, exact=&do_exact, lxpmle=&do_lxpmle, FUT=1);
%end;





