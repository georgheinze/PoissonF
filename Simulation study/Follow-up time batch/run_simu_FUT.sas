*** run simulation for varying follow-up times;




PROC IMPORT OUT= WORK.simdata 
	            DATAFILE= "E:\1Projects\2018-08 NegBinom FC\Simu Data\Alldata FUT\Simdata_FUT_&ncov._&nbeta._&n..csv" 
	            DBMS=CSV REPLACE;
	     GETNAMES=YES;
	     DATAROW=2; 
	RUN;



data simdata;
set simdata;
logFU=log(FUT);
run;

* sanity check: analyze all data sets together;
*%flacpoisson(data=simdata, y=y, varlist=x1 x2, by=, maxiter=50, offset=logFU, print=1, nonotes=0);

data simdata10;
set simdata;
if isim<=10000;
run;

options nonotes;
%flacpoisson(data=simdata10, y=y, varlist=x1 x2, by=isim, maxiter=100, offset=logFU, print=0, odsselect=all, nonotes=1);
options notes;

data poislrci.FLAC_FUT_&ncov._&nbeta._&n.;
set _FLACparms;
run;

data poislrci.FIRTH_FUT_&ncov._&nbeta._&n.;
set _FIRTHparms;
run;


ods select none;
proc genmod data=simdata10;
ods output parameterestimates=poislrci.ML_FUT_&ncov._&nbeta._&n.;
model y=x1 x2/dist=poisson offset=logFU;
by isim;
run;
ods select all;

/*
ods select none;
proc genmod data=simdata10;
ods output exactparmest=poislrci.EX_FUT_&ncov._&nbeta._&n.;
model y=x1 x2/dist=poisson offset=logFU;
exact x1/estimate cltype=exact;       ********** offset not available with exact model!
by isim;
run;
ods select all;


ods select none;
proc genmod data=simdata10;
ods output exactparmest=poislrci.MIDP_FUT_&ncov._&nbeta._&n.;
model y=x1 x2/dist=poisson offset=logFU;
exact x1/estimate cltype=midp;********** offset not available with exact model!
by isim;
run;
ods select all;
*/



