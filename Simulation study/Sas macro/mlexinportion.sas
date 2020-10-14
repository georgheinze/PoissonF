%macro MLexinportion(portions=10, nsim=10000, y_2go=y, exact=0, lxpmle=0, fut=0);


* performs ML and exact computation (exact and midp CI) and Firth PMLE (logxact implementation with PLCI) in portions to avoid memory overflow;
* analyses y as outcome variable and &xlist as independent variables;

*** sas data sets will be written to library poislrci;
*** below, sas data sets will then be exported to csv to githubpath;

data poislrci.ML_&ncov._&nbeta._&n.;
run;

data poislrci.ML_PRED_&ncov._&nbeta._&n.;
run;


%if &lxpmle=1 %then %do;
	data poislrci.LXPMLE_&ncov._&nbeta._&n.;
	run;
%end;

%if &exact=1 %then %do;
	data poislrci.logXACT_&ncov._&nbeta._&n.;
	run;
%end;



%let portionsize=%sysevalf(&nsim / &portions);

%do iportion=1 %to &portions;

	data simdata_portion;
	set simdata10;
	if isim > (&iportion -1) * &portionsize AND isim <= &iportion * &portionsize;
	run;

	ods select none;
	proc genmod data=simdata_portion;
*	ods output parameterestimates=ML_parm_portion exactparmest=EXACT_parm_portion;
	ods output parameterestimates=ML_parm_portion;
	model &y_2go = &xlist /dist=poisson %if &fut=1 %then %do; offset=logFUT %end; ;
*	exact x1 /cltype=exact estimate;
*	exact x1 /cltype=midp estimate;
*	exact x2 /cltype=exact estimate;
*	exact x2 / cltype=midp estimate;
*	exactoptions method=direct;
	output out=ML_Pred_portion predicted=_MLPred;
	by isim;
	run;
	ods select all;


	data poislrci.ML_&ncov._&nbeta._&n.;
	set poislrci.ML_&ncov._&nbeta._&n. ML_parm_portion;
	run;

	data poislrci.ML_PRED_&ncov._&nbeta._&n.;
	set poislrci.ML_PRED_&ncov._&nbeta._&n. ML_Pred_portion;
	run;


	%if &exact=1 %then %do;
		proc logxact data=simdata_portion noprint;
		model &y_2go = &xlist / link=poisson;
		rate FUT;
		es /ex estimatefile=LogXACT_parm_portion x1 x2;
		by isim;
		run;
		data poislrci.logXACT_&ncov._&nbeta._&n.;
		set poislrci.logXACT_&ncov._&nbeta._&n. logXACT_parm_portion;
		run;
	%end;
	
	%if &lxpmle=1 %then %do;
		proc logxact data=simdata_portion noprint;
		model &y_2go = &xlist / link=poisson;
		rate FUT;
		es /as Intercept estimatefile=LXPMLE_parm_portion PMLE PLCI  x1 x2;
		by isim;
		run;


		data poislrci.LXPMLE_&ncov._&nbeta._&n.;
		set poislrci.LXPMLE_&ncov._&nbeta._&n. LXPMLE_parm_portion;
		run;
	%end;


%end;

*** export tables into githubpath;

%data2r(data=poislrci.ML_&ncov._&nbeta._&n., dir=%str(&localpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=ML_FUT_&ncov._&nbeta._&n..csv);
%data2r(data=poislrci.ML_PRED_&ncov._&nbeta._&n., dir=%str(&localpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=ML_FUT_PRED_&ncov._&nbeta._&n..csv);

%if &exact=1 %then %do;
	%data2r(data=poislrci.logXACT_&ncov._&nbeta._&n., dir=%str(&localpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=logXACT_FUT_&ncov._&nbeta._&n..csv);
%end;

%if &lxpmle=1 %then %do;
	%data2r(data=poislrci.LXPMLE_&ncov._&nbeta._&n., dir=%str(&localpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=LXPMLE_FUT_&ncov._&nbeta._&n..csv);
%end;

%mend;


