%macro MLLXinportion(portions=10, nsim=10000, y_2go=y);


* performs ML and Firth PMLE (logxact implementation) computation (PLCI) in portions to avoid memory overflow;
* analyses y as outcome variable and &xlist as independent variables;

*** sas data sets will be written to library poislrci;
*** below, sas data sets will then be exported to csv to githubpath;

data poislrci.ML_&ncov._&nbeta._&n.;
run;


data poislrci.LXPMLE_&ncov._&nbeta._&n.;
run;


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
	model &y_2go = &xlist /dist=poisson;
*	exact x1 /cltype=exact estimate;
*	exact x1 /cltype=midp estimate;
*	exact x2 /cltype=exact estimate;
*	exact x2 / cltype=midp estimate;
*	exactoptions method=direct;
	by isim;
	run;
	ods select all;


	proc logxact data=simdata_portion noprint;
	model &y_2go = &xlist / link=poisson;
	es /as Intercept estimatefile=LXPMLE_parm_portion PMLE PLCI  x1 x2;
	by isim;
	run;


	data poislrci.ML_&ncov._&nbeta._&n.;
	set poislrci.ML_&ncov._&nbeta._&n. ML_parm_portion;
	run;


	data poislrci.LXPMLE_&ncov._&nbeta._&n.;
	set poislrci.LXPMLE_&ncov._&nbeta._&n. LXPMLE_parm_portion;
	run;


%end;

*** export tables into githubpath;

%data2r(data=poislrci.LXPMLE_&ncov._&nbeta._&n., dir=&githubpath., file=LXPMLE_&ncov._&nbeta._&n..csv);
%data2r(data=poislrci.ML_&ncov._&nbeta._&n., dir=&githubpath., file=ML_&ncov._&nbeta._&n..csv);

%mend;


