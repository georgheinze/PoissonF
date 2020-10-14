%macro MLinportion(portions=10, nsim=10000);


data poislrci.ML_FUT_&ncov._&nbeta._&n.;
run;


%let portionsize=%sysevalf(&nsim / &portions);

%do iportion=1 %to &portions;

	data simdata_portion;
	set simdata10;
	if isim > (&iportion -1) * &portionsize AND isim <= &iportion * &portionsize;
	run;

	ods select none;
	proc genmod data=simdata_portion;
	ods output parameterestimates=parm_portion;
	model y_fut=&xlist /dist=poisson offset=logFU;
	by isim;
	run;
	ods select all;


	data poislrci.ML_FUT_&ncov._&nbeta._&n.;
	set poislrci.ML_FUT_&ncov._&nbeta._&n. parm_portion;
	run;


%end;

%mend;


