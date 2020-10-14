%macro flacinportion(portions=10, nsim=10000, y_2go=y_fut, FUT=1, keep=);

%if &FUT=1 %then %let INSERT=_FUT;
%else %let INSERT=;

data poislrci.FLAC&INSERT._&ncov._&nbeta._&n.;
run;

data poislrci.FIRTH&INSERT._&ncov._&nbeta._&n.;
run;
data poislrci.FLAC&INSERT._PRED_&ncov._&nbeta._&n.;
run;

data poislrci.FIRTH&INSERT._PRED_&ncov._&nbeta._&n.;
run;

%let portionsize=%sysevalf(&nsim / &portions);

%do iportion=1 %to &portions;

	data simdata_portion;
	set simdata10;
	if isim > (&iportion -1) * &portionsize AND isim <= &iportion * &portionsize;
	run;
	options nonotes;
	%flacpoisson(data=simdata_portion, y=&y_2go, varlist=&xlist , by=isim, maxiter=100, offset=logFUT, print=0, odsselect=all, nonotes=1, keep=&keep);  
	** for simus without follow-upt time, logFUT=0;
	options notes;

	data poislrci.FLAC&INSERT._&ncov._&nbeta._&n.;
	set poislrci.FLAC&INSERT._&ncov._&nbeta._&n. _FLACparms;
	run;

	data poislrci.FLAC&INSERT._PRED_&ncov._&nbeta._&n.;
	set poislrci.FLAC&INSERT._PRED_&ncov._&nbeta._&n. _FLACPredictions;
	run;

	data poislrci.FIRTH&INSERT._&ncov._&nbeta._&n.;
	set poislrci.FIRTH&INSERT._&ncov._&nbeta._&n. _FIRTHparms;
	run;

	data poislrci.FIRTH&INSERT._PRED_&ncov._&nbeta._&n.;
	set poislrci.FIRTH&INSERT._PRED_&ncov._&nbeta._&n. _FIRTHPredictions;
	run;


%end;





*** export to csv;

%data2r(data=poislrci.FLAC&INSERT._&ncov._&nbeta._&n., dir=%str(&githubpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=FLAC&INSERT._&ncov._&nbeta._&n..csv);
%data2r(data=poislrci.FIRTH&INSERT._&ncov._&nbeta._&n., dir=%str(&githubpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=FIRTH&INSERT._&ncov._&nbeta._&n..csv);

%data2r(data=poislrci.FLAC&INSERT._PRED_&ncov._&nbeta._&n., dir=%str(&githubpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=FLAC&INSERT._PRED_&ncov._&nbeta._&n..csv);
%data2r(data=poislrci.FIRTH&INSERT._PRED_&ncov._&nbeta._&n., dir=%str(&githubpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=FIRTH&INSERT._PRED_&ncov._&nbeta._&n..csv);



%mend;
