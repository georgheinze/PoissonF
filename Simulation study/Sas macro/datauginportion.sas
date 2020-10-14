%macro datauginportion(portions=10, nsim=10000, y_2go=y_fut, FUT=1, keep=);

%if &FUT=1 %then %let INSERT=_FUT;
%else %let INSERT=;

data poislrci.DATAUG&INSERT._&ncov._&nbeta._&n.;
run;

data poislrci.DATAUG&INSERT._PRED_&ncov._&nbeta._&n.;
run;


%let portionsize=%sysevalf(&nsim / &portions);

%do iportion=1 %to &portions;

	data simdata_portion;
	set simdata10;
	if isim > (&iportion -1) * &portionsize AND isim <= &iportion * &portionsize;
	run;
	options nonotes;
	%dataugpoisson(data=simdata_portion, y=&y_2go, varlist=&xlist , by=isim, offset=logFUT, print=0, odsselect=all, id=&keep,
				priorinterval=1000 1000 1000 1000 1000 1000 100 100 100 100);  
	** for simus without follow-upt time, logFUT=0;
	options notes;

	data poislrci.DATAUG&INSERT._&ncov._&nbeta._&n.;
	set poislrci.DATAUG&INSERT._&ncov._&nbeta._&n. _DATAUGparms;
	run;

	data poislrci.DATAUG&INSERT._PRED_&ncov._&nbeta._&n.;
	set poislrci.DATAUG&INSERT._PRED_&ncov._&nbeta._&n. _DATAUGPredictions;
	run;



%end;





*** export to csv;

%data2r(data=poislrci.FLAC&INSERT._&ncov._&nbeta._&n., dir=%str(&githubpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=FLAC&INSERT._&ncov._&nbeta._&n..csv);
%data2r(data=poislrci.FIRTH&INSERT._&ncov._&nbeta._&n., dir=%str(&githubpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=FIRTH&INSERT._&ncov._&nbeta._&n..csv);

%data2r(data=poislrci.FLAC&INSERT._PRED_&ncov._&nbeta._&n., dir=%str(&githubpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=FLAC&INSERT._PRED_&ncov._&nbeta._&n..csv);
%data2r(data=poislrci.FIRTH&INSERT._PRED_&ncov._&nbeta._&n., dir=%str(&githubpath.PoissonF\Simulation study\Follow-up time batch\Results\), file=FIRTH&INSERT._PRED_&ncov._&nbeta._&n..csv);



%mend;
