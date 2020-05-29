%let Githubpath = %str(E:\1Projects\2018-08 NegBinom FC\Github\);

%inc "&githubpath.PoissonF\Simulation study\Sas macro\array.sas";
%inc "&githubpath.PoissonF\Simulation study\Sas macro\do_over.sas";

%array(ncov, values=2 2 2 2 2 2 2 2 2 2);
%array(epv, values=3 3 3 3 3 10 10 10 10 10);
%array(nbeta, values=1 2 3 4 5 1 2 3 4 5);
%array(what, values=Firth Firth Firth Firth Firth Firth Firth Firth Firth Firth);



%macro exp_parms(ncov, nbeta, epv, what);
	%let n=%sysevalf(&epv*&ncov*10);
				PROC EXPORT DATA= POISLRCI.&what._fut_&ncov._&nbeta._&n. 
				            OUTFILE= "&githubpath.PoissonF\Simulation study\Follow-up time batch\&what._FUT_&ncov._&nbeta._&n..csv" 
				            DBMS=CSV REPLACE;
				     PUTNAMES=YES;
				RUN;
%mend;

options mprint mlogic macrogen;
%do_over(ncov nbeta epv what, macro=exp_parms); 

%macro summary(ncov, nbeta, epv, what);
	%if %upcase(&what)=FIRTH %then %let CI=LRCL;
	%else %if %upcase(&what)=FLAC %then %let CI=LRCL;
	%else %if %upcase(&what)=ML %then %let CI=WaldCL;
	%let n=%sysevalf(&epv*&ncov*10);

	data summary_i;
	set poislrci.&what._FUT_&ncov._&nbeta._&n.;
	array beta(5);
	beta(1)=-2.77;
	beta(2)=-0.69;
	beta(3)=0;
	beta(4)=0.69;
	beta(5)=2.77;
	if Parameter="x1";                      ****** <<<<< change for different beta >>>>>> *******;
	bias=estimate-beta(&nbeta.);
	MSE=(estimate-beta(&nbeta.))**2;
	lowercover=(lower&CI. < beta&nbeta.);
	uppercover=(upper&CI. > beta&nbeta.);
	power=(lower&CI. > 0 ) | (upper&CI. < 0);
	run;
	proc means data=summary_i noprint;
	var bias MSE lowercover uppercover power;
	output out=poislrci.summary_&what._FUT_&ncov._&nbeta._&n. mean=bias MSE lowercover uppercover power;
	run;
	data poislrci.summary_&what._FUT_&ncov._&nbeta._&n.;
	set poislrci.summary_&what._FUT_&ncov._&nbeta._&n.;
	ncov=&ncov;
	nbeta=&nbeta;
	epv=&epv;
	n=&n;
	run;
%mend;

%do_over(ncov nbeta epv what, macro=summary);
 

%macro exp_summary(what);

data summary_FUT_2_x_n;
set   poislrci.summary_&what._FUT_2_1_60 
poislrci.summary_&what._FUT_2_2_60 
poislrci.summary_&what._FUT_2_3_60 
poislrci.summary_&what._FUT_2_4_60 
poislrci.summary_&what._FUT_2_5_60 
poislrci.summary_&what._FUT_2_1_200
poislrci.summary_&what._FUT_2_2_200
poislrci.summary_&what._FUT_2_3_200
poislrci.summary_&what._FUT_2_4_200
poislrci.summary_&what._FUT_2_5_200
;
run;




data poislrci.summary_&what._FUT;
set summary_FUT_2_x_n;
RMSE=sqrt(MSE);
RMSE_SQRTN=RMSE*sqrt(n);
run;



PROC EXPORT DATA= POISLRCI.summary_&what._FUT 
				            OUTFILE= "&githubpath.PoissonF\Simulation study\Follow-up time batch\summary_&what._FUT.csv" 
				            DBMS=CSV REPLACE;
				     PUTNAMES=YES;
				RUN;


%mend;

%exp_summary(what=Firth);


%array(what, values=Ml Ml Ml Ml Ml Ml Ml Ml Ml Ml);
%do_over(ncov nbeta epv what, macro=exp_parms); 
%do_over(ncov nbeta epv what, macro=summary);


%exp_summary(what=Ml);
