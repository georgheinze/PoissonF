%macro do_sim_ex(all_ncov=2, all_epv=3 5 10, all_nbeta=1 2 3 4 5 6 7 8 9, exact=1, lxpmle=1, ml=1);

%let n_all_ncov=0;
%do %while(%scan(&all_ncov,&n_all_ncov+1)~=);
 %let n_all_ncov=%eval(&n_all_ncov+1);
 %let ncov&n_all_ncov=%scan(&all_ncov,&n_all_ncov);
%end;


%let n_all_epv=0;
%do %while(%scan(&all_epv,&n_all_epv+1)~=);
 %let n_all_epv=%eval(&n_all_epv+1);
 %let epv&n_all_epv=%scan(&all_epv,&n_all_epv);
%end;



%let n_all_nbeta=0;
%do %while(%scan(&all_nbeta,&n_all_nbeta+1)~=);
 %let n_all_nbeta=%eval(&n_all_nbeta+1);
 %let nbeta&n_all_nbeta=%scan(&all_nbeta,&n_all_nbeta);
%end;



%do j_all_ncov = 1 %to &n_all_ncov;
	%let ncov = &&ncov&j_all_ncov;
	%if &ncov<=5 %then %do; %let do_exact=&exact;   * will be used as input in %mlexinportion(); %end;
	%else %do; %let do_exact=0; %end;
	%let xlist =;
	%do j_x = 1 %to &ncov;
		%let xlist= &xlist x&j_x;
	%end;
	%do j_all_epv = 1 %to &n_all_epv;
		%let n=%sysevalf(&&epv&j_all_epv * 10 * &ncov);
		%if &n > 250 %then %do; %let do_exact=0; %end;   *** do not compute logXact for n>250;
		%let do_lxpmle=&lxpmle;
		%if &n > 500 & &ncov > 5 %then %do; %let do_lxpmle=0; %end; * do not compute LXPMLE for N>500 and K>5;
		%do j_all_nbeta=1 %to &n_all_nbeta;
			%let nbeta= &&nbeta&j_all_nbeta;
			%include "&githubpath.PoissonF\Simulation study\Follow-up time batch\run_simu_ex.sas";
		%end;
	%end;
%end;
%mend;
