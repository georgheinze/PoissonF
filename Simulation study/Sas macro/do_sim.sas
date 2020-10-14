%macro do_sim(all_ncov=2, all_epv=3 5 10, all_nbeta=1 2 3 4 5 6 7 8 9);

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
	%let xlist =;
	%do j_x = 1 %to &ncov;
		%let xlist= &xlist x&j_x;
	%end;
	%do j_all_epv = 1 %to &n_all_epv;
		%let n=%sysevalf(&&epv&j_all_epv * 10 * &ncov);
		%do j_all_nbeta=1 %to &n_all_nbeta;
			%let nbeta= &&nbeta&j_all_nbeta;
			%include "&githubpath.PoissonF\Simulation study\Follow-up time batch\run_simu_fut.sas";
		%end;
	%end;
%end;
%mend;
