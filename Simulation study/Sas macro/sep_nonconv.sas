**** extract separation rates and estimation failures from individual results;
**** GH 13/10/2020;

%macro sep_nonconv(all_ncov=2 5 10, all_epv=3 5 10, all_nbeta=1 2 3 4 5 6 7 8 9, methods=ml lxpmle logxact, se_thresh=10);

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


data _bigresults;
run;

%do j_all_ncov = 1 %to &n_all_ncov;
	%let ncov = &&ncov&j_all_ncov;
	%do j_all_epv = 1 %to &n_all_epv;
		%let n=%sysevalf(&&epv&j_all_epv * 10 * &ncov);
		%let do_logxact=1;
		%let do_lxpmle=1;
		%if &n > 250 %then %do; %let do_logxact=0; %end;   *** did not compute logXact for n>250;
		%if &n > 500 & &ncov > 5 %then %do; %let do_lxpmle=0; %end; * did not compute LXPMLE for N>500 and K>5;
		%do j_all_nbeta=1 %to &n_all_nbeta;
			%let nbeta= &&nbeta&j_all_nbeta;
			%let postfix=_&ncov._&nbeta._&n;
			
			%put "NOTE: &postfix";

			* ML;
			%put "NOTE: ...ML";
			data _tmp;
			set poislrci.ML&postfix;
			truebeta1 = (&nbeta - 5)*log(2);
			truebeta2 = log(2);
			if Parameter="x1" then do; 
				if stderr>&SE_thresh then separation1=1;
				else separation1=0;
			end;
			else if Parameter="x2" then do;
				if stderr>&SE_thresh then separation2=1;
				else separation2=0;
			end;
			run;

			proc means data=_tmp noprint;
			var separation1 separation2;
			output out=_tmp_mean mean=separation1 separation2;
			run;

		
			data _tmp_merge&postfix;
			set _tmp_mean;
			n=&n;
			ncov=&ncov;
			epv=&&epv&j_all_epv;
			truebeta1 = (&nbeta - 5)*log(2);
			truebeta2 = log(2);
			method="ML        ";
			run;

			data _bigresults;
			set _bigresults _tmp_merge&postfix;
			run;

			* LXPMLE;
			%put "NOTE: ...LXPMLE";

			%if &do_lxpmle = 1 %then %do;
				data _tmp;
				set poislrci.LXPMLE&postfix;
				truebeta1 = (&nbeta - 5)*log(2);
				truebeta2 = log(2);
				if Effect="x1" then do; 
					if type2="AS" then do; 
						if beta = . then BetaNA1=1;
						else BetaNA1=0;
					end;
					if type2="Profile" then do;
						if LowerCI = . then NoLowerCI1=1;
						else NoLowerCI1=0;
						if UpperCI = . then NoUpperCI1=1;
						else NoUpperCI1=0;
					end;
				end;
				else if Effect="x2" then do;
					if type2="AS" then do; 
						if beta = . then BetaNA2=1;
						else BetaNA2=0;
					end;
					if type2="Profile" then do;
						if LowerCI = . then NoLowerCI2=1;
						else NoLowerCI2=0;
						if UpperCI = . then NoUpperCI2=1;
						else NoUpperCI2=0;
					end;
				end;
				run;

				proc means data=_tmp noprint;
				var BetaNA1 BetaNA2  noLowerCI1 NoUpperCI1 NoLowerCI2 NoUpperCI2;
				output out=_tmp_mean mean=BetaNA1 BetaNA2  noLowerCI1 NoUpperCI1 NoLowerCI2 NoUpperCI2;
				run;


		
				data _tmp_merge&postfix;
				set _tmp_mean;
				n=&n;
				ncov=&ncov;
				epv=&&epv&j_all_epv;
				truebeta1 = (&nbeta - 5)*log(2);
				truebeta2 = log(2);
				method="LXPMLE    ";
				run;

				data _bigresults;
				set _bigresults _tmp_merge&postfix;
				run;
			%end;


			* LogXact;
			%if &do_LogXact = 1 %then %do;
				* Exact CI;
				%put "NOTE: ...Exact";

				data _tmp;
				set poislrci.LogXact&postfix;
				truebeta1 = (&nbeta - 5)*log(2);
				truebeta2 = log(2);
				if substr(Effect,1,2)="x1" then do; 
					if type2="Exact" then do; 
						if beta = . then BetaNA1=1;
						else BetaNA1=0;
						if Effect="x1                      CMLE" then MUE1=1;
						else MUE1=0;
					end;
				end;
				else if substr(Effect,1,2)="x2" then do;
					if type2="Exact" then do; 
						if beta = . then BetaNA2=1;
						else BetaNA2=0;
						if Effect="x1                      CMLE" then MUE2=1;
						else MUE2=0;
					end;
				end;
				run;

				proc means data=_tmp noprint;
				var BetaNA1 MUE1 BetaNA2 MUE2;
				output out=_tmp_mean mean=BetaNA1 MUE1 BetaNA2 MUE2;
				run;



				data _tmp_merge&postfix;
				set _tmp_mean;
				n=&n;
				ncov=&ncov;
				epv=&&epv&j_all_epv;
				truebeta1 = (&nbeta - 5)*log(2);
				truebeta2 = log(2);
				method="Exact    ";
				run;

				data _bigresults;
				set _bigresults _tmp_merge&postfix;
				run;


			%end;


		%end;
	%end;
%end;


%mend;
