**** extract simulation results from individual results;
**** GH 09/10/2020;

%macro extr_results(all_ncov=2 5 10, all_epv=3 5 10, all_nbeta=1 2 3 4 5 6 7 8 9, methods=ml dataug firth flac lxpmle logxact);

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
		%if &n > 250 %then %do; %let do_logxact=0; %end;   *** do not compute logXact for n>250;
		%if &n > 500 & &ncov > 5 %then %do; %let do_lxpmle=0; %end; * do not compute LXPMLE for N>500 and K>5;
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
				bias1 = Estimate - truebeta1;
				MSEN1 = (Estimate - truebeta1)**2;
				leftCov1 = (LowerWaldCL < truebeta1) & (LowerWaldCL ne .);
				rightCov1 = (UpperWaldCL > truebeta1) & (UpperWaldCL ne .);
				power1 = (lowerWaldCL ne .) & (UpperWaldCL ne .) & ((lowerWaldCL > 0) | (upperwaldCL < 0));
				width1 = UpperWaldCL-LowerWaldCL;
			end;
			else if Parameter="x2" then do;
				bias2 = Estimate - truebeta2;
				MSEN2 = (Estimate - truebeta2)**2;
				leftCov2 = (LowerWaldCL < truebeta2) & (LowerWaldCL ne .);
				rightCov2 = (UpperWaldCL > truebeta2) & (UpperWaldCL ne .);
				power2 = (lowerWaldCL ne .) & (UpperWaldCL ne .) & ((lowerWaldCL > 0) | (upperwaldCL < 0));
				width2 = UpperWaldCL-LowerWaldCL;
			end;
			run;

			proc means data=_tmp noprint;
			var bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
			output out=_tmp_mean mean=bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
			run;

			proc means data=_tmp noprint;
			var width1 width2;
			output out=_tmp_median median=width1 width2;
			run;

			data _tmp2;
			set poislrci.ML_pred&postfix;
			MSPEN = (_MLPred - mu_FUT)**2;
			run;

			proc means data=_tmp2 noprint;
			var MSPEN;
			output out=_tmp2_agg sum=MSPEN;
			by isim;
			where isim ne .;
			run;

			proc means data=_tmp2_agg noprint;
			var MSPEN;
			output out=_tmp2_agg2 mean=MSPEN;
			run;

			data _tmp_merge&postfix;
			merge _tmp_mean _tmp_median _tmp2_agg2;
			rmsen1=sqrt(msen1 * &n);
			rmsen2=sqrt(msen2 * &n);
			rmspen = sqrt(mspen);
			n=&n;
			ncov=&ncov;
			epv=&&epv&j_all_epv;
			method="ML        ";
			truebeta1 = (&nbeta - 5)*log(2);
			truebeta2 = log(2);
			run;

			data _bigresults;
			set _bigresults _tmp_merge&postfix;
			run;



			* Firth;
			%put "NOTE: ...Firth";

			data _tmp;
			set poislrci.Firth_fut&postfix;
			truebeta1 = (&nbeta - 5)*log(2);
			truebeta2 = log(2);
			if Parameter="x1" then do; 
				bias1 = Estimate - truebeta1;
				MSEN1 = (Estimate - truebeta1)**2;
				leftCov1 = (LowerLRCL < truebeta1) & (LowerLRCL ne .);
				rightCov1 = (UpperLRCL > truebeta1) & (UpperLRCL ne .);
				power1 = (lowerLRCL ne .) & (UpperLRCL ne .) & ((lowerLRCL > 0) | (upperLRCL < 0));
				width1 = UpperLRCL-LowerLRCL;
			end;
			else if Parameter="x2" then do;
				bias2 = Estimate - truebeta2;
				MSEN2 = (Estimate - truebeta2)**2;
				leftCov2 = (LowerLRCL < truebeta2) & (LowerLRCL ne .);
				rightCov2 = (UpperLRCL > truebeta2) & (UpperLRCL ne .);
				power2 = (lowerLRCL ne .) & (UpperLRCL ne .) & ((lowerLRCL > 0) | (upperLRCL < 0));
				width2 = UpperLRCL-LowerLRCL;
			end;
			run;

			proc means data=_tmp noprint;
			var bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
			output out=_tmp_mean mean=bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
			run;

			proc means data=_tmp noprint;
			var width1 width2;
			output out=_tmp_median median=width1 width2;
			run;

			data _tmp2;
			set poislrci.Firth_fut_pred&postfix;
			MSPEN = (_FIRTHPred - mu_FUT)**2;
			run;

			proc means data=_tmp2 noprint;
			var MSPEN;
			output out=_tmp2_agg sum=MSPEN;
			by isim;
			where isim ne .;
			run;

			proc means data=_tmp2_agg noprint;
			var MSPEN;
			output out=_tmp2_agg2 mean=MSPEN;
			run;

			data _tmp_merge&postfix;
			merge _tmp_mean _tmp_median _tmp2_agg2;
			rmsen1=sqrt(msen1 * &n);
			rmsen2=sqrt(msen2 * &n);
			rmspen = sqrt(mspen);
			n=&n;
			ncov=&ncov;
			epv=&&epv&j_all_epv;
			truebeta1 = (&nbeta - 5)*log(2);
			truebeta2 = log(2);
			method="Firth     ";
			run;

			data _bigresults;
			set _bigresults _tmp_merge&postfix;
			run;




			* FLAC;
			%put "NOTE: ...FLAC";

			data _tmp;
			set poislrci.FLAC_fut&postfix;
			truebeta1 = (&nbeta - 5)*log(2);
			truebeta2 = log(2);
			if Parameter="x1" then do; 
				bias1 = Estimate - truebeta1;
				MSEN1 = (Estimate - truebeta1)**2;
				leftCov1 = (LowerLRCL < truebeta1) & (LowerLRCL ne .);
				rightCov1 = (UpperLRCL > truebeta1) & (UpperLRCL ne .);
				power1 = (lowerLRCL ne .) & (UpperLRCL ne .) & ((lowerLRCL > 0) | (upperLRCL < 0));
				width1 = UpperLRCL-LowerLRCL;
			end;
			else if Parameter="x2" then do;
				bias2 = Estimate - truebeta2;
				MSEN2 = (Estimate - truebeta2)**2;
				leftCov2 = (LowerLRCL < truebeta2) & (LowerLRCL ne .);
				rightCov2 = (UpperLRCL > truebeta2) & (UpperLRCL ne .);
				power2 = (lowerLRCL ne .) & (UpperLRCL ne .) & ((lowerLRCL > 0) | (upperLRCL < 0));
				width2 = UpperLRCL-LowerLRCL;
			end;
			run;

			proc means data=_tmp noprint;
			var bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
			output out=_tmp_mean mean=bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
			run;

			proc means data=_tmp noprint;
			var width1 width2;
			output out=_tmp_median median=width1 width2;
			run;

			data _tmp2;
			set poislrci.FLAC_fut_pred&postfix;
			MSPEN = (_FLACPred - mu_FUT)**2;
			run;

			proc means data=_tmp2 noprint;
			var MSPEN;
			output out=_tmp2_agg sum=MSPEN;
			by isim;
			where isim ne .;
			run;

			proc means data=_tmp2_agg noprint;
			var MSPEN;
			output out=_tmp2_agg2 mean=MSPEN;
			run;

			data _tmp_merge&postfix;
			merge _tmp_mean _tmp_median _tmp2_agg2;
			rmsen1=sqrt(msen1 * &n);
			rmsen2=sqrt(msen2 * &n);
			rmspen = sqrt(mspen);
			n=&n;
			ncov=&ncov;
			epv=&&epv&j_all_epv;
			truebeta1 = (&nbeta - 5)*log(2);
			truebeta2 = log(2);
			method="FLAC      ";
			run;

			data _bigresults;
			set _bigresults _tmp_merge&postfix;
			run;



			* DatAug;
			%put "NOTE: ...DatAug";

			data _tmp;
			set poislrci.Dataug_fut&postfix;
			truebeta1 = (&nbeta - 5)*log(2);
			truebeta2 = log(2);
			if Parameter="x1" then do; 
				bias1 = Estimate - truebeta1;
				MSEN1 = (Estimate - truebeta1)**2;
				leftCov1 = (LowerLRCL < truebeta1) & (LowerLRCL ne .);
				rightCov1 = (UpperLRCL > truebeta1) & (UpperLRCL ne .);
				power1 = (lowerLRCL ne .) & (UpperLRCL ne .) & ((lowerLRCL > 0) | (upperLRCL < 0));
				width1 = UpperLRCL-LowerLRCL;
			end;
			else if Parameter="x2" then do;
				bias2 = Estimate - truebeta2;
				MSEN2 = (Estimate - truebeta2)**2;
				leftCov2 = (LowerLRCL < truebeta2) & (LowerLRCL ne .);
				rightCov2 = (UpperLRCL > truebeta2) & (UpperLRCL ne .);
				power2 = (lowerLRCL ne .) & (UpperLRCL ne .) & ((lowerLRCL > 0) | (upperLRCL < 0));
				width2 = UpperLRCL-LowerLRCL;
			end;
			run;

			proc means data=_tmp noprint;
			var bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
			output out=_tmp_mean mean=bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
			run;

			proc means data=_tmp noprint;
			var width1 width2;
			output out=_tmp_median median=width1 width2;
			run;

			data _tmp2;
			set poislrci.Dataug_fut_pred&postfix;
			MSPEN = (_DATAUGPred - mu_FUT)**2;
			run;

			proc means data=_tmp2 noprint;
			var MSPEN;
			output out=_tmp2_agg sum=MSPEN;
			by isim;
			where isim ne .;
			run;

			proc means data=_tmp2_agg noprint;
			var MSPEN;
			output out=_tmp2_agg2 mean=MSPEN;
			run;

			data _tmp_merge&postfix;
			merge _tmp_mean _tmp_median _tmp2_agg2;
			rmsen1=sqrt(msen1 * &n);
			rmsen2=sqrt(msen2 * &n);
			rmspen = sqrt(mspen);
			n=&n;
			ncov=&ncov;
			epv=&&epv&j_all_epv;
			truebeta1 = (&nbeta - 5)*log(2);
			truebeta2 = log(2);
			method="DatAug    ";
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
						bias1 = beta - truebeta1;
						MSEN1 = (beta - truebeta1)**2;
					end;
					if type2="Profile" then do;
						leftCov1 = (LowerCI < truebeta1) & (LowerCI ne .);
						rightCov1 = (UpperCI > truebeta1) & (UpperCI ne .);
						power1 = (lowerCI ne .) & (UpperCI ne .) & ((lowerCI > 0) | (upperCI < 0));
						width1 = UpperCI-LowerCI;
					end;
				end;
				else if Effect="x2" then do;
					if type2="AS" then do; 
						bias2 = beta - truebeta2;
						MSEN2 = (beta - truebeta2)**2;
					end;
					if type2="Profile" then do;
						leftCov2 = (LowerCI < truebeta2) & (LowerCI ne .);
						rightCov2 = (UpperCI > truebeta2) & (UpperCI ne .);
						power2 = (lowerCI ne .) & (UpperCI ne .) & ((lowerCI > 0) | (upperCI < 0));
						width2 = UpperCI-LowerCI;
					end;
				end;
				run;

				proc means data=_tmp noprint;
				var bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
				output out=_tmp_mean mean=bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
				run;

				proc means data=_tmp noprint;
				var width1 width2;
				output out=_tmp_median median=width1 width2;
				run;

		
				data _tmp_merge&postfix;
				merge _tmp_mean _tmp_median;
				rmsen1=sqrt(msen1 * &n);
				rmsen2=sqrt(msen2 * &n);
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
						bias1 = beta - truebeta1;
						MSEN1 = (beta - truebeta1)**2;
					end;
					if type2="Exact" then do;
						leftCov1 = (LowerCI < truebeta1) & (LowerCI ne .);
						rightCov1 = (UpperCI > truebeta1) & (UpperCI ne .);
						power1 = (lowerCI ne .) & (UpperCI ne .) & ((lowerCI > 0) | (upperCI < 0));
						width1 = UpperCI-LowerCI;
					end;
				end;
				else if substr(Effect,1,2)="x2" then do;
					if type2="Exact" then do; 
						bias2 = beta - truebeta2;
						MSEN2 = (beta - truebeta2)**2;
					end;
					if type2="Exact" then do;
						leftCov2 = (LowerCI < truebeta2) & (LowerCI ne .);
						rightCov2 = (UpperCI > truebeta2) & (UpperCI ne .);
						power2 = (lowerCI ne .) & (UpperCI ne .) & ((lowerCI > 0) | (upperCI < 0));
						width2 = UpperCI-LowerCI;
					end;
				end;
				run;

				proc means data=_tmp noprint;
				var bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
				output out=_tmp_mean mean=bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
				run;

				proc means data=_tmp noprint;
				var width1 width2;
				output out=_tmp_median median=width1 width2;
				run;


				data _tmp_merge&postfix;
				merge _tmp_mean _tmp_median;
				rmsen1=sqrt(msen1 * &n);
				rmsen2=sqrt(msen2 * &n);
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


				* Mid-P CI;
				%put "NOTE: ...Mid-P";
	
				data _tmp;
				set poislrci.LogXact&postfix;
				truebeta1 = (&nbeta - 5)*log(2);
				truebeta2 = log(2);
				if substr(Effect,1,2)="x1" then do; 
					if type2="Exact" then do; 
						bias1 = beta - truebeta1;
						MSEN1 = (beta - truebeta1)**2;
					end;
					if type2="Ex-MidP" then do;
						leftCov1 = (LowerCI < truebeta1) & (LowerCI ne .);
						rightCov1 = (UpperCI > truebeta1) & (UpperCI ne .);
						power1 = (lowerCI ne .) & (UpperCI ne .) & ((lowerCI > 0) | (upperCI < 0));
						width1 = UpperCI-LowerCI;
					end;
				end;
				else if substr(Effect,1,2)="x2" then do;
					if type2="Exact" then do; 
						bias2 = beta - truebeta2;
						MSEN2 = (beta - truebeta2)**2;
					end;
					if type2="Ex-MidP" then do;
						leftCov2 = (LowerCI < truebeta2) & (LowerCI ne .);
						rightCov2 = (UpperCI > truebeta2) & (UpperCI ne .);
						power2 = (lowerCI ne .) & (UpperCI ne .) & ((lowerCI > 0) | (upperCI < 0));
						width2 = UpperCI-LowerCI;
					end;
				end;
				run;

				proc means data=_tmp noprint;
				var bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
				output out=_tmp_mean mean=bias1 msen1 leftcov1 rightcov1 power1 bias2 msen2 leftcov2 rightcov2 power2;
				run;

				proc means data=_tmp noprint;
				var width1 width2;
				output out=_tmp_median median=width1 width2;
				run;


				data _tmp_merge&postfix;
				merge _tmp_mean _tmp_median;
				rmsen1=sqrt(msen1 * &n);
				rmsen2=sqrt(msen2 * &n);
				n=&n;
				ncov=&ncov;
				epv=&&epv&j_all_epv;
				truebeta1 = (&nbeta - 5)*log(2);
				truebeta2 = log(2);
				method="Mid-P    ";
				run;

				data _bigresults;
				set _bigresults _tmp_merge&postfix;
				run;
			%end;


		%end;
	%end;
%end;


%mend;
