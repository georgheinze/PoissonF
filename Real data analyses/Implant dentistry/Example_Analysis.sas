*** Analysis of Implant Dentistry study;
*** The study is described in detail, including results from a GEE analysis, in Feher et al, Clin Oral Implants Res 2020, doi:10.1111/clr.13636;
*** Here we use only some of the variables on an aggregated level;

*** Outcome variable: 
	Hema ... number of hematological complications;
*** Independent variables: 
	Light_vs_no ... Dummy variable for light or heavy smoking (1/0)
	Heavy_vs_light ... dummy variable for heavy smoking (1/=0)
 	* if used together, these two variables encode the effects of light vs. no smoking, and of heavy vs. light smoking in ordinal coding *
    Diabetes ... Patient has diabetes (1/0)
	Age_decade ... age in full decades, centered to 50 years (-4 to 3);
*** Rate multiplier:
	Implants ... number of implantations
	log_implants ... logarithm of number of implantations (to be used as offset);

*** Each line of data may refer to up to 199 patients;
*** This information has been removed to protect privacy. For Poisson regression, data on the aggregated level carries the
	same information as data on the individual-patient level.;




	

*** read SAS macros;

filename dataug url 'https://raw.githubusercontent.com/georgheinze/flicflac/master/PoissonRegression/dataugpoisson.sas' ;
%include dataug;

filename flacpois url 'https://raw.githubusercontent.com/georgheinze/flicflac/master/PoissonRegression/flacpoisson.sas' ;
%include flacpois;


*** read data set;

filename impdent temp;
proc http
url="https://raw.githubusercontent.com/georgheinze/PoissonF/master/Real data analyses/Implant dentistry/Impdent.csv" method=get out=impdent;
run;

proc import file=impdent out=work.Impdent replace dbms=csv;
run;


* Maximum likelihood analysis and exact Poisson regression;

proc genmod data=Impdent;
model Hema = Light_vs_no Heavy_vs_light Diabetes Age_decade / dist=poisson offset=log_Implants;
exact Light_vs_no / cltype=exact estimate;
exact Heavy_vs_light /cltype=exact estimate;
exact Diabetes /cltype=exact estimate;
exact Age_decade /cltype=exact estimate;
exact Light_vs_no / cltype=midp estimate;
exact Heavy_vs_light /cltype=midp estimate;
exact Diabetes /cltype=midp estimate;
exact Age_decade /cltype=midp estimate;
run;

* Analysis with FLAC ;

%flacpoisson(data=Impdent, y=Hema, offset=log_Implants, varlist= Light_vs_no Heavy_vs_light Diabetes Age_decade);


* Analysis with Bayesian data augmentation, setting prior intervals to [1/1000, 1000] for dichotomous variables and [1/100, 1/100] for age in decades;

%dataugpoisson(data=Impdent, y=Hema, offset=log_Implants, varlist= Light_vs_no Heavy_vs_light Diabetes Age_decade, priorinterval=1000 1000 1000 100);

* Analysis with Bayesian data augmentation, setting prior intervals to [1/50, 50] for dichotomous variables and [1/5, 1/5] for age in decades;

%dataugpoisson(data=Impdent, y=Hema, offset=log_Implants, varlist= Light_vs_no Heavy_vs_light Diabetes Age_decade, priorinterval=50 50 50 5);

* Analysis with exact Poisson regression;


proc genmod data=Impdent;
model Hema = Light_vs_no Heavy_vs_light Diabetes Age_decade / dist=poisson offset=log_Implants;
run;

