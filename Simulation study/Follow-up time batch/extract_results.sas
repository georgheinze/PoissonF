%let Githubpath = %str(C:\Users\Georg\Documents\Github\);
%let Localpath = %str(C:\Users\Georg\Documents\Local\);

libname poislrci "&Localpath.PoissonF\Simulation study\Follow-up time batch\Results";

%include "&githubpath.PoissonF\Simulation study\Sas macro\data2r.sas";
%include "&githubpath.PoissonF\Simulation study\Sas macro\extr_results.sas";
%include "&githubpath.PoissonF\Simulation study\Sas macro\sep_nonconv.sas";


options nomprint nomlogic nomacrogen nonotes nosource;

%extr_results();

options notes source;

data poislrci.bigresults;
set _bigresults;
run;


%data2r(data=poislrci.bigresults, dir=%str(&githubpath.PoissonF\Simulation study\Follow-up time batch\Summaries\), file=PoissonFResultsTable.csv);

options nomprint nomlogic nomacrogen nonotes nosource;

%sep_nonconv();

options notes source;

data poislrci.sep_nonconv;
set _bigresults;
run;

%data2r(data=poislrci.sep_nonconv, dir=%str(&githubpath.PoissonF\Simulation study\Follow-up time batch\Summaries\), file=PoissonF_sep_nonconv.csv);
