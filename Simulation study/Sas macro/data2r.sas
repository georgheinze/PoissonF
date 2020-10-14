%macro data2R(data=_last_,  varlist=, file=, dir=%str(E:/work/), if=%str(1=1), cco=0);

%if &file= %then %let file=&data;
%if %upcase(&file)=_LAST_ %then %let file=data2r.txt;

%if &varlist ne %then %do;

 %let nvar=0;
 %do %while(%scan(&varlist,&nvar+1)~=);
  %let nvar=%eval(&nvar+1);
  %let var&nvar=%scan(&varlist,&nvar);
 %end;


 data _work;
 set &data;
 %if &cco=1 %then %do;
  %do j=1 %to &nvar;
   if &&var&j ne .;
  %end;
%end;
 if &if;
 keep &varlist;
 run;

%end;
%else %do;
 data _work;
 set &data;
  if &if;
 run;
%end;

PROC EXPORT DATA= WORK._work 
            OUTFILE= "&dir.&file" 
            DBMS=TAB 
			dbms=dlm replace;
      delimiter=';' ;
RUN;




%put NOTE: Open the file in R using:;
%put;
%put &file <- read.table("&dir.&file", header=T, sep=";") ;
%put head(&file) ;
%put;


%mend;

