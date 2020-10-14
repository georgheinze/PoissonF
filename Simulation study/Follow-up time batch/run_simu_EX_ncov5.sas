%let Githubpath = %str(C:\Users\Georg\Documents\Github\);
%let Localpath = %str(C:\Users\Georg\Documents\Local\);

%include "&githubpath.PoissonF\Simulation study\Sas macro\do_sim_ex.sas";


libname poislrci "&Localpath.PoissonF\Simulation study\Follow-up time batch\Results";

%do_sim_ex(all_ncov=5, all_epv=3 5 10, all_nbeta=1 2 3 4 5 6 7 8 9);
*%do_sim_ex(all_ncov=5, all_epv=3 5, all_nbeta=1 2 3 4 5 6 7 8 9, exact=0, ml=0);
*%do_sim_ex(all_ncov=5, all_epv=10, all_nbeta=1 2 3 4 5 6 7 8 9, exact=0, ml=1);

