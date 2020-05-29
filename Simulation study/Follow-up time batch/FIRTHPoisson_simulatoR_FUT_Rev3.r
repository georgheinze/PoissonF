
### This R code perform simulations for 10 combinations of number of covariates, values of beta1, sample sizes using SimulatoR package
### Part 1 is simulation setup with SimulatoR
### Part 2 performs simulations for 2 covariates, 5 values of beta1, 2 different sample sizes(for 2 values of Expected number of events)
### Simulations are performed for data with follow-up time.
### Simulation data will be saved in SimDataDir. Please specify the path where you want to save the simulation data.

library("SimulatoR")
SimDataDir =  "E:/1Projects/2018-08 NegBinom FC/Github/PoissonF/Simulation study/Follow-up time batch/"

################################ Part 1 #######################################
transforms_pois<-list(v1 = function(z) (as.integer(z[,1] > 1.28)),    # change to Puhr, to achieve balance of 10:90
                      v2 = function(z) (as.integer(z[,2] > 0.35)))

relations_pois = cor_from_upper(2, rbind(c(1,2,0.5)))
sim_design_pois = design(relations=relations_pois, transform=transforms_pois, names_final=paste("x", 1:2, sep=""))
# truncate=c(NA, NA, NA),


################################ Function to specify condition on simulated data ##########
f1<-function(x){return(length(which(x!=0)))}
f2<-function(x){return(length(unique(x)))}

MoreThan3<-function (x) 
{
  x = as.matrix(x)
  val1 = any(apply(x, 2, f1) < 3)
  val2 = any(apply(x, 2, f2) < 1)
  if (val1) 
    warning("At least one variable had less than 3 non-zero values.\n")
  if (val2) 
    warning("At least one variable has no variation.\n")
  #print(val)
  val = val1 + val2
  return(val)
}


################################################################################################


simdata<-simulate_data(sim_design_pois, n_obs=10000)
# simdata<-conditional_simulate_data(relations_pois, n_obs=10000, reject = list(MoreThan3,is_collinear),
#                                    reject_max_iter = 100, on_reject = "stop",transform=transforms_pois, 
#                                    
#                                    names_final=paste("x", 1:10, sep="") )
# plot_initial_cor_network(sim_design_pois)   # plots network for the normal deviates
# plot_final_cor_network(sim_design_pois)  # plots network for the simulated 'real' covariates


################################ Part 2 #######################################
# start simulations
nsim =10000 # specify number of simulations
RNGkind(kind ="Mersenne-Twister", normal.kind ="Inversion", sample.kind ="Rejection" )
set.seed(79123)
RejSim =0
Interceptfile = paste(SimDataDir,"Intercept.csv",sep="")
write.table(t(c("ncov","Beta1","Intercept")),Interceptfile,sep=",",row.names=FALSE, col.names = FALSE)

for (ncov in c(2)){
  
  B1cnt =0 
  for(Beta1 in c(-log(16),-log(2),0, log(2),log(16))){
    B1cnt = B1cnt +1
    betaV<-numeric(2)
    a<-1   # in Puhr2017, was set to 0, 0.5, 1
    
    betaV[1] = Beta1
    betaV[2] = log(2)
    
    #Target function to determine intercept parameter:
    tar_fun<-function(Intercept, mu=0.1) {
      simdata<-simulate_data(sim_design_pois, n_obs=1000000)
      simdata<- simdata[,1:ncov]
      FUT <- round(runif(1000000, 1E-2,2),2)
      simdata$FUT<- FUT
      linpred<- as.matrix(simdata[,1:ncov]) %*% betaV + Intercept+log(FUT)    
      simdata$y <- rpois(length(linpred), exp(linpred))
      while( any(is.na(simdata$y)==TRUE)){
        # simdata<-conditional_simulate_data(relations_pois, n_obs=1000000, reject = list(MoreThan3,is_collinear),
        #                                    reject_max_iter = 100, on_reject = "stop",transform=transforms_pois, 
        #                                    names_final=paste("x", 1:10, sep="") )
        simdata<-simulate_data(sim_design_pois, n_obs=1000000)
        simdata<- simdata[,1:ncov]
        FUT <- round(runif(1000000, 1E-2,2),2)
        simdata$FUT<- FUT
        linpred<- as.matrix(simdata[,1:ncov]) %*% betaV + Intercept+log(FUT)     
        simdata$y <- rpois(length(linpred), exp(linpred))
      }
      my<-mean(simdata$y)
      return((my-mu)**2)
    }
    # find Intercept for an outcome event rate of approximately 10%
    op<-optimize(tar_fun, interval=c(-20, 0),mu=0.1, tol=0.005)
    Intercept <- op$minimum
    print(c("ncov",ncov,Beta1,"Intercept",Intercept))
    write.table(t(c(ncov,Beta1,Intercept)),Interceptfile,sep=",",row.names=FALSE,col.names = FALSE,append=TRUE)
    
    for(EPV in c(3, 10)){
      nobs= EPV*ncov/0.1
      bigsimdata <- matrix(NA, nsim*nobs, 1+ncov+2)
      colnames(bigsimdata)<-c("isim", paste("x", 1:ncov, sep=""), "FUT","y")
      file1 = paste(SimDataDir,"Simdata_FUT_",ncov,"_",nbeta,"_",nobs,".csv",sep="")
      #Now simulate with this design, e.g. with n_obs=nobs
      path1=paste(SimDataDir,"Data_",ncov,"_",round(Beta1,2),"_",nobs,sep="")
      dir.create(path1)
      succ=totEvent=RejSim=0
      rejflag<-c()
      while(succ<nsim){
        rejflag=0
        simdata<-simulate_data(sim_design_pois, n_obs=nobs)
        
        simdata<- simdata[,1:ncov]
        FUT <- round(runif(nobs, 1E-2,2),2)
        simdata$FUT<- FUT
        linpred<- as.matrix(simdata[,1:ncov]) %*% betaV + Intercept+log(FUT)     ### for the 2 covariates designs!
        simdata$y <- rpois(length(linpred), exp(linpred))
        if(any(is.na(simdata$y)==TRUE)){
          RejSim=RejSim+1
          rejflag=1
        }else{
          if((sum(simdata$y)>ncov) & (MoreThan3(simdata[,1:ncov])< 1) &  !(is_collinear(simdata[,1:ncov])) ){
            
            succ = succ+1
            
#            file1 = paste(SimDataDir,"Data_",ncov,"_",round(Beta1,2),"_",nobs,"/Simdata_",ncov,"_",round(Beta1,2),"_",nobs,"_",succ,".csv",sep="")
#            write.table(simdata,file1,sep=",",row.names=FALSE)
            bigsimdata[(succ-1)*nobs+(1:nobs),]<-cbind(isim=succ, as.matrix(simdata[,1:ncov]), 
                                                       FUT=as.numeric(simdata$FUT), y=as.numeric(simdata[,"y"]))
            totEvent<- totEvent+sum(simdata$y)
          }else{
            
            RejSim=RejSim+1
            rejflag=1
          }
        }
        
      }
      print(c("Avg. number of Events ",totEvent/succ))
      print(c("No. of Rejected simulation ",RejSim))
      write.table(bigsimdata,file1,sep=",",row.names=FALSE)
      
    }
  }
  
}


