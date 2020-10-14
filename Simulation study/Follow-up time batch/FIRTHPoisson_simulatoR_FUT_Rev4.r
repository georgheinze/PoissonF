### need to be rerun (mu_FUT was actually mu - but does not affect SAS calculations as long as adequately handled)
### mu_FUT could be simply computed as mu*FUT
### Follow-up 'time' FUT is actually simulated from a truncated (positive) Poisson distribution to better match the data example

### This R code generates simulated data for 91 combinations of number of covariates (3), values of beta1 (9), sample sizes (3) using SimulatoR package
### Part 1 is simulation setup with SimulatoR
### Part 2 performs simulations for 2,5,10 covariates, 9 values of beta1, 3 different sample sizes(for 3 values of Expected number of events per variable: 3, 5, 10)
### Simulations are performed for data with follow-up time.
### Simulation data will be saved in SimDataDir. Please specify the path where you want to save the simulation data.

library(simdata)
library(extraDistr)
library(truncnorm)

githubpath = "E:/1Projects/2018-08 NegBinom FC/Github/"
SimDataDir =  paste(githubpath,"PoissonF/Simulation study/Follow-up time batch/Data/", sep="")

################################ Part 1 #######################################
transforms_pois<-function_list(v1 = function(z) (z[,1] > 1.28),    # change to Puhr, to achieve balance of 10:90
                      v2 = function(z) (z[,2] > 0.35),
                      v3 = function(z) (z[,3] > 0),
                      v4 = function(z) (z[,4] > 0),
                      v5 = function(z) (z[,5] >= -1.2) + (z[,5] >= 0.75),
                      v6 = function(z) (z[,6] >= 0.5)+(z[,6] >= 1.5),
                      v7 = function(z) floor(10*z[,7]+55),
                      v8 = function(z) floor(pmax(0, 100*exp(qtruncnorm(pnorm(z[,8]), b=qnorm(0.99)))-20)),
                      v9 = function(z) floor(pmax(0, 2*exp(qtruncnorm(pnorm(z[,9]), b=qnorm(0.99)))-1)),
                      v10 = function(z) floor(10*z[,10]+120)
)
relations_pois = cor_from_upper(10, 
                                rbind(c(1,2,0.5), c(1,3,0.5), c(1,7, 0.5),
                                      c(3,4,-0.5), c(3,5,-0.3),
                                      c(4,5,0.5), c(4,7,0.3),
                                      c(4,8,0.5), c(4,9,0.3),
                                      c(5,8,0.3), c(5,9,0.3),
                                      c(6,7,-0.3), c(6,8,0.3),
                                      c(8,9,0.5), c(7,10,0.5)))
sim_design_pois = mvtnorm_simdesign(relations=relations_pois, transform=transforms_pois, 
                         truncate=c(NA, NA, NA, NA, NA, NA, 5, 5, 10, 5),
                         names_final=paste("x", 1:10, sep=""))
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

## define betas


sext.diff<-function(x) diff(quantile(x,c(1/6,5/6)))


theta<-numeric(10)
names(theta)<-colnames(simdata)
theta[1]<-0

# binary covariates

theta[2]<-0.69
theta[3]<- -0.69
theta[4]<-0.69

# ordinal covariates (0,1,2)

theta[5] <- -0.345
theta[6]<-0.345

# continuous covariates

theta[7]<- -log(2)/sext.diff(simdata$x7)
theta[8]<-log(2)/sext.diff(simdata$x8)
theta[9]<- -log(2)/sext.diff(simdata$x9)
theta[10]<-log(2)/sext.diff(simdata$x10)

theta[2:10]<-round(theta[2:10],3)


################################ Part 2 #######################################
# start simulations
nsim =10000 # specify number of simulations
RNGkind(kind ="Mersenne-Twister", normal.kind ="Inversion", sample.kind ="Rejection" )
set.seed(79123)
RejSim =0
Interceptfile = paste(SimDataDir,"Intercept.csv",sep="")
write.table(t(c("ncov","Beta1","nbeta","nobs","nsim","Intercept","Events", "Events(FUT)", "RejSim")),Interceptfile,sep=",",row.names=FALSE, col.names = FALSE)

#for (ncov in c(2, 5, 10)){
#for (ncov in c(5, 10)){  
for (ncov in c(10)){  
  nbeta =0 
  for(Beta1 in c(-log(16),-log(8), -log(4), -log(2),0, log(2), log(4), log(8),log(16))){
    nbeta = nbeta +1
    betaV<-rep(0,ncov)
    a<-1   # in Puhr2017, was set to 0, 0.5, 1
    
    betaV[1] = Beta1
    betaV[2:ncov]<-theta[2:ncov] * a
    
    #Target function to determine intercept parameter:
    tar_fun<-function(Intercept, mu=0.1) {
      simdata<-simulate_data(sim_design_pois, n_obs=1000000)
      simdata<- simdata[,1:ncov]
#      FUT <- round(runif(1000000, 1E-2,2),2)
      FUT <-rtpois(1000000, 1.6, a=0, b=Inf)
      simdata$FUT<- FUT
      linpred<- as.matrix(simdata[,1:ncov]) %*% betaV + Intercept+log(FUT)    
      simdata$y <- rpois(length(linpred), exp(linpred))
      while( any(is.na(simdata$y)==TRUE)){
        # simdata<-conditional_simulate_data(relations_pois, n_obs=1000000, reject = list(MoreThan3,is_collinear),
        #                                    reject_max_iter = 100, on_reject = "stop",transform=transforms_pois, 
        #                                    names_final=paste("x", 1:10, sep="") )
        simdata<-simulate_data(sim_design_pois, n_obs=1000000)
        simdata<- simdata[,1:ncov]
#        FUT <- round(runif(1000000, 1E-2,2),2)
        FUT <-rtpois(1000000, 1.6, a=0, b=Inf)
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
#    print(c("ncov",ncov,Beta1,"Intercept",Intercept))
    
    for(EPV in c(3, 5, 10)){
      nobs= EPV*ncov*10
      bigsimdata <- matrix(NA, nsim*nobs, 1+ncov+5)
      colnames(bigsimdata)<-c("isim", paste("x", 1:ncov, sep=""), "FUT","y", "mu", "y_FUT", "mu_FUT")
      file1 = paste(SimDataDir,"Simdata_FUT_",ncov,"_",nbeta,"_",nobs,".csv",sep="")
      #Now simulate with this design, e.g. with n_obs=nobs
#      path1=paste(SimDataDir,"Data_",ncov,"_",round(Beta1,2),"_",nobs,sep="")
#      dir.create(path1)
      succ=totEvent=totEvent_FUT=RejSim=0
      rejflag<-c()
      while(succ<nsim){
        rejflag=0
        simdata<-simulate_data(sim_design_pois, n_obs=nobs)
        
        simdata<- simdata[,1:ncov]
#        FUT <- round(runif(nobs, 1E-2,2),2)
        FUT <-rtpois(nobs, 1.6, a=0, b=Inf)
        simdata$FUT<- FUT
        linpred<- as.matrix(simdata[,1:ncov]) %*% betaV + Intercept    
        linpred_FUT<- linpred + log(FUT)
        simdata$y <- rpois(length(linpred), exp(linpred))
        simdata$y_FUT <- rpois(length(linpred), exp(linpred_FUT))
        if(any(is.na(simdata$y_FUT)==TRUE)){
          RejSim=RejSim+1
          rejflag=1
        }else{
          if((sum(simdata$y)>ncov) & (sum(simdata$y_FUT)>ncov) & (MoreThan3(simdata[,1:ncov])< 1) &  !(is_collinear(simdata[,1:ncov])) ){
            
            succ = succ+1
            
#            file1 = paste(SimDataDir,"Data_",ncov,"_",round(Beta1,2),"_",nobs,"/Simdata_",ncov,"_",round(Beta1,2),"_",nobs,"_",succ,".csv",sep="")
#            write.table(simdata,file1,sep=",",row.names=FALSE)
            bigsimdata[(succ-1)*nobs+(1:nobs),]<-cbind(isim=succ, as.matrix(simdata[,1:ncov]), 
                                                       FUT=as.numeric(simdata$FUT), y=as.numeric(simdata[,"y"]), mu=exp(linpred), y_FUT=as.numeric(simdata[,"y_FUT"]), mu_FUT=exp(linpred_FUT))
            totEvent<- totEvent+sum(simdata$y)
            totEvent_FUT<- totEvent_FUT+sum(simdata$y_FUT)
          }else{
            
            RejSim=RejSim+1
            rejflag=1
          }
        }
        
      }
      print(cbind(ncov,Beta1,nobs,nsim,Intercept,events_per_dataset=totEvent/succ,events_with_FUT=totEvent_FUT/succ, RejSim))
      write.table(t(c(ncov,Beta1,nobs,nsim,Intercept,totEvent/succ, totEvent_FUT/succ, RejSim)),Interceptfile,sep=",",row.names=FALSE,col.names = FALSE,append=TRUE)
      write.table(bigsimdata,file1,sep=",",row.names=FALSE)
      
    }
  }
  
}


