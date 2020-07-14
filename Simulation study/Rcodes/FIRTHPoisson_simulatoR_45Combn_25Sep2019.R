
### This R code perform simulations for 45 combinations of number of covariates, values of beta1, sample sizes using SimulatoR package
### Part 1 is simulation setup with SimulatoR
### Part 2 performs simulations for 3 different number of covariates, 5 values of beta1, 3 different sample sizes(for 3 Expected number of events)
### Simulation data will be saved in SimDataDir. Please specify the path where you want to save the simulation data.
### Author: Ashwini Joshi
### Date: 25-Sep-2019
####################################################################################

library("SimulatoR")
SimDataDir = "<Path>"

################################ Part 1 #######################################
transforms_pois<-list(v1 = function(z) (z[,1] > 1.28),    # change to Puhr, to achieve balance of 10:90
                      v2 = function(z) (z[,2] > 0.52),    # change  to achieve balance of 30:70
                      v3 = function(z) (z[,3] > 0),
                      v4 = function(z) (z[,4] > 0),
                      v5 = function(z) (z[,5] >= -1.2) + (z[,5] >= 0.75),
                      v6 = function(z) (z[,6] >= 0.5)+(z[,6] >= 1.5),
                      v7 = function(z) floor(10*z[,7]+55),
                      v8 = function(z) floor(pmax(0, 100*exp(z[,8])-20)),
                      v9 = function(z) floor(pmax(0, 80*exp(z[,9])-20)),
                      v10 = function(z) floor(10*z[,10]+55)
)
relations_pois = cor_from_upper(10, 
                                rbind(c(1,2,0.5), c(1,3,0.5), c(1,7, 0.5),
                                      c(3,4,-0.5), c(3,5,-0.3),
                                      c(4,5,0.5), c(4,7,0.3),
                                      c(4,8,0.5), c(4,9,0.3),
                                      c(5,8,0.3), c(5,9,0.3),
                                      c(6,7,-0.3), c(6,8,0.3),
                                      c(8,9,0.5)))
sim_design_pois = design(relations=relations_pois, transform=transforms_pois, 
                         truncate=c(NA, NA, NA, NA, NA, NA, 5, 5, 5, 5),
                         names_final=paste("x", 1:10, sep=""))

################################ Function to specify conditions on simulated data ##########
f1<-function(x){return(length(which(x!=0)))}
f2<-function(x){return(length(unique(x)))}

MoreThan3<-function (x) 
{
  x = as.matrix(x)
  val1 = any(apply(x, 2, f1) < 3)
  val2 = any(apply(x, 2, f2) < 2)
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

plot_initial_cor_network(sim_design_pois)   # plots network for the normal deviates
plot_final_cor_network(sim_design_pois)  # plots network for the simulated 'real' covariates

# Specify true values of theta
theta<-numeric(10)
names(theta)<-colnames(simdata)
theta[1]<-log(16)
theta[2]<-theta[3]<-theta[4]<-0.69
theta[5]<-theta[6]<-0.345

sext.diff<-function(x) diff(quantile(x,c(1/6,5/6)))

theta[7]<-log(2)/sext.diff(simdata$x7)
theta[8]<-log(2)/sext.diff(simdata$x8)
theta[9]<-log(2)/sext.diff(simdata$x9)
theta[10]<-log(2)/sext.diff(simdata$x10)

################################ Part 2 #######################################
# start simulations
nsim =10000 # specify number of simulations

set.seed(79123)
RejSim =0 # Initialize number of rejected simulation to 0
Interceptfile = paste(SimDataDir,"Intercept.csv",sep="")
write.table(t(c("ncov","Beta1","Intercept")),Interceptfile,sep=",",row.names=FALSE, col.names = FALSE)

# For loop of number of covariates
for (ncov in c(2, 5, 10)){
  B1cnt =0 
  # For loop of Beta1 values 
  for(Beta1 in c(-log(16),-log(2),0, log(2),log(16))){
    B1cnt = B1cnt +1
    betaV<-numeric(10)
    a<-1   # in Puhr2017, was set to 0, 0.5, 1
    for(i in 1:10){
      betaV[i]<-theta[i]*(-1)^i # Every alternate beta value is negative
    }
    betaV[1] = Beta1
    betaV=betaV[1:ncov]
    
    #Target function to determine intercept parameter:
    tar_fun<-function(Intercept, mu=0.1) {
      
      simdata<-simulate_data(sim_design_pois, n_obs=1000000)
      simdata<- simdata[,1:ncov]
      linpred<- as.matrix(simdata) %*% betaV + Intercept    
      simdata$y <- rpois(length(linpred), exp(linpred))
      
      while( any(is.na(simdata$y)==TRUE)){
        simdata<-simulate_data(sim_design_pois, n_obs=1000000)
        simdata<- simdata[,1:ncov]
        linpred<- as.matrix(simdata) %*% betaV + Intercept    
        simdata$y <- rpois(length(linpred), exp(linpred))
      }
      my<-mean(simdata$y)
      return((my-mu)**2)
    }
    # find Intercept for an outcome event rate of approximately 10%
    op<-optimize(tar_fun, interval=c(-10, 10),mu=0.1, tol=0.005)
    Intercept <- op$minimum
    # print(c("ncov",ncov,"Beta1",B1cnt,"Intercept",Intercept))
    write.table(t(c(ncov,Beta1,Intercept)),Interceptfile,sep=",",row.names=FALSE,col.names =FALSE,append=TRUE)
   
    # For loop of Events per variable
    for(EPV in c(3, 5, 10)){
      nobs= EPV*ncov/0.1
      #Now simulate with this design, e.g. with n_obs=nobs
      path1=paste(SimDataDir,"Data_",ncov,"_",round(Beta1,2),"_",nobs,sep="")
      dir.create(path1)
      succ=totEvent=RejSim=0
      rejflag<-c()
      
      while(succ<nsim){
        rejflag=0
        simdata<-simulate_data(sim_design_pois, n_obs=nobs)
        simdata<-simdata[,1:ncov]
        linpred<- as.matrix(simdata) %*% betaV + Intercept    ### for the 10 covariates designs!
        simdata$y <- rpois(length(linpred), exp(linpred))
        if(any(is.na(simdata$y)==TRUE)){
          RejSim=RejSim+1
          rejflag=1
        }else{
          # If conditions for discarding extreme data
          if(((sum(simdata$y)/nobs) >0.01) & ((sum(simdata$y)/nobs) < 2) & (MoreThan3(simdata[,1:ncov])< 1) &  !(is_collinear(simdata[,1:ncov]))){
            succ = succ+1
            # Specify filename to save the simulation data
            file1 = paste(SimDataDir,"Data_",ncov,"_",round(Beta1,2),"_",nobs,"/Simdata_",ncov,"_",round(Beta1,2),"_",nobs,"_",succ,".csv",sep="")
            write.table(simdata,file1,sep=",",row.names=FALSE)
            totEvent<- totEvent+sum(simdata$y)
          }else{
            RejSim=RejSim+1
            rejflag=1
          }
        }
      }
      print(c("Avg. number of Events ",totEvent/succ))
      print(c("No. of Rejected simulation ",RejSim))
    }
  }
}




