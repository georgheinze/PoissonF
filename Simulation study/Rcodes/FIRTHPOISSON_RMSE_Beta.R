
### This R code perform RMSE computation for 4 methods - MLE,  Firth corrected MLE (using modified score function) and Firth corrected MLE (using augmented data)
### Directory of input data is 'SimDataDir'. Please specify the path where you have saved the simulation data.
### Also specify working directory as 'workdir'. All the output files will be saved in this directory.
### Author: Ashwini Joshi
### Date: 25-Sep-2019
####################################################################################################################
workdir= "<Path>"
SimDataDir = "<Path>"
setwd(workdir)

Nummethod="GLM"
EstMethod= "MLE"  ## Method can be one of "MLE", "Firth","ModFirth"

nsim =10000 # specify number of simulations

theta<-c()
theta[1] =log(16)
theta[2]<-theta[3]<-theta[4]<-0.69
theta[5]<-theta[6]<-0.345
theta[7]<-0.036481
theta[8]<-0.0029559
theta[9]<-0.0038941
theta[10]<-0.036481

betaV<-numeric(10)
for(i in 1:10){
  betaV[i]<-theta[i]*(-1)^i
}
RMSEFilenm<- c()

IntcFile= read.csv(paste(SimDataDir,"Intercept.csv",sep=""),header= TRUE)

# for loop of number of covariates
for (ncov in c(2, 5, 10)){

  # for loop of Beta1 values
    for(Beta1 in c(1:5)){
      if(Beta1==1){betaV[1]=-log(16)}else{if(Beta1==2){betaV[1]= -log(2)}else{if(Beta1==3){betaV[1]= 0}else{if(Beta1==4){betaV[1]= log(2)}else{betaV[1]= log(16)}}}}
      
      betaV=betaV[1:ncov]
      IntcVal=IntcFile[which(IntcFile$ncov== ncov & (abs(round(IntcFile$Beta1,2)-betaV[1])<1E-2)),"Intercept"]
      
      # for loop of events per variable
      for(EPV in c(3, 5, 10)){
        
        nobs= EPV*ncov/0.1
        Est_BetaMLE <-Var_BetaMLE <- Bias_BetaMLE <- RMSE_BetaMLE <- c()
        
        # Create column header array
        colHead<-ans<-c()
        colHead=c("Num.ofCov","Beta1","nobs","Beta0", paste("Beta0_", EstMethod,"_RMSE" ,sep=""), paste("Beta0_", EstMethod,"_Var" ,sep=""), paste("Beta0_", EstMethod,"_Bias" ,sep="") )
   
        CoefData<- read.csv(paste(workdir,"Output/Coeff_",EstMethod,"_",Nummethod,"_",ncov,"_",round(betaV[1],2),"_",nobs,".csv",sep=""),sep=",",header = TRUE)
        Est_BetaMLE<- sapply(CoefData[,4:(4+ncov)],mean)
        Var_BetaMLE <-sapply(CoefData[,4:(4+ncov)],var)
        
        Bias_BetaMLE<- Est_BetaMLE- as.numeric( c(IntcVal,betaV))
        RMSE_BetaMLE<- sqrt(Var_BetaMLE *(nsim-1)/nsim + (Bias_BetaMLE)^2)
        
        ans=c(ncov,betaV[1],nobs,IntcVal,RMSE_BetaMLE[1],Var_BetaMLE[1],Bias_BetaMLE[1])
        for (i in 2:(ncov+1)){
          ans=c(ans,betaV[i-1],RMSE_BetaMLE[i],Var_BetaMLE[i],Bias_BetaMLE[i])
          colHead=c(colHead,paste("Beta",(i-1),sep=""),paste("Beta",(i-1),"_", EstMethod,"_RMSE" ,sep="") ,paste("Beta",(i-1),"_", EstMethod,"_Var" ,sep="") ,paste("Beta",(i-1),"_", EstMethod,"_Bias" ,sep="")  )
        }
        RMSEFilenm<- paste(workdir,"Output/", EstMethod, "_RMSE_Beta_",ncov,".csv",sep="")
        if(file.exists(RMSEFilenm)){
          write.table(t(ans),RMSEFilenm,sep=",",append = TRUE,col.names = FALSE,row.names = FALSE)
        }else{
          write.table(t(colHead),RMSEFilenm,sep=",",append = FALSE,col.names = FALSE,row.names = FALSE)
          write.table(t(ans),RMSEFilenm,sep=",",append = TRUE,col.names = FALSE,row.names = FALSE)
        }      
    }# for loop of EPV ends
  }# for loop of Beta1 ends
}# for loop of ncov ends



