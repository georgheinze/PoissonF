
### This R code perform RMSE computation for 3 methods - MLE, Firth corrected MLE (using modified score function) and Modified Firth corrected MLE (using augmented data)
### Directory of input data is 'SimDataDir'. Please specify the path where you have saved the simulation data.
### Also specify working directory as 'workdir'. All the output files will be saved in this directory.
### Author: Ashwini Joshi
### Date: 25-Sep-2019
####################################################################################################################

workdir= "<Path>"
SimDataDir = "<Path>"
setwd(workdir)

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

# Create a file to save RMSE values
write.table(t(c("Num.ofCov","Beta1","nobs","MLE_RMSE","Firth_RMSE","ModFirth_RMSE")),paste(workdir,"64CoreOutput/RMSE_",EstMethod,"_Prediction.csv",sep=""),sep=",",append = FALSE,col.names = FALSE,row.names = FALSE)

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
      MLE_RMSE =Firth_RMSE = ModFirth_RMSE = MLEFirth_RMSE = 0
      CoefData<- read.csv(paste(workdir,"64CoreOutput/Coeff_",EstMethod,"_",Nummethod,"_",ncov,"_",round(betaV[1],2),"_",nobs,".csv",sep=""),sep=",",header = TRUE)
      lapply(CoefData[,4:(4+ncov)],mean)
      
      # for loop for simulations 
      for(succ in 1:nsim){
        file1 = paste(SimDataDir,"Data_",ncov,"_",round(betaV[1],2),"_",nobs,"/Simdata_",ncov,"_",round(betaV[1],2),"_",nobs,"_",succ,".csv",sep="")
        data1<- read.csv(file1,sep=",")
        if (ncov==10){
          for(scaling in c(1:4)){
            data1[,(6+scaling)]= (data1[,(6+scaling)]/sd( data1[,(6+scaling)]))
          }
        }
        
        MLECoef<- CoefData[which(CoefData$SimID==succ),4:(4+ncov)]
        MLE_RMSE= MLE_RMSE + sum((exp(as.matrix(cbind(rep(1,dim(data1)[1]),data1[,-(ncov+1)]))%*%t(as.matrix(MLECoef)))-exp(as.matrix(cbind(rep(1,dim(data1)[1]),data1[,-(ncov+1)]))%*%as.matrix(c(IntcVal,betaV))))^2)
        
        FirthCoef<- CoefData[which(CoefData$SimID==succ),(8+2*ncov):(8+3*ncov)]
        Firth_RMSE= Firth_RMSE + sum((exp(as.matrix(cbind(rep(1,dim(data1)[1]),data1[,-(ncov+1)]))%*%t(as.matrix(FirthCoef)))-exp(as.matrix(cbind(rep(1,dim(data1)[1]),data1[,-(ncov+1)]))%*%as.matrix(c(IntcVal,betaV))))^2)

        ModFirthCoef<- CoefData[which(CoefData$SimID==succ),(10+3*ncov):(10+4*ncov)]
        ModFirth_RMSE= ModFirth_RMSE + sum((exp(as.matrix(cbind(rep(1,dim(data1)[1]),data1[,-(ncov+1)]))%*%t(as.matrix(ModFirthCoef)))-exp(as.matrix(cbind(rep(1,dim(data1)[1]),data1[,-(ncov+1)]))%*%as.matrix(c(IntcVal,betaV))))^2)
        
      }# For loop of simulation ends
      MLE_RMSE = sqrt(MLE_RMSE/(nobs*nsim))
      Firth_RMSE = sqrt(Firth_RMSE/(nobs*nsim))
      ModFirth_RMSE = sqrt(ModFirth_RMSE/(nobs*nsim))
      
      write.table(t(c(ncov,round(betaV[1],2),nobs,MLE_RMSE,Firth_RMSE,ModFirth_RMSE)),paste(workdir,"Output/RMSE_Prediction.csv",sep=""),sep=",",append = TRUE,col.names = FALSE,row.names = FALSE)
     
    }# for loop of EPV ends
  }# for loop of Beta1 ends
}# for loop of ncov ends



