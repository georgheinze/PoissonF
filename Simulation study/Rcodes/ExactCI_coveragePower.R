
### This R code reads exact and mid-p CI values from LogXact output.
### Specify working directory where Logxact outputs are saved as 'workdir'. All the output files will be saved in this directory.
### Author: Ashwini Joshi
### Date: 25-Sep-2019
####################################################################################################################

library("MASS")
library("stringr")
library("parallel")
library("data.table")

workdir= "<Path>"
setwd(workdir)

Nummethod="GLM"
EstMethod="Exact"  

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

# for loop for number of covariates
for (ncov in c(2,5)){
    write.table(t(c("ncov","nbos","Beta1","Beta2",paste0(EstMethod,"CI"),"PowerEx1", "ExCoverage1", "LowExCoverage1","UpExCoverage1","PowerEx2", "ExCoverage2", "LowExCoverage2","UpExCoverage2",
                    paste0(EstMethod,"MidpCI"),"PowerMidp1", "MidpCoverage1", "LowMidpCoverage1","UpMidpCoverage1", "PowerMidp2", "MidpCoverage2", "LowMidpCoverage2","UpMidpCoverage2")),
                paste(workdir,"ExactCICoverage_",ncov,".csv",sep=""),sep=",",append = FALSE,col.names = FALSE,row.names=FALSE)
  
  # for loop for Beta1 values
  for(Beta1 in c(1:5)){
    if(Beta1==1){betaV[1]=-log(16)}else{if(Beta1==2){betaV[1]= -log(2)}else{if(Beta1==3){betaV[1]= 0}else{if(Beta1==4){betaV[1]= log(2)}else{betaV[1]= log(16)}}}}
    
    betaV=betaV[1:ncov]
    
    # for loop for events per variable
    for(EPV in c(3, 5)){
      nobs= EPV*ncov/0.1
     
          file1 =paste0("ExCI_",ncov,"_",round(betaV[1],2),"_",nobs,".csv")
          file1=paste(workdir,file1,sep="")
         
          if (! file.exists(file1))
             stop(paste(file1, "does not exist"))
          
          data1<- read.csv(file1,sep=",")
          
          PowerExact1<- CoverageEx1<- LowCoverageEx1<-UpCoverageEx1<-PowerEx2<- CoverageEx2<- LowCoverageEx2<-UpCoverageEx2<-0
          PowerMidp1<- CoverageMidp1<- LowCoverageMidp1<-UpCoverageMidp1<-PowerMidp2<- CoverageMidp2<- LowCoverageMidp2<-UpCoverageMidp2<-0

          # Power computation for Beta1  Exact and Midp CI
          PowerExact1<- round((nsim-length(which(data1[,5]<=0 & data1[,6]>=0)))/10,0)
          PowerMidp1<- round((nsim-length(which(data1[,7]<=0 & data1[,8]>=0)))/10,0)
          
          # Power computation for Beta2  Exact and Midp CI
          PowerExact2<- round((nsim-length(which(data1[,12]<=0 & data1[,13]>=0)))/10,0)
          PowerMidp2<- round((nsim-length(which(data1[,14]<=0 & data1[,15]>=0)))/10,0)
          
          # CI coverage for Beta1  Exact and Midp CI
          CoverageEx1<- round(length(which(data1[,5]<=betaV[1] & data1[,6]>=betaV[1]))/10,0)
          CoverageMidp1<- round(length(which(data1[,7]<=betaV[1] & data1[,8]>=betaV[1]))/10,0)
          
          LowCoverageEx1<- round(length(which(data1[,5]<=betaV[1]))/10,0)
          LowCoverageMidp1<- round(length(which(data1[,7]<=betaV[1] ))/10,0)
          
          UpCoverageEx1<- round(length(which(data1[,6]>=betaV[1]))/10,0)
          UpCoverageMidp1<- round(length(which( data1[,8]>=betaV[1]))/10,0)
          
          # CI coverage for Beta2  Exact and Midp CI
          CoverageEx2<- round(length(which(data1[,12]<=betaV[2] & data1[,13]>=betaV[2]))/10,0)
          CoverageMidp2<- round(length(which(data1[,14]<=betaV[2] & data1[,15]>=betaV[2]))/10,0)
          
          LowCoverageEx2<- round(length(which(data1[,12]<=betaV[2]))/10,0)
          LowCoverageMidp2<- round(length(which(data1[,14]<=betaV[2] ))/10,0)
          
          UpCoverageEx2<- round(length(which(data1[,13]>=betaV[2]))/10,0)
          UpCoverageMidp2<- round(length(which( data1[,15]>=betaV[2]))/10,0)
         
          CIVal<- c(ncov,nobs,betaV[1],betaV[2],"",PowerExact1, CoverageEx1, LowCoverageEx1,UpCoverageEx1, PowerExact2, CoverageEx2, LowCoverageEx2,UpCoverageEx2,"", 
                    PowerMidp1, CoverageMidp1, LowCoverageMidp1,UpCoverageMidp1,PowerMidp2, CoverageMidp2, LowCoverageMidp2,UpCoverageMidp2)
        
          write.table(t(CIVal),paste(workdir,"ExactCICoverage_",ncov,".csv",sep=""),sep=",",append = TRUE,col.names = FALSE, row.names=FALSE)

      }# for loop of Beta1 ends
    }# for loop of EPV ends
  }# for loop of ncov ends
  
#}
