### This R code counts number of non-convergence in MLE computation
### Specify working directory where MLE outputs are saved as 'workdir'. All the output files will be saved in this directory.
### Author: Ashwini Joshi
### Date: 25-Sep-2019
####################################################################################################################
library("data.table")

workdir="<Path>"
nsim =10000 # specify number of simulations

Nummethod="GLM"

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
EstMethod="MLE"

SeparationFile <- paste(workdir,"SeparationFile_",EstMethod,"_",Nummethod,"_.csv",sep="")
write.table(t(c("ncov","Beta1","nobs","Number of Separations")), SeparationFile,sep=",",append = FALSE,col.names = FALSE,row.names=FALSE)

# for loop of number of covariates
for(ncov in c(2,5,10)){

  # for loop of Beta1 values
  for(Beta1 in c(1:5)){
    if(Beta1==1){betaV[1]=-log(16)}else{if(Beta1==2){betaV[1]= -log(2)}else{if(Beta1==3){betaV[1]= 0}else{if(Beta1==4){betaV[1]= log(2)}else{betaV[1]= log(16)}}}}
    betaV=betaV[1:ncov]
    
    # for loop of events per variable
    for(EPV in c(3,5,10) ){
      nobs= EPV*ncov/0.1
      
      file1 = paste(workdir,"Coeff_",EstMethod,"_",Nummethod,"_",ncov,"_",round(betaV[1],2),"_",nobs,".csv",sep="")
      data1<- data.frame(read.csv(file1,sep=","))
      data2<-data1[, (5+ncov+1):(6+ncov+ncov)]
      maxarray<-apply(data2,1,max)
      sepno<- length( which(maxarray>10))
      write.table(t(c(ncov,betaV[1], nobs,sepno )),SeparationFile,append=TRUE,sep=",",row.names = FALSE,col.names = FALSE)
    }
  }
}
