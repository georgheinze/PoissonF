
### This R code perform MLE computation, Firth corrected MLE (using modified score function) and Modified Firth corrected MLE (using augmented data)
### It saves coefficients obtained from these three estimation methods in 45 files depending upon the combination of number of covariates, values of beta1, sample sizes.
### Directory of input data is 'SimDataDir'. Please specify the path where you have saved the simulation data.
### Also specify working directory as 'workdir'. All the output files will be saved in 'Output' subdirectory of the working directory.
### For MLE method, regression coefficients, standard errors and wald CI are computed.
### For Firth method, regression coefficients, standard errors and penalized profile likelihood CI are computed.
### For Modified Firth method, regression coefficients and standard errors are computed.
### Author: Ashwini Joshi
### Date: 25-Sep-2019
####################################################################################################################
library("MASS")
library("stringr")
library("parallel")
library("data.table")
library("likelihoodAsy")
library("brglm2")

workdir="<Path>" # Specify the working directory. Create 'Output' folder inside the directory to save all the outputs.
SimDataDir = "<Path>" # Specify the directory/folder where simulated data is saved.

setwd(workdir)

Nummethod="GLM" 

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

EstMethod="Firth" # Method can be "MLE", "Firth" or "ModFirth"
#EvalEst<- function(workdir,SimDataDir, betaV,EstMethod,CovList, EPVList, Beta1List ){
  # for loop of number of covariates
  for (ncov in CovList){
    
  	# for loop for number of events per variable
    for(EPV in EPVList){
      nobs= EPV*ncov/0.1
      
    	# for loop of Beta1 values
    	for(Beta1 in Beta1List){
    		if(Beta1==1){betaV[1]=-log(16)}else{if(Beta1==2){betaV[1]= -log(2)}else{if(Beta1==3){betaV[1]= 0}else{if(Beta1==4){betaV[1]= log(2)}else{betaV[1]= log(16)}}}}
    
    		betaVnew=betaV[1:ncov]
  
    		
    		write.table(t(c("SimID","Simdata",paste0(EstMethod,"Error"),paste0(EstMethod, 0:ncov),paste0(EstMethod,"SE"),paste0(EstMethod, 0:ncov),
    		                paste0(EstMethod,"WaldCI"),"LowWaldCI1", "UpWaldCI1","LowWaldCI2", "UpWaldCI2","PowerWald1", "WaldCoverage1", "LowWaldCoverage1","UpWaldCoverage1","PowerWald2", "WaldCoverage2", "LowWaldCoverage2","UpWaldCoverage2",
    		                paste0(EstMethod,"PLCI"),"PowerPL1", "PLCoverage1", "LowPLCoverage1","UpPLCoverage1", "PowerPL2", "PLCoverage2", "LowPLCoverage2","UpPLCoverage2")),
    		            paste(workdir,"Output/Coeff_",EstMethod,"_",Nummethod,"_",ncov,"_",round(betaV[1],2),"_",nobs,".csv",sep=""),sep=",",append = FALSE,col.names = FALSE,row.names=FALSE)
    		
        SimDataDir1=paste(SimDataDir,"Data_",ncov,"_",round(betaV[1],2),"_",nobs,"/",sep="")
        filenames=list.files(path=SimDataDir1)
       
        # Function to compute regression coefficients, standard error and CI.
        # This function can be called in parallel threads.
        SimCompute<- function(succ,filenames,EstMethod, Nummethod){
          
          file1=filenames[succ]
          file1=paste(SimDataDir1,file1,sep="")
         
          if (! file.exists(file1))
             stop(paste(file1, "does not exist"))
          
          data1<- as.data.frame (read.csv(file1,sep=","))
          
          if (ncov==10){
            for(scaling in c(1:4)){ 
              data1[,(6+scaling)]= (data1[,(6+scaling)]/sd( data1[,(6+scaling)])) # To scale the continuous variables.
            }
          }
         
          formula1<- paste("y~",paste(paste0("x", 1:ncov), collapse= "+"))
         
          # Estimation #########################################################################################
          MLECoef <- MLESE <- confcoef <- MLE_WaldCI <- maxloglkhd <- c()
          PowerWald1<- CoverageWald1<- LowCoverageWald1<-UpCoverageWald1<-PowerWald2<- CoverageWald2<- LowCoverageWald2<-UpCoverageWald2<-0
          PowerPL1<- CoveragePL1<- LowCoveragePL1<-UpCoveragePL1<-PowerPL2<- CoveragePL2<- LowCoveragePL2<-UpCoveragePL2<-0

          .hat2<-function(beta,data1){
            x<- data1[1:ncov]
            x<- as.matrix(cbind(rep(1,nobs),x))
            w <- rep(1,nobs)
            mu<-exp(x %*% beta)*w
            
            mumat<-matrix(0,length(mu), length(mu))
            
            for(i in 1:length(mu)) mumat[i,i]<-mu[i]
            if(is.finite(det(t(x) %*% mumat %*% x))){
              wmat<- solve(t(x) %*% mumat %*% x)
              #  wmat<- ginv(t(x) %*% mumat %*% x)
              hatmatdiag<- c()
              hatmatdiag =rep(0,length(mu))
              for ( i in 1:dim(x)[2]){
                for ( j in 1:dim(x)[2]){
                  hatmatdiag=hatmatdiag+(mu*x[,i]*x[,j]*wmat[i,j])
                }
              }
              return(hatmatdiag)
            }else{return(rep(Inf,length(mu)))}
          }
          
          if(EstMethod=="MLE"){
            res_Est<-glm(formula=formula1, data=data1, family = poisson(link= "log"),control=list(maxit=500))
            
          }else{
            
              if(EstMethod=="Firth"){
                res_Est<- glm(formula=formula1, data=data1, family = poisson(link= "log"),control=list(maxit=500), method = "brglmFit")
              }else{  
                # Estmethod is ModFirth
                x<- data1[1:ncov]
                x<- as.matrix(cbind(rep(1,nobs),x))
                w <- rep(1,nobs)
                y<- data1$y
                
                # Get hat diagonal elements for MLE 
                simid<- substr(filenames[succ],(str_locate_all(filenames[succ], "_")[[1]][length(str_locate_all(filenames[succ], "_")[[1]])] +1),(str_locate_all(filenames[succ], ".c")[[1]][2] -2))
                MLE1<-c()
                aa<- as.vector(CoefData[which(CoefData$SimID== simid),4:(4+ncov)])
                for(i in 1:(ncov+1)){MLE1<- c(MLE1,aa[[i]])}
                InitHat<- .hat2(MLE1, data1)
                
                # create augmented data
                x_aug<- rbind(x,x)
                y_aug<- c(y,0.5*InitHat)
                w_aug<- c(w,w)
                ghost1<- c(rep(0,length(InitHat)),rep(1,length(InitHat)))
                x_aug<- cbind(x_aug,ghost1)
                
                dat_aug<- data.frame(cbind(x_aug,y_aug),row.names = NULL,check.rows = FALSE,check.names = FALSE)
                
                # Create formula to be passed to GLM
                formula1 = as.formula(paste("y_aug~",paste(colnames(x_aug)[-1], collapse= "+")),env = parent.frame())
                res_Est<- glm(formula=formula1, data=dat_aug, family = poisson(link= "log"),control=list(maxit=500))
              }
          }

          # returns log likelhood value
          .likelihood3<-function(beta, data1){
            x<- data1[1:ncov]
            x<- as.matrix(cbind(rep(1,nobs),x))
            w <- rep(1,nobs)
            y<- data1$y
            eta<-x %*% beta
            mu<-exp(eta)*w #let mu denote the rate
            loglike<-sum((y*log(mu)-mu))
            if(is.finite(loglike)){
              if(is.finite(loglike)){
                attr(loglike,"gradient") <- .gradient(beta,data1)
                attr(loglike,"hessian") <- .hessian(beta,data1)
                return(loglike)
              }else{return(Inf)}
            }else{return(Inf)}
          }
          
          # returns gradient
          .gradient<-function(beta, data1){
            x<- data1[1:ncov]
            x<- as.matrix(cbind(rep(1,nobs),x))
            w <- rep(1,nobs)
            y<- data1$y
            mu<-exp(x %*% beta)*w
            hdiag<-rep(0,length(y))
            g<-unlist(sapply(1:length(beta), function(i) sum(((y+hdiag/2)-mu)*x[,i])))
            if(is.na(g)){g<- rep(0,(ncov+1))}
            return(g)
          }
          
          # returns hessian
          .hessian<-function(beta, data1){
            x<- data1[1:ncov]
            x<- as.matrix(cbind(rep(1,nobs),x))
            w <- rep(1,nobs)
            y<- data1$y
            mu<-exp(x %*% beta)*w
            mumat<-matrix(0,length(mu), length(mu))
            for(i in 1:length(mu)) mumat[i,i]<-mu[i]
            return(   -t(x) %*%  mumat %*% x)
          }
          
          # returns  log likelhood with Firth correction
          .likelihood3Firth<-function(beta, data1){
            x<- data1[1:ncov]
            x<- as.matrix(cbind(rep(1,nobs),x))
            w <- rep(1,nobs)
            y<- data1$y
            eta<-x %*% beta
            mu<-exp(eta)*w  #let mu denote the rate
            loglike<- sum(((y+.hat2(beta,data1)/2)*log(mu)-mu))
            if(is.finite(loglike)){
                attr(loglike,"gradient") <- .gradientFirth(beta,data1)
                attr(loglike,"hessian") <- .hessian(beta,data1)
                return(loglike)
            }else{return(Inf)}
          }
          
          # returns  gradient with Firth correction
          .gradientFirth<-function(beta, data1){
            x<- data1[1:ncov]
            x<- as.matrix(cbind(rep(1,nobs),x))
            w <- rep(1,nobs)
            y<- data1$y
            mu<-exp(x %*% beta)*w
            hdiag<- .hat2(beta,data1)
            g<-unlist(sapply(1:length(beta), function(i) sum(((y+hdiag/2)-mu)*x[,i])))
            if(is.na(g)){g<- rep(0,(ncov+1))}
            return(g)
          }
        
        
          if(EstMethod!="ModFirth"){
            # for MLE and Firth, point and interval estimates are computed.
            if(res_Est$converged==0){MLEErr="NonConv"}else{MLEErr=""}
            sum1<-summary(res_Est)
            MLECoef <- as.vector(sum1$coefficients[,1])
            MLESE =as.vector(sum1$coefficients[,2])
            confcoef<- qnorm(1- 0.05/2)
            MLE_WaldCI = cbind(MLECoef[2:3]-MLESE[2:3]*confcoef, MLECoef[2:3]+MLESE[2:3]*confcoef)
            #maxloglkhd =sum1$deviance
           
            if(res_Est$converged==0){MLE_PLCI<-matrix(c(0,0,0,0),nrow =2) # for non-convergence
              B1LogLike=B2LogLike=0 
              maxloglkhd=Inf
            }else{
              # when GLM is converged 
              if(EstMethod=="MLE"){
                maxloglkhd=.likelihood3(MLECoef,data1)
                list.psi1=c(0,betaV[1])
                # B1LogLike stores values of log-likelihood for beta1 = 0 and actual value of beta1
                B1LogLike <- sapply(list.psi1, logPL, data=data1, thetainit=MLECoef, floglik=.likelihood3, 
                                    fscore=.gradient, indpsi=2, trace=FALSE) 
                list.psi2=c(0,betaV[2])
                # B2LogLike stores values of log-likelihood for beta2 = 0 and actual value of beta2
                B2LogLike <- sapply(list.psi2, logPL, data=data1, thetainit=MLECoef, floglik=.likelihood3, 
                                    fscore=.gradient, indpsi=3, trace=FALSE) 
              }else{
                maxloglkhd=.likelihood3Firth(MLECoef,data1)
                list.psi1=c(0,betaV[1])
                # B1LogLike stores values of log-likelihood for beta1 = 0 and actual value of beta1
                B1LogLike <- sapply(list.psi1, logPL, data=data1, thetainit=MLECoef, floglik=.likelihood3Firth, 
                                    fscore=.gradientFirth, indpsi=2, trace=FALSE) 
                list.psi2=c(0,betaV[2])
                # B2LogLike stores values of log-likelihood for beta2 = 0 and actual value of beta2
                B2LogLike <- sapply(list.psi2, logPL, data=data1, thetainit=MLECoef, floglik=.likelihood3Firth, 
                                    fscore=.gradientFirth, indpsi=3, trace=FALSE) 
              }
               Qchi2<- qchisq(0.95,1)/2 # for alpha = 0.05
               
               # Power computation for Beta1  PLCI
               if(B1LogLike[1] < (maxloglkhd-Qchi2) ){PowerPL1 = 1;}

               # Power computation for Beta2  PLCI
               if(B2LogLike[1] < (maxloglkhd-Qchi2) ){PowerPL2 = 1;}

               # CI coverage for Beta1  PLCI
               if(B1LogLike[2] < (maxloglkhd-Qchi2) ){
                 CoveragePL1 = 0;
                 if(betaV[1]>MLECoef[2]){LowCoveragePL1=1;UpCoveragePL1=0}else{LowCoveragePL1=0;UpCoveragePL1=1}
               }else{CoveragePL1 = 1;LowCoveragePL1=1;UpCoveragePL1=1}

               # CI coverage for Beta2  PLCI
                if(B2LogLike[2] < (maxloglkhd-Qchi2) ){
                  CoveragePL2 = 0;
                  if(betaV[2]>MLECoef[3]){LowCoveragePL2=1;UpCoveragePL2=0}else{LowCoveragePL2=0;UpCoveragePL2=1}
                }else{CoveragePL2 = 1;LowCoveragePL2=1;UpCoveragePL2=1}
            }
            
            # Power computation for Beta1  Wald CI
            if(MLE_WaldCI[1,1]<= 0 & MLE_WaldCI[1,2] >= 0){PowerWald1 = 0}else{PowerWald1 = 1}

            # Power computation for Beta2  Wald CI
            if(MLE_WaldCI[2,1]<= 0 & MLE_WaldCI[2,2] >= 0){PowerWald2 = 0}else{PowerWald2 = 1}

            # CI coverage for Beta1  Wald CI
            if(MLE_WaldCI[1,1]<= betaV[1] ){LowCoverageWald1 = 1}else{LowCoverageWald1 = 0}
            if(MLE_WaldCI[1,2]>= betaV[1] ){UpCoverageWald1 = 1}else{UpCoverageWald1 = 0}
            if(MLE_WaldCI[1,1]<= betaV[1] & MLE_WaldCI[1,2] >= betaV[1]){CoverageWald1 = 1}else{CoverageWald1 = 0}

            # CI coverage for Beta2  Wald CI
            if(MLE_WaldCI[2,1]<= betaV[2] ){LowCoverageWald2 = 1}else{LowCoverageWald2 = 0}
            if(MLE_WaldCI[2,2]>= betaV[2] ){UpCoverageWald2 = 1}else{UpCoverageWald2 = 0}
            if(MLE_WaldCI[2,1]<= betaV[2] & MLE_WaldCI[2,2] >= betaV[2]){CoverageWald2 = 1}else{CoverageWald2 = 0}
            
            simid<- substr(filenames[succ],(str_locate_all(filenames[succ], "_")[[1]][length(str_locate_all(filenames[succ], "_")[[1]])] +1),(str_locate_all(filenames[succ], ".c")[[1]][2] -2))
            CoefArr <- c(simid,file1,MLEErr,MLECoef,"",MLESE)
           
            CIVal<- c("",MLE_WaldCI[1,1],MLE_WaldCI[1,2],MLE_WaldCI[2,1],MLE_WaldCI[2,2],PowerWald1, CoverageWald1, LowCoverageWald1,UpCoverageWald1,PowerWald2, CoverageWald2, LowCoverageWald2,UpCoverageWald2,
                  "",PowerPL1, CoveragePL1, LowCoveragePL1,UpCoveragePL1,PowerPL2, CoveragePL2, LowCoveragePL2,UpCoveragePL2)
         }else{
            
            # For EstMethod = ModFirth only beta estimates and SE are computed for augmented data.
            if(res_Est$converged==0){MLEErr="NonConv"}else{MLEErr=""}
            sum1<-summary(res_Est)
            MLECoef <- as.vector(sum1$coefficients[,1])[-(ncov+2)]
            MLESE =as.vector(sum1$coefficients[,2])[-(ncov+2)]
    
            CoefArr <- c(simid,file1,MLEErr,MLECoef,"",MLESE)
            CIVal<- c()
          }
          return(t(c(CoefArr,CIVal)))
         
        }
        
        # Number of datasets for each combination 
  	    seqSims =1:length(filenames)
        
        # create a vector of variables to be sent in to the function that runs model
        varVec <- c("ncov","betaV","nobs","Nummethod","SimDataDir1","EstMethod","filenames")

        # create clusters
        cl <- makeCluster(7) # set number of parallel threads 

        clusterCall(cl, function(){
          # specify the functions/libraries required by each paralle thread
          library("MASS")
          library("stringr")
          library("data.table")
          library("stats")
          library("utils")
          library("likelihoodAsy")
          library("brglm2")
        })

        clusterExport(cl, varVec)
        z <- clusterApplyLB(cl = cl, x = seqSims, fun = SimCompute, filenames,EstMethod, Nummethod)

        stopCluster(cl)

        df <- data.frame(lapply(data.frame(t(sapply(z, `[`)), stringsAsFactors = FALSE), unlist))
        fwrite(df,paste(workdir,"Output/Coeff_",EstMethod,"_",Nummethod,"_",ncov,"_",round(betaV[1],2),"_",nobs,".csv",sep=""),sep=",",append = TRUE,col.names = FALSE, row.names=FALSE)
     
      }# for loop of Beta1 ends
    }# for loop of EPV ends
  }# for loop of ncov ends
  
#}
#EvalEst(workdir,SimDataDir, betaV,EstMethod="MLE",CovList=(2),EPVList=c(3), Beta1List=c(1,2) )
#debug(EvalEst)
