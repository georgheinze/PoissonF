### This R code nested loops from output of simulation study performed for Modified Firth method for Poisson regression
### Specify working directory where the output files are kelp as 'workdir'. All the plots will be saved in this directory.
### Author: Ashwini Joshi
### Date: 1-Mar-2020
####################################################################################################################

library(ggplot2)
library(loopR)

workdir= "<Path>"  # Specify working directory path here.
setwd(workdir)

#################################################################################

# Read the csv file of prediction RMSE
RMSE_Prediction <- read.csv("RMSE_Prediction.csv",sep=",",header=TRUE)

RMSE_Prediction<-RMSE_Prediction[,-5] # remove column of MLE-Firth method

design1<-RMSE_Prediction[,1:6]

design1$Beta1<-round(design1$Beta1,2)
# Call nested loop function
nested_loop(design1, x="Beta1", steps="nobs", grid_cols="Num.ofCov",steps_base = min(design1[,4:6])-70, steps_height=min(design1[,4:6])/4, colors=c("black", "blue", "red"), steps_names_annotate=FALSE,  steps_values_annotate=TRUE,line_size = 0.5, y_name=expression(RMSE(mu)), x_name=expression(beta[1]), design="fractional", 
            legend_labels = c("MLE","Firth","ModFirth"),theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0,linetype=3)

ggsave("prdiction.png",width = 7.4, height = 4.2 ) # save the plot
#################################################################################

RMSE_Beta1 <- read.csv("RMSE_Beta1.csv",sep=",",header=TRUE) # Read the csv file of Beta1 RMSE

design2<- RMSE_Beta1[,1:6]
design2<-design2[,-5] # remove column of MLE-Firth method
design2$Beta1<-round(design2$Beta1,2)
nested_loop(design2, x="Beta1", steps="nobs", grid_cols="Num.ofCov", steps_base = min(design2[,4:5]),steps_height=max(min(design2[,4:5]/4,10),10), colors=c("black", "blue"),steps_names_annotate=FALSE, steps_values_annotate=TRUE,
            y_name=expression(RMSE(beta[1])), x_name=expression(beta[1]), design="fractional", legend_labels = c("MLE","Firth"),theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0,linetype=3)
 
ggsave("RMSE.png",width = 7.4, height = 4.2 )# save the plot
 
design2Exact<-RMSE_Beta1[which(RMSE_Beta1[,1]< 6 & RMSE_Beta1[,3]<500),1:7]
design2Exact <-design2Exact[,-5] # remove column of MLE-Firth method
design2Exact$Beta1<-round(design2Exact$Beta1,2)
nested_loop(design2Exact, x="Beta1", steps="nobs", grid_cols="Num.ofCov", steps_base = min(design2Exact[,4:6])-10,steps_height=10, colors=c("black", "blue", "red"), steps_names_annotate=FALSE, steps_values_annotate=TRUE,
             y_name=expression(RMSE(beta[1])), x_name=expression(beta[1]), design="fractional", legend_labels = c("MLE","Firth","Exact"),theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0,linetype=3)
 
ggsave("RMSE_Ex.png",width = 7.4, height = 4.2 ) # save the plot

#################################################################################

Bias_Beta1 <- read.csv("Bias_Beta1.csv",sep=",",header=TRUE) # Read the csv file of bias in Beta1 estimation

design3<- Bias_Beta1[,1:6]
design3<-design3[,-5] # remove column of MLE-Firth method
design3$Beta1<- round(design3$Beta1,2)
nested_loop(design3, x="Beta1", steps="nobs", grid_cols="Num.ofCov", steps_base = min(design3[,4:5]),steps_height=10, colors=c("black", "blue"),steps_names_annotate=FALSE, steps_values_annotate=TRUE,
            y_name=expression(Bias(Beta[1])), x_name=expression(beta[1]), design="fractional", legend_labels = c("MLE","Firth"),theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0,linetype=3)
ggsave("Bias.png",width = 7.4, height = 4.2 ) # save the plot

design3Exact<- Bias_Beta1[which(Bias_Beta1[,1]< 6 & Bias_Beta1[,3]<500),1:7]
design3Exact<- design3Exact[,-5] # remove column of MLE-Firth method
design3Exact$Beta1<-round(design3Exact$Beta1,2)
nested_loop(design3Exact, x="Beta1", steps="nobs", grid_cols="Num.ofCov", steps_base = min(design3Exact[,4:6])-10 ,steps_height=10, colors=c("black", "blue", "red"), steps_names_annotate=FALSE, steps_values_annotate=TRUE,
            y_name=expression(Bias(Beta[1])), x_name=expression(beta[1]), design="fractional", legend_labels = c("MLE","Firth","Exact"),theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0,linetype=3)
ggsave("Bias_Ex.png",width = 7.4, height = 4.2 ) # save the plot

#################################################################################

PowerCalcn<- read.csv("PowerCalcn.csv",sep=",",header=TRUE) # Read the csv file of Power of CI


design<-PowerCalcn[,1:6]
design4<- design[,-5] # remove column of MLE-Firth method
design4$Beta1<-round(design4$Beta1,2)
design4$MLE<-round(design4$MLE/1000,3)
design4$Firth<-round(design4$Firth/1000,3)
nested_loop(design4, x="Beta1", steps="nobs", grid_cols="Num.ofCov", steps_base = min(design4[,4:5]),steps_height=0.01, colors=c("black", "blue"), steps_names_annotate=FALSE,steps_values_annotate=TRUE,
            y_name="Power", x_name=expression(beta[1]), design="fractional",legend_labels = c("MLE","Firth"), theme_add=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0,linetype=3)
ggsave("Power.png",width = 7.4, height = 4.2 ) # save the plot

design4Exact<-PowerCalcn[which(PowerCalcn[,1]< 6 & PowerCalcn[,3]<500),1:8]
design4Exact <- design4Exact[,-5] # remove column of MLE-Firth method
design4Exact$Beta1<-round(design4Exact$Beta1,2)
design4Exact$MLE<-round(design4Exact$MLE/1000,3)
design4Exact$Firth<-round(design4Exact$Firth/1000,3)
design4Exact$CCMLE..MUE.<-round(design4Exact$CCMLE..MUE./1000,3)
design4Exact$Mid.p<-round(design4Exact$Mid.p/1000,3)
nested_loop(design4Exact, x="Beta1", steps="nobs", grid_cols="Num.ofCov", steps_base = min(design4Exact[,4:7]),steps_height=0.01, colors=c("black", "blue", "red", "green"),steps_names_annotate=FALSE, steps_values_annotate=TRUE,
            y_name="Power", x_name=expression(beta[1]), design="fractional",legend_labels = c("MLE","Firth","Exact","Mid-p"), theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0,linetype=3)
ggsave("Power_Ex.png",width = 7.4, height = 4.2 ) # save the plot

#################################################################################

CICalcn<- read.csv("CIValues.csv",sep=",",header=TRUE) # Read the csv file of coverage of CI

design<-CICalcn[,1:3]
design5<-CICalcn[,1:6]
design5<-design5[,-5] # remove column of MLE-Firth method
design5$Beta1<-round(design5$Beta1,2)
design5$MLE_Wald_Low<-round(design5$MLE_Wald_Low/1000,3)
design5$Firth_PPL_Low<-round(design5$Firth_PPL_Low/1000,3)
nested_loop(design5, x="Beta1", steps="nobs", grid_cols="Num.ofCov", steps_base = min(design5[,4:5],0.8),steps_height= 0.01, colors=c("black", "blue"),steps_names_annotate=FALSE, steps_values_annotate=TRUE,
            y_name="LowerCICoverage", x_name=expression(beta[1]), design="fractional",legend_labels = c("MLE","Firth"), theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0.975,linetype=3)
ggsave("LowCI.png",width = 7.4, height = 4.2 ) # save the plot

design6<-cbind(design,CICalcn[,9:11])
design6<-design6[,-5] # remove column of MLE-Firth method
design6$Beta1<-round(design6$Beta1,2)
design6$MLE_Wald_Upper<-round(design6$MLE_Wald_Upper/1000,3)
design6$Firth_PPL_Upper<-round(design6$Firth_PPL_Upper/1000,3)
nested_loop(design6, x="Beta1", steps="nobs", grid_cols="Num.ofCov", steps_base = min(design6[,4:5],0.8) ,steps_height=0.01, colors=c("black", "blue"),steps_names_annotate=FALSE, steps_values_annotate=TRUE,
            y_name="UpperCICoverage", x_name=expression(beta[1]), design="fractional",legend_labels = c("MLE","Firth"), theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0.975,linetype=3)

ggsave("UppCI.png",width = 7.4, height = 4.2 ) # save the plot


design5Exact<-CICalcn[which(CICalcn[,1]< 6 & CICalcn[,3]<500),1:8]
design5Exact <- design5Exact[,-5] # remove column of MLE-Firth method
design5Exact$Beta1<-round(design5Exact$Beta1,2)
design5Exact$MLE_Wald_Low <-round(design5Exact$MLE_Wald_Low/1000,3)
design5Exact$Firth_PPL_Low<-round(design5Exact$Firth_PPL_Low/1000,3)
design5Exact$CMLE_MUE_Low  <-round(design5Exact$CMLE_MUE_Low/1000,3)
design5Exact$MidP_Low<-round(design5Exact$MidP_Low/1000,3)
nested_loop(design5Exact, x="Beta1", steps="nobs", grid_cols="Num.ofCov", steps_base = min(design5Exact[,4:7],0.8),steps_height=0.01, colors=c("black", "blue", "red", "green"),steps_names_annotate=FALSE, steps_values_annotate=TRUE,
            y_name="LowerCICoverage", x_name=expression(beta[1]), design="fractional", legend_labels = c("MLE","Firth","Exact","Mid-p"),theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0.975,linetype=3)
ggsave("LowCI_Ex.png",width = 7.4, height = 4.2 ) # save the plot

design6 <- cbind(design,CICalcn[,9:13])
design6Exact<-design6[which(design6[,1]< 6 & design6[,3]<500),1:8]
design6Exact <- design6Exact[,-5] # remove column of MLE-Firth method
design6Exact$Beta1<-round(design6Exact$Beta1,2)
design6Exact$MLE_Wald_Upper<-round(design6Exact$MLE_Wald_Upper/1000,3)
design6Exact$Firth_PPL_Upper<-round(design6Exact$Firth_PPL_Upper/1000,3)
design6Exact$CMLE_MUE_Upper<-round(design6Exact$CMLE_MUE_Upper/1000,3)
design6Exact$MidP_Upper <-round(design6Exact$MidP_Upper/1000,3)

nested_loop(design6Exact, x="Beta1", steps="nobs", grid_cols="Num.ofCov", steps_base = min(design6Exact[,4:7],0.8),steps_height=0.01, colors=c("black", "blue", "red", "green"),steps_names_annotate=FALSE, steps_values_annotate=TRUE,
            y_name="UpperCICoverage", x_name=expression(beta[1]), design="fractional", legend_labels = c("MLE","Firth","Exact","Mid-p"),theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "top")) + geom_hline(yintercept = 0.975,linetype=3)

ggsave("UppCI_Ex.png",width = 7.4, height = 4.2 ) # save the plot
######################################################################################################

SepData <- read.csv("SeparationData.csv",sep=",",header=TRUE) # Read the csv file of number of separations

design1<- SepData
design1$Beta1<-round(design1$Beta1,2)
design1$Sep<-round(design1$Sep / 100,2)
nested_loop(design1, x="Beta1", steps="N", grid_cols="Num.ofCov",steps_base = -10, steps_height=2, colors=c("black"),
            steps_names_annotate=FALSE,  steps_values_annotate=TRUE,line_size = 0.5, y_name="Separation %", x_name=expression(beta[1]), design="fractional", theme_add=theme(axis.text.x = element_text(angle = 270, hjust = 1),legend.position = "none")) + geom_hline(yintercept = 0,linetype=3)
ggsave("Separation.png",width = 7.4, height = 4.2 ) # save the plot

###########################################################################################################