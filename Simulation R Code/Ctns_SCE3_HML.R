
vi S3PP_ctns_wt_CV.R


setwd("~/trim_05")
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
#args accepts arguments for simulation settings
print(args)

wt=args[1]

V_lmd2=seq(as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]))
V_lmd1=seq(as.numeric(args[5]),as.numeric(args[6]),as.numeric(args[7]))

T_rep=as.numeric(args[8])
b_num1=as.numeric(args[9])
b_num2=as.numeric(args[10])
b_tag1=args[11]



set.seed(142)
library(Matrix)
library(glmnet)
library(reshape2)
library(dplyr)
library(gdata)
library(foreach)
library(doParallel)

cl <- 5 #not to overload Your computer
registerDoParallel(cl)

#y_hat=X*beta_est=c(1,X)*(b0,b)
#include intercept in X
Linear_Pred<-function(X, beta){
    beta<-data.matrix(beta)
    linear_pred=data.matrix(X)%*%beta
        return(linear_pred)
}

#solution path for a lambda2 and a series of lambda1 
Cnt_Spth<-function(X, Y, A_matrix, lmd2, V_lmd1){
      
    #save solution path for all lambda1
    beta_lmd1<-matrix(nrow=p+1)[,-1]
    
    #design matrix X*
    sqrt_lmd2<-sqrt(2^(lmd2))
      
    Xp1<-cbind(1,X)
    Ap0<-cbind(0,sqrt_lmd2*A_matrix)
    colnames(Xp1)<-colnames(Ap0)
    
    X_aug<-rbind(Xp1,Ap0)

    #Y*
    Y_aug<-append(Y,rep.int(0,dim(A_matrix)[1]))

    #Note that any variable with penalty.factor equal to zero is not penalized at all
    #all penalty factor set to 1 (lenght equal to # columns of Xaugmented)
    p.fac <- rep(1, dim(X_aug)[2])
    #set the first penalty factor as 0
    p.fac[1] <- 0
    #exclude ***intercept(the first column of Xaugmented) from lasso penalty  
    fit_tr <- glmnet(X_aug,Y_aug, alpha = 1, nlambda = 100, intercept=FALSE, penalty.factor=p.fac)

    #intercept=FALSE gives a coef with 0 at the first row, exclude it in the reported values 
    beta_lmd1<-coef(fit_tr, s = 2^V_lmd1, exact = FALSE)[-1,]

    return(beta_lmd1)
} 

#common CNV, freq>0.05
#read in files 
# output files
in_label=c("05")
out_label=c(paste0("~/ctns",b_tag1),b_tag1)

#data CNV, covariants, and Y  
Xmrg<-read.table(paste0("Xmrg_",in_label,"parmCHR10.txt"),header=TRUE)
dim(Xmrg)

A_matrix=read.table(paste0(wt,"_parmCHR10.txt"),header=TRUE)

N=dim(Xmrg)[1]
p=dim(Xmrg)[2]-2 # do not include ID and Y




X<-data.matrix(Xmrg[,2:(p+1)]) # all CNV and covariants as X 

#causal segments
CNV_p1<-c(145:147,225:227)
p1=length(CNV_p1)
CNV_p2<-c(150:152,231:233)
p2=length(CNV_p2)




# True effect size
b_0<-(-2)



beta_True<-matrix(0,nrow=(p+1),ncol=1)
beta_True[1,1]=b_0
beta_True[1+CNV_p1,1]=b_1
beta_True[1+CNV_p2,1]=b_2



#simulated Y_bar 
xb=Linear_Pred(cbind(1,X),beta_True)
print(sd(xb)^2)


print(length(xb[xb!=b_0]))
print("length of xb==b_0")
print(length(xb[which(xb==b_0)]))
print("length of xb")
print(length(xb))

#save results for final output
beta_cv<-matrix(0,ncol=1,nrow=p+1)[,-1]
beta_cv<-data.frame(beta_cv)



for(rep in c(1:T_rep)){

   #simulate Y
    Xmrg$Y<-data.matrix(rnorm(length(xb), mean=xb, sd=1))
    print("summary of Y")
    print(summary(Xmrg$Y))
    print("summary of Y: xb>b_0 (with causal segments)")
    print(summary(Xmrg$Y[xb!=b_0]))  
    print("summary of Y: xb<=b_0 (no causal segments)")
    print(summary(Xmrg$Y[xb==b_0]))
    print("length of Y")
    print(length(Xmrg$Y))
   
    #generate uniform data, but quantile considers xb<=b_0 and xb>b_0 seperately
    #ensure samples with causal segments are randomly distributed to all folds
    Xmrg$tr<-runif(dim(Xmrg)[1],0,1)
    Xmrg$tr<- case_when(
         Xmrg$tr<=quantile(Xmrg$tr[xb==b_0],0.2)~1,
        (Xmrg$tr<=quantile(Xmrg$tr[xb==b_0],0.4)&Xmrg$tr>quantile(Xmrg$tr[xb==b_0],0.2))~2,
        (Xmrg$tr<=quantile(Xmrg$tr[xb==b_0],0.6)&Xmrg$tr>quantile(Xmrg$tr[xb==b_0],0.4))~3,
        (Xmrg$tr<=quantile(Xmrg$tr[xb==b_0],0.8)&Xmrg$tr>quantile(Xmrg$tr[xb==b_0],0.6))~4,
         Xmrg$tr>quantile(Xmrg$tr[xb==b_0],0.8)~5
        )
    #generate uniform for samples with xb!=b_0, and find quantile 
    Xmrg$tr[xb!=b_0]<-runif(length(xb[xb!=b_0]),0,1)
    Xmrg$tr<- case_when(
        Xmrg$tr==1~1,
        Xmrg$tr==2~2,
        Xmrg$tr==3~3,
        Xmrg$tr==4~4,
        Xmrg$tr==5~5,
        Xmrg$tr<=quantile(Xmrg$tr[xb!=b_0],0.2)~1,
        (Xmrg$tr<=quantile(Xmrg$tr[xb!=b_0],0.4)&Xmrg$tr>quantile(Xmrg$tr[xb!=b_0],0.2))~2,
        (Xmrg$tr<=quantile(Xmrg$tr[xb!=b_0],0.6)&Xmrg$tr>quantile(Xmrg$tr[xb!=b_0],0.4))~3,
        (Xmrg$tr<=quantile(Xmrg$tr[xb!=b_0],0.8)&Xmrg$tr>quantile(Xmrg$tr[xb!=b_0],0.6))~4,
        Xmrg$tr<1&Xmrg$tr>quantile(Xmrg$tr[xb!=b_0],0.8)~5
        )
  
    print("table of training set for all samples")
    print(table(Xmrg$tr))
    print("table of training set for xb!=b_0")
    print(table(Xmrg$tr[xb!=b_0]))
    print("table of training set for xb==b_0")
    print(table(Xmrg$tr[xb==b_0]))
   # print("table of training set for xb>0")
   # print(table(Xmrg$tr[xb>0]))

    #prepare data for regression
    Y=Xmrg$Y
    colnames(Y)<-"Y"

    
    ###Model fitting#############
    #save RSS loss for each fold, lambda1, and lambda2 
    
    #five fold cross-validation, the same fold as in the lasso_cv method
    loss_matrix<-foreach(i=c(1:5), .combine=rbind)%dopar%{
        loss_fold<-matrix(nrow=1)[-1,]
        #set training data 
        X_tr<-X[Xmrg["tr"]!=i,]
        Y_tr<-data.matrix(Y[Xmrg["tr"]!=i,])
          
        #predict in x_val and Y_val data.
        X_val<-X[Xmrg["tr"]==i,]
        Y_val<-data.matrix(Y[Xmrg["tr"]==i,])

        Xval_N=dim(X_val)[1]
        for (lmd2 in V_lmd2){

            #beta coef for one lmd2
            
            beta_lmd1<-Cnt_Spth(X_tr,Y_tr,A_matrix,lmd2,V_lmd1)
            
            #predict Y_hat for RSS calculation in validation set 

            Y_hat_val<-Linear_Pred(cbind(1,X_val),beta_lmd1)

            #loss info : fold, lmd2, lmd1, 1/n sum(Y_hat-Y_val)
            loss_fold<-rbind(loss_fold,cbind(i,lmd2,V_lmd1,colSums((Y_hat_val-matrix(Y_val, nrow=dim(Y_hat_val)[1],ncol=dim(Y_hat_val)[2]))^2)/length(Y_val)))
        }
    loss_fold
    }

    loss_matrix<-data.frame(loss_matrix)
    colnames(loss_matrix)<-c("fold","lmd2","lmd1","rss")
    #average rss over 5 fold
    loss_sum<-data.frame(summarise(loss_matrix%>%group_by(lmd2,lmd1),
                rss_avg=mean(rss),.groups = 'drop'))

    #print(dim(loss_sum))


    idx<-which(loss_sum$rss_avg==min(loss_sum$rss_avg))

    wide_rd<-dcast(loss_sum,lmd2~lmd1,value.var="rss_avg", fun.aggregate = median)
    print(wide_rd)
      
    #minimum loss
    Mloss=loss_sum[idx,]
  

    #tunning parameters and beta coefficients
    lmd2=as.numeric(Mloss[1,1])
    lmd1=as.numeric(Mloss[1,2])
    print(c("rep","lambda2","lambda1"))
 
    print(c(rep,lmd2,lmd1))
    beta_Y<-Cnt_Spth(X,Y,A_matrix,lmd2,lmd1)
    

    beta_cv<-cbind(beta_cv,beta_Y)  

    if(rep%%50==0){
        write.table(beta_cv,file = paste0(out_label[1],wt,out_label[2],".txt"), sep = "\t", row.names = FALSE)
    }  
    
}
    

T_rep=100


for (b_tag in c("4P4", "4P2", "4P1")){
    if(b_tag=="4P4"){
       
        b_2=0.4
    }else if(b_tag=="4P2"){
        b_2=0.2
    }else{
        b_2=0.1
    }  
 for (wt in c("pw_1","pw_cs")){      
    if(!file.exists(paste0("ctns", wt, b_tag, ".bash"))){
        file.create(paste0("ctns",wt,b_tag, ".bash"))
    }
          
        fileConn<-file(paste0("ctns",wt,b_tag, ".bash"))

        writeLines(c("#!/bin/bash",
                     paste0("#SBATCH --output=S3PP_ctns", wt, b_tag, ".out"),
                    "hostname",
                    "pwd",
        #fit weighted smoothness
        paste0("/opt/R/4.1.3/bin/Rscript S3PP_ctns_wt.R ",wt," ",V2_1," ", V2_2, " ", V2_p," ", V1_1, " ",V1_2, " ", V1_p, " ", T_rep," 0.4 ",b_2," ",b_tag)
        ), fileConn)

    close(fileConn)
 }

}

for (b_tag in c("4P4", "4P2", "4P1")){
    for (wt in c("pw_1","pw_cs")){

     print(paste0("sbatch ctns",wt, b_tag,".bash"))
    }
}

