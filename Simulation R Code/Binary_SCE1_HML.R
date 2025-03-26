## code used on HPC to simulate a binary scenario 1 with causal CNV
## in the middle of a consective region.
vi bnry_wt_S1TX.R

#Check correlation: 

# data file: XY and weight 
setwd("~/trim_05")
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
#
print(args)

wt=args[1]

V_lmd2=seq(as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]))
V_lmd1=seq(as.numeric(args[5]),as.numeric(args[6]),as.numeric(args[7]))

T_rep=as.numeric(args[8])
###set simulation signal level with b_num1
b_num1=as.numeric(args[9])

b_tag1=args[10]

#packages
set.seed(142)
library(Matrix)
library(glmnet)
library(reshape2)
library(dplyr)
library(gdata)
library(foreach)
library(doParallel)

cl <- 5 #not to overload your computer
registerDoParallel(cl)

#y_hat=X*beta_est=c(1,X)*(b0,b)
#include intercept in X
Linear_Pred<-function(X, beta){
    beta<-data.matrix(beta)
    linear_pred=data.matrix(X)%*%beta
        return(linear_pred)
}

#calculate probability based on design*beta
Prob_Pred<-function(X){
    p=1/(exp(-X)+1)
    return(p)
}

#original** log likelihood(not negative )
LogLH<-function(X,beta,Y){
    linear_pred<-Linear_Pred(X,beta)
    ll<-colSums(linear_pred*as.vector(Y)-log(1+exp(linear_pred)))
    return(ll)
}


Rwls_Spth<-function(X, Y, A_matrix, lmd2, V_lmd1){
   # print(lmd2)
   #save beta coefficient for a series of lambda 1 (for a given lambda 2)
    beta_lmd1<-matrix(nrow=p+1)[,-1]
    #iteration info
    lmd_sum21<-matrix(ncol=5)[-1,]
    colnames(lmd_sum21)<-c("lmd2","lmd1","iter","loss_dif","beta_dif")

    ##initial beta values, default type.measure="deviance"   
    XN=dim(X)[1]
    Xp=dim(X)[2]
    print(c(XN,Xp))
    AN=dim(A_matrix)[1]

    fit_yb_init<- glmnet(data.matrix(X), Y, family="binomial", alpha = 1, standardize = TRUE)
    #initial coef at candidate lambda 1 values
    beta_0<-data.matrix(coef(fit_yb_init,s=2^V_lmd1))

    #iteratively update each initial coef
    for(N_lmd1 in c(1:length(V_lmd1))){
        #iteration 0_initial
        iter=0  
        lmd1=V_lmd1[N_lmd1]
        
       # beta_matrix<-matrix(nrow=p+1)[,-1]
        beta_cur<-data.matrix(beta_0[,N_lmd1])

        #stop criteria
        #   1.coef change: set initial change of beta values as max(abs(beta))
        maxdif<-max(abs(beta_cur))
        #beta_matrix<-cbind(beta_matrix,beta_cur)
       
        #2. loss change: Negtive log likelihood,set initial as current loss 
        loss<-(-1/XN)*LogLH(cbind(1,X),beta_cur,Y)
        loss_dif<-loss

        #print(c(iter,loss,maxdif))

        #for data augmentation
        sqrt_term<-sqrt(2*(XN+AN)*(2^(lmd2)))
        while((iter<8)&&(maxdif>10^(-3))&&(loss_dif>10^(-6))){
            ts_iter<-Sys.time()

            #prepare linear prediction and probability prediction for iteration update
            linear_pred<-Linear_Pred(cbind(1,X),beta_cur)
            prob_pred<-Prob_Pred(linear_pred)

#the predicted value are all equal no need to update
            if(var(linear_pred)<10^(-4)){
                break
            }

            #update u*, v* 
            iter=iter+1
            u=data.matrix(linear_pred+(Y-prob_pred)/(prob_pred*(1-prob_pred)))

            #sqrt_v<-diag(as.vector(sqrt(prob_pred*(1-prob_pred))),nrow=length(u),ncol=length(u))
            sqrt_v<-as.matrix(sqrt(prob_pred*(1-prob_pred)),ncol=1)
            #sqrt_v[1:10]
            #update X_augmented               
            Xp1=cbind(1,X)
            colnames(Xp1)<-c(paste0("X",c(0:(dim(Xp1)[2]-1))))
            Ap0<-cbind(0,A_matrix)
            colnames(Ap0)<-c(paste0("X",c(0:(dim(Ap0)[2]-1))))
            X_aug<-rbind(Xp1*as.vector(sqrt_v),sqrt_term*Ap0)
            #X_aug<-rbind(sqrt_v%*%Xp1,sqrt_term*Ap0)
              

            ##update Y_augmented 
            y_app<-data.matrix(rep.int(0,dim(A_matrix)[1]))
            colnames(y_app)<-"y"
            #y_org<-data.matrix(sqrt_v%*%u)
            y_org<-data.matrix(u*sqrt_v)
            colnames(y_org)<-"y"
            Y_aug<-rbind(y_org,y_app)
           
            #betaprint(dim(X_aug))
            #betaprint(length(Y_aug))
            #betaX_aug<-X_aug[which(!is.na(Y_aug)),]
            #betaY_aug<-Y_aug[which(!is.na(Y_aug))]
            #betaprint(dim(X_aug))
            #betaprint(length(Y_aug))
                                  
            #Note that any variable with penalty.factor equal to zero is not penalized at all
            #all penalty factor set to 1 (lenght equal to # columns of Xaugmented)
            p.fac <- rep(1, dim(X_aug)[2])
            #set the first penalty factor as 0
            p.fac[1] <- 0  
            #exclude ***intercept(the first column of Xaugmented) from lasso penalty                  
            fit_tr<-glmnet(data.matrix(X_aug),data.matrix(Y_aug), alpha = 1, nlambda = 100, intercept=FALSE, penalty.factor=p.fac)
            
            #update beta coefficient
            #intercept=FALSE gives a coef with 0 at the first row, exclude it in the reported values 
            beta_next<-data.matrix(coef(fit_tr,s=2^(lmd1), exact = FALSE)[-1,])                
            
            #stop criteria  beta               
            beta_dif<-abs(beta_next-beta_cur)
            maxdif<-max(abs(beta_dif))
                  
            #stop criteria loss function 
            loss_next<-(-1/XN)*LogLH(cbind(1,X),beta_next,Y)
            loss_dif<-loss-loss_next
            
            if(is.infinite(loss)){
                loss_dif=0
                iter=iter-1
                break
            }
             
            if((loss_dif<0)){
                
                loss_dif=abs(loss_dif)
                #iter=iter-1
                #break
            }

            loss=loss_next
            #prepare for next iteration
            beta_cur<-data.matrix(beta_next)
           # beta_matrix<-cbind(beta_matrix,beta_cur)

            #print(c(iter,loss,maxdif))
            #print(beta_matrix)

            te_iter<-Sys.time()
            t_iter<-te_iter-ts_iter
        }
        #beta_matrix

        beta_lmd1<-cbind(beta_lmd1,beta_cur)
        lmd_sum21<-rbind(lmd_sum21,c(lmd2,lmd1,iter,loss_dif,maxdif))
        print(c(lmd2,lmd1,iter,loss_dif,maxdif))
    }

    return(list(beta_lmd1,lmd_sum21))
}

#common CNV, freq>0.05
#read in files 
# output files
in_label=c("05")
out_label=c(paste0("~/out",b_tag1,"/bnr",b_tag1),b_tag1)

#data CNV, covariants, and Y  
Xmrg<-read.table(paste0("Xmrg_",in_label,"parmCHR10.txt"),header=TRUE)
dim(Xmrg)

A_matrix=read.table(paste0(wt,"_parmCHR10.txt"),header=TRUE)

N=dim(Xmrg)[1]
p=dim(Xmrg)[2]-2 # do not include ID and Y


X<-data.matrix(Xmrg[,2:(p+1)]) 

#causal segments
CNV_p1<-c(145:152,225:233)
p1=length(CNV_p1)



# True effect size
b_0<-(-2)
b_1=b_num1

beta_True<-matrix(0,nrow=(p+1),ncol=1)
beta_True[1,1]=b_0
beta_True[1+CNV_p1,1]=b_1





#simulated linear predictor
xb=Linear_Pred(cbind(1,X),beta_True)
print(sd(xb)^2)


p_b=Prob_Pred(xb)
print(summary(p_b))


print(length(xb[xb!=b_0]))
print("length of xb==b_0")
print(length(xb[which(xb==b_0)]))
print("length of xb")
print(length(xb))

#save results for final output
beta_cv<-matrix(0,ncol=1,nrow=p+1)[,-1]
beta_cv<-data.frame(beta_cv)


#######################
#simulate Y repeats
for(rep in c(1:T_rep)){

   #generate Y once,split the data
    Xmrg$Y<-data.matrix(rbinom(n=length(p_b)[1], size=1, prob=p_b))
    print("table of Y")
    print(table(Xmrg$Y))
    print("table of Y: xb>0 (with causal segments)")
    print(table(Xmrg$Y[xb!=b_0])) 
    print(table(Xmrg$Y[xb==b_0])) 

#generate uniform data, but quantile considers xb==b_0 and xb!=b_0 seperately
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
   # print("table of training set for xb<=0")
   # print(table(Xmrg$tr[xb<=0]))

    #prepare data for regression
    Y=Xmrg$Y
    colnames(Y)<-"Y"


    ###Model fitting#############
    #save -LogL loss for each fold, lambda1, and lambda2
    
    loss_matrix<-foreach(i=c( 1:5), .combine=rbind)%dopar%{

        loss_fold<-matrix(nrow=1)[-1,]
        #set training data 
        X_tr<-X[Xmrg["tr"]!=i,]
        Y_tr<-data.matrix(Y[Xmrg["tr"]!=i,])
          
        #predict in x_test and y_test data.
        X_val<-X[Xmrg["tr"]==i,]
        Y_val<-data.matrix(Y[Xmrg["tr"]==i,])

        Xval_N=dim(X_val)[1]
        for (lmd2 in V_lmd2){
            
            lmd2_spth<-Rwls_Spth(X_tr,Y_tr,A_matrix,lmd2,V_lmd1)
            beta_lmd1<-lmd2_spth[[1]]
            lmd_sum21<-lmd2_spth[[2]]

            loss<-(-1/Xval_N)*LogLH(cbind(1,X_val),beta_lmd1,Y_val)
            loss_fold<-rbind(loss_fold,cbind(i,lmd_sum21[,1:3],loss))
        }

        loss_fold
    }
#average loss of 5 folds
    loss_matrix<-data.frame(loss_matrix)
    colnames(loss_matrix)<-c("fold","lmd2","lmd1","iter","loss")

    loss_sum<-data.frame(summarise(loss_matrix%>%group_by(lmd2,lmd1),
                loss_avg=mean(loss),.groups = 'drop'))

    #print(dim(loss_sum))

    #best lmd1+lmd2 and coefficients
    idx<-which(loss_sum$loss==min(loss_sum$loss))

    wide_rd<-dcast(loss_sum,lmd2~lmd1,value.var="loss_avg", fun.aggregate = median)
    print(wide_rd)
    
    #minimum loss
    Mloss=loss_sum[idx,]
  
    #tunning parameters and beta coefficients
    lmd2=as.numeric(Mloss[1,1])
    lmd1=as.numeric(Mloss[1,2])
    print(c("rep","lambda2","lambda1"))

    print(c(rep,lmd2,lmd1))

    rst_reg<-Rwls_Spth(X,Y,A_matrix,lmd2,lmd1)
    beta_y=data.matrix(rst_reg[[1]])

    beta_cv<-cbind(beta_cv,beta_y)  

    if(rep%%20==0){
    write.table(beta_cv,file = paste0(out_label[1],wt,out_label[2],".txt"), sep = "\t", row.names = FALSE)
    }  
}

