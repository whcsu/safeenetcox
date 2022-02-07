rm(list = ls())
library(MASS)
library(glmnet)
library(survival)
library(reshape2)
library(ggplot2)
library(survcomp)
library(microbenchmark)

getdata<-function(n,p,rho){
  sig<-matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      sig[i,j]<-rho^(abs(i-j))
    }
  }
  X <- mvrnorm(n,rep(0,p),sig)
  s<-runif(n,0,1)
  beta<-c(5,5,5,-15*rho,rep(0,(p-4)))  ###真实beta1~beta4
  muu<-10*exp(X%*%beta)
  t<-(-1)*log(s)/muu
  c<-rexp(n,1/10)
  status <- as.numeric(t<=c)
  time <- pmin(t,c)   ###右删�?
  return(list(X=X,time=time,status=status))
}


###
set.seed(5)
data<-getdata(200,7000,0.50)
data<-data.frame(data$time,data$status,data$X)
data<-data[order(data$data.time),]
X<-as.matrix(data[,3:length(data)])
Y<-data[,1:2]
colnames(Y)[1]<-'time'
colnames(Y)[2]<-'status'
n=dim(X)[1]
p=dim(X)[2]
data<-data.frame(Y,X)

#################
encox_safer<-function(X,Y){
  n=dim(X)[1]
  p=dim(X)[2]
  c<-numeric(p)
  smax<-0
  smin<-0
  for(i in 1:(n-1)){
    if(Y[i,2]==1){
      c<-c+X[i,]
      smax<-smax+apply(X[i:n,],2,max)
      smin<-smin+apply(X[i:n,],2,min)
    }
  }
  if(Y[n,2]==1){
    c<-c+X[i,]
    smax<-smax+X[n,]
    smin<-smin+X[n,]
  }
  r1<-c-smin
  r2<-smax-c
  r<-pmax(r1,r2)
  return(r)
}


r<-encox_safer(X,Y)

dloglik = function(X,Y,beta){
  p = length(beta)
  L = numeric(p)
  for(i in 1:n){
    temp = 0.0
    temp1 = numeric(p)
    for(j in i:n){
      temp = temp + exp(sum(beta*X[j,]))
      temp1 = temp1 + X[j,]*exp(sum(beta*X[j,]))
    }
    L = L - Y[i,2]*(X[i,] - temp1/temp)
  }
  return(L) 
}

lammax<-function(X,Y,k,g){
  p=dim(X)[2]
  beta0<-rep(0,p)
  df0<-dloglik(X,Y,beta0)
  tmp<-sqrt(t(df0)%*%df0)
  lambda_max<-max(tmp)
  lambda_min<-g*lambda_max
  lambda<-numeric(100)
  for(i in 1:k){
    lambda[i]<-lambda_max*(lambda_min/lambda_max)^(i/k)
  }
  return(list(lambda=lambda,lambda_min=lambda_min,lambda_max=lambda_max))
}

safe_rrsr<-function(X,Y,r,lambda,alpha){
  n=dim(X)[1]
  p=dim(X)[2]
  rej_radio<-numeric(100)
  scr_radio<-numeric(100)
  for(k in 1:100){
    print(k)
    fit_elasn<-cv.glmnet(X,Surv(Y[,1],Y[,2]),family = "cox",alpha =alpha)
    elasn<-coef(fit_elasn,s=lambda[k]*0.001)
    elasn_logi<-elasn==0
    count_elasn<-length(which(elasn_logi[1:length(elasn_logi)]))
    delj_count<-0
    for(j in 1:p){
      if(r[j]<alpha*lambda[k]){
        delj_count=delj_count+1
      }
    }
    rej_radio[k]<-delj_count/count_elasn
    scr_radio[k]<-1-delj_count/dim(X)[2]
  }
  return(list(rej_radio=rej_radio,scr_radio=scr_radio))
}

index<-function(X,Y,r,lambda,alpha){
  n=dim(X)[1]
  p=dim(X)[2]
  delj<-c()
  for(j in 1:p){
    if(r[j]<alpha*lambda){
      delj<-rbind(delj,j)
    }
  }
  ind<-rep(0,p)
  ind[delj]=1
  return(ind)
}

alam<-lammax(X,Y,100,0.25)
lambda<-alam$lambda
lambda0<-alam$lambda_min
lambda1<-alam$lambda_max
lam_radio<-lambda/lambda1

k=100
indexx<-matrix(0,ncol=k,nrow=p)
m<-rep(0,k)
for(i in 1:k){
  indexx[,i]<-index(X,Y,r,lambda[i],1)
  m[i]<-sum(indexx[,i])
}

alpha=0.4
Rn=100
totalfold=2
ci_coxlasso<-c(rep(0,Rn))
ci_coxelastic_lasso<-c(rep(0,Rn))
ci_coxsafe_elastic<-c(rep(0,Rn))
ibs_coxlasso<-c(rep(0,Rn))
ibs_coxelastic_lasso<-c(rep(0,Rn))
ibs_coxsafe_elastic<-c(rep(0,Rn))
data<-data.frame(Y,X)

lambda_best_index<-70
lambda.best<-lambda[lambda_best_index]*0.001/alpha

s.data<-data.frame(Y,X[,which(indexx[,lambda_best_index]!=1)])#after safe data

for (ii in seq(1, Rn, by=totalfold)){
  print(ii)
  set.seed(5)
  data<-data[sample(nrow(data)),]
  s.data<-s.data[sample(nrow(s.data)),]
  folds <- cut(seq(1,nrow(data)),breaks=totalfold,labels=FALSE)
  for(k in 1:totalfold){
    testIndexes <- which(folds==k,arr.ind=TRUE)
    y1<-data[,1:2]
    y2<-s.data[,1:2]
    names(y1)<-c("time","status")
    names(y2)<-c("time","status")
    data<-data.frame(y1,data[,-c(1:2)])
    ss.data<-data.frame(y2,s.data[,-c(1:2)])
    teset <- data[testIndexes, ]
    trset <- data[-testIndexes, ]
    s.teset <- ss.data[testIndexes, ]
    s.trset <- ss.data[-testIndexes, ]
    trset_r<-encox_safer(s.trset[,3:dim(s.trset)[2]],s.trset[,1:2])
    trset_index<-index(s.trset[,3:dim(s.trset)[2]],s.trset[,1:2],trset_r,lambda.best,alpha)
    s.trset<-s.trset[,which(trset_index!=1)]
    s.teset<-s.teset[,which(trset_index!=1)]
    print("coxlass")
    tryCatch({
      coxglmfit=cv.glmnet(x=as.matrix(trset[,-c(1:2)]),y=Surv(trset$time, trset$status),family="cox")
      coxglmpre=predict(coxglmfit,as.matrix(teset[,-c(1:2)]), type="response",s=lambda.best*0.4)
      ci_coxlasso[ii+k-1]=unlist(concordance.index(coxglmpre,teset$time,teset$status)[1])
      if(!is.na(ci_coxlasso[ii+k-1])){
        #cal ibs
        uniquetimes=sort(unique(trset$time))
        coxglmpre2=as.numeric(predict(coxglmfit,as.matrix(teset[,-c(1:2)]),s=lambda.best, type="link"))
        basesurvresult <- basesurv(Surv(teset$time,teset$status), coxglmpre2,uniquetimes)
        cox_surv_prob <- exp(exp(coxglmpre2) %*% -t(basesurvresult$cumBaseHaz))
        ibs_coxlasso[ii+k-1]=isbrier(Surv(teset$time,teset$status),t(cox_surv_prob),uniquetimes)
      }
      
      if (is.na(ci_coxlasso[ii+k-1]))
        ibs_coxlasso[ii+k-1]=NA
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    print("coxelastic_lasso")
    tryCatch({
      coxelafit=cv.glmnet(x=as.matrix(trset[,-c(1,2)]),y=Surv(trset$time, trset$status),family="cox",alpha=alpha)
      coxglmpre=predict(coxelafit,as.matrix(teset[,-c(1,2)]), type="response",s=lambda.best)
      ci_coxelastic_lasso[ii+k-1]=unlist(concordance.index(coxglmpre,teset$time,teset$status)[1])
      
      if(!is.na(ci_coxelastic_lasso[ii+k-1])){
        #cal ibs
        uniquetimes=sort(unique(trset$time))
        coxglmpre2=as.numeric(predict(coxelafit,as.matrix(teset[,-c(1,2)]),s=lambda.best, type="link"))
        basesurvresult <- basesurv(Surv(teset$time,teset$status), coxglmpre2,uniquetimes)
        cox_surv_prob  <- exp(exp(coxglmpre2) %*% -t(basesurvresult$cumBaseHaz))
        ibs_coxelastic_lasso[ii+k-1]=isbrier(Surv(teset$time,teset$status),t(cox_surv_prob),uniquetimes)
        
      }
      if (is.na(ci_coxelastic_lasso[ii+k-1]))
        ibs_coxelastic_lasso[ii+k-1]=NA
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    print("coxsafe_elastic")
    #logrank
    tryCatch({
      coxelafit=cv.glmnet(x=as.matrix(s.trset[,-c(1,2)]),y=Surv(s.trset$time, s.trset$status),family="cox",alpha=alpha)
      coxglmpre=predict(coxelafit,as.matrix(s.teset[,-c(1,2)]), type="response",s=lambda.best)
      ci_coxsafe_elastic[ii+k-1]=unlist(concordance.index(coxglmpre,s.teset$time,s.teset$status)[1])
      
      if(!is.na(ci_coxsafe_elastic[ii+k-1])){
        #cal ibs
        uniquetimes=sort(unique(trset$time))
        coxglmpre2=as.numeric(predict(coxelafit,as.matrix(s.teset[,-c(1,2)]),s=lambda.best, type="link"))
        basesurvresult <- basesurv(Surv(s.teset$time,s.teset$status), coxglmpre2,uniquetimes)
        cox_surv_prob  <- exp(exp(coxglmpre2) %*% -t(basesurvresult$cumBaseHaz))
        ibs_coxsafe_elastic[ii+k-1]=isbrier(Surv(s.teset$time,s.teset$status),t(cox_surv_prob),uniquetimes)
        
      }
      if (is.na(ci_coxsafe_elastic[ii+k-1]))
        ibs_coxsafe_elastic[ii+k-1]=NA
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

#######
CI<-c(ci_coxsafe_elastic,ci_coxelastic_lasso,ci_coxlasso)
CMethods<-c(rep("Safe_ENCox",100),rep("ENCox",100),rep("Lasso_Cox",100))
CMethod<- factor(CMethods, levels=c("Safe_ENCox","ENCox","Lasso_Cox"), ordered=TRUE)
xdataC<-data.frame(CI,CMethod)
##
mytheme<- theme(axis.title=element_text(face="plain",size=12, color="black"),
                axis.text=element_text(face="plain", size=12,color="black"), 
                panel.background=element_rect(color="white"),
                panel.grid.minor.y=element_blank(),
                panel.grid.minor.x=element_blank())

boxplot_CI<-ggplot(data=xdataC,aes(x=CMethod,y=CI,fill=CMethods))+
  geom_boxplot(width=0.5,outlier.colour="black", outlier.shape=1.5, outlier.size=3)+
  labs(title="",x="",y="C-index")+theme(legend.position="none")+mytheme+
  scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  theme(axis.title.y=element_text(size=16),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),legend.position="none")
boxplot_CI

###b-score
IBS<-c(ibs_coxsafe_elastic,ibs_coxelastic_lasso,ibs_coxlasso)
IBMethods<-c(rep("Safe_ENCox",100),rep("ENCox",100),rep("Lasso_Cox",100))
IBMethod<- factor(CMethods, levels=c("Safe_ENCox","ENCox","Lasso_Cox"), ordered=TRUE)
xdataIB<-data.frame(IBS,IBMethods)

boxplot_IBS<-ggplot(data=xdataIB,aes(x=IBMethod,y=IBS,fill=IBMethods))+
  geom_boxplot(width=0.5,outlier.colour="black", outlier.shape=1.5, outlier.size=3)+
  labs(title="",x="",y="B-score")+theme(legend.position="none")+mytheme+
  scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  theme(axis.title.y=element_text(size=16),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),legend.position="none")

boxplot_IBS

#####time
Rn=10
totalfold=2
for (ii in seq(1, Rn, by=totalfoldalfold)){
  print(ii)
  data<-data[sample(nrow(data)),]
  s.data<-s.data[sample(nrow(s.data)),]
  folds <- cut(seq(1,nrow(data)),breaks=totalfold,labels=FALSE)
  for(k in 1:totalfold){
    print(k)
    testIndexes <- which(folds==k,arr.ind=TRUE)
    y1<-data[,1:2]
    y2<-s.data[,1:2]
    names(y1)<-c("time","status")
    names(y2)<-c("time","status")
    data<-data.frame(y1,data[,-c(1:2)])
    ss.data<-data.frame(y2,s.data[,-c(1:2)])
    teset <- data[testIndexes, ]
    trset <- data[-testIndexes, ]
    s.teset <- ss.data[testIndexes, ]
    s.trset <- ss.data[-testIndexes, ]
    mbm<- microbenchmark(
      "Lasso_Cox"={
        coxglmfit=cv.glmnet(x=as.matrix(trset[,-c(1:2)]),y=Surv(trset$time, trset$status),family="cox",alpha=1)
        coxglmpre=predict(coxglmfit,as.matrix(teset[,-c(1:2)]), type="response",s=coxglmfit$lambda.min)
      },
      "ENCox"={
        coxglmfit=cv.glmnet(x=as.matrix(trset[,-c(1,2)]),Surv(trset$time, trset$status),type.measure = "deviance",nfold = 10,family = "cox",alpha = alpha)
        coxglmpre=predict(coxglmfit,as.matrix(teset[,-c(1,2)]), type="response",s=coxglmfit$lambda.min)
      },
      "Safe_ENCox"={
        coxglmfit=cv.glmnet(x=as.matrix(s.trset[,-c(1,2)]),Surv(s.trset$time,s.trset$status),type.measure = "deviance",nfold = 10,family = "cox",alpha =alpha)
        coxglmpre=predict(coxglmfit,as.matrix(s.teset[,-c(1,2)]), type="response",s=coxglmfit$lambda.min)
      }, times=100L)
  }
}

mbm0<-data.frame(mbm$expr,mbm$time/1000000000)
colnames(mbm0)<-c('expr','time')

p<-ggplot(mbm0, aes(x=expr, y = time))
p+geom_violin()+coord_flip()+labs(x='',y='Time[seconds]')