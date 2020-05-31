#QTE_MR3 
#Bootstrapping
QTE_MR_BOOT<-function(ID,odata,ToC=1,prob=0.5){
  #data<-Datasimulate(1000,setting="continuous",indep=TRUE,psm=1,outm="A")
  
  data<-odata[ID,]
  y1<-data$y1;y0<-data$y0
  data<-data[,-c(11,12)]
  n<-dim(data)[1];n1<-sum(data$z)
  n0<-n-n1
  if(ToC==1){
    m=n1
  }else{
    m=n0
  }
  Treat<-data$z
  #PS models
  ps1<-glm(z~w1 + w2 + w4 + w5 + w7 + w8-1,data=data,family=binomial)
  ps2<-glm(z~I(-exp(w1/3)) + w2 + I(w1-w4+3) + I(w1*w5/10+0.5) + w7 + I(-exp(w8/2)) -1,data=data,family=binomial)
  
  #Outcome models
  lr1<-lm(y~w1 + w2 + w3 + w4 + w5 + w6,data=data[data$z==ToC,])
  lr2<-lm(y~w1 + I(exp(w2/2)) + w3 + I(exp(w4/3)) + I(w5*w6) + w6 + I(w1^2),data=data[data$z==ToC,])
  
  beta1<-lr1$coefficients
  beta2<-lr2$coefficients
  
  model.X1<-as.matrix(cbind(rep(1,n),data[,1:6]))
  model.X2<-as.matrix(cbind(rep(1,n),data$w1,exp(data$w2/2),data$w3,exp(data$w4/3),data$w5*data$w6,data$w6,data$w1^2))
  
  #outcome model1
  ocm1<-function(x){
    G<-pnorm(x-as.numeric(model.X1%*%beta1))
    mean(G)-prob
  }
  #q0<-quantile(data$y[data$z==0],0.5)
  lb <- min(as.numeric(model.X1%*%beta1))+qnorm(prob)
  ub <- max(as.numeric(model.X1%*%beta1))+qnorm(prob)
  q_1<-uniroot(ocm1,interval=c(lb,ub))$root
  #uniroot(ocm1,interval=c(-10,15))$root
  #outcome model 2
  ocm2<-function(x){
    G<-pnorm(x-as.numeric(model.X2%*%beta2))
    mean(G)-prob
  }
  #q0<-quantile(data$y[data$z==0],0.5)
  lb <- min(as.numeric(model.X2%*%beta2))+qnorm(prob)
  ub <- max(as.numeric(model.X2%*%beta2))+qnorm(prob)
  q_2<-uniroot(ocm2,interval=c(lb,ub))$root
  G1.i<-pnorm(q_1-as.numeric(model.X1%*%beta1))
  mG1<-sum(G1.i)/n
  G2.i<-pnorm(q_2-as.numeric(model.X2%*%beta2))
  mG2<-sum(G2.i)/n
  G3.i<-ps1$fitted.values
  mG3<-sum(G3.i)/n
  G4.i<-ps2$fitted.values
  mG4<-sum(G4.i)/n
  
  ########################################MR with 0011
  LG_0011<-function(r){
    # G1.i<-pnorm(q_1-as.numeric(model.X1%*%beta1))
    # mG1<-sum(G1.i)/n
    # G2.i<-pnorm(q_2-as.numeric(model.X2%*%beta2))
    # mG2<-sum(G2.i)/n
    -sum( log( 1+r[1]*(mG1-G1.i[Treat==ToC])+r[2]*(mG2-G2.i[Treat==ToC]) ) ) 
  }
  #opt_0011<-optim(c(0,0), LG_0011,control=list("maxit"=100,"trace"=FALSE), method="L-BFGS-B")
  constr0<-cbind(mG1-G1.i[Treat==ToC],mG2-G2.i[Treat==ToC])
  opt_0011<-constrOptim(c(0,0), LG_0011,grad=NULL,ui=constr0,ci=rep(-1,m),control=list("maxit"=100,"trace"=FALSE))
  r0<-opt_0011$par
  
  w0<-1/(1+r0[1]*(mG1-G1.i[Treat==ToC])+r0[2]*(mG2-G2.i[Treat==ToC]))/m
  w0<-w0/sum(w0)
  #wtd.quantile(data$y[data$z==1],weights=w0,prob)
  #q_2
  #q_true
  # f3.c<-function(q){
  #   G.i<-pnorm(q-as.numeric(model.X%*%beta0))
  #   sum(G.i[Treat==ToC]*w2)-prob
  # }
  ########################################MR with 1110
  LG_1110<-function(r){
    
    -sum( log( 1+r[1]*(mG1-G1.i[Treat==ToC])+r[2]*(mG3-G3.i[Treat==ToC])+r[3]*(mG4-G4.i[Treat==ToC]) ) ) 
  }
  #opt_1110<-optim(c(0,0,0), LG_1110,control=list("maxit"=100,"trace"=FALSE), method="BFGS")
  constr1<-cbind(mG1-G1.i[Treat==ToC],mG3-G3.i[Treat==ToC],mG4-G4.i[Treat==ToC])
  opt_1110<-constrOptim(c(0,0,0), LG_1110,grad=NULL,ui=constr1,ci=rep(-1,m),control=list("maxit"=100,"trace"=FALSE))
  r1<-opt_1110$par
  
  w1<-1/(1+r1[1]*(mG1-G1.i[Treat==ToC])+r1[2]*(mG3-G3.i[Treat==ToC])+r1[3]*(mG4-G4.i[Treat==ToC]))/m
  w1<-w1/sum(w1)
  ########################################MR with 1101
  LG_1101<-function(r){
    
    -sum( log( 1+r[1]*(mG2-G2.i[Treat==ToC])+r[2]*(mG3-G3.i[Treat==ToC])+r[3]*(mG4-G4.i[Treat==ToC]) ) ) 
  }
  #opt_1101<-optim(c(0,0,0), LG_1101,control=list("maxit"=100,"trace"=FALSE), method="BFGS")
  constr2<-cbind(mG2-G2.i[Treat==ToC],mG3-G3.i[Treat==ToC],mG4-G4.i[Treat==ToC])
  opt_1101<-constrOptim(c(0,0,0), LG_1101,grad=NULL,ui=constr2,ci=rep(-1,m),control=list("maxit"=100,"trace"=FALSE))
  r2<-opt_1101$par
  
  w2<-1/(1+r2[1]*(mG2-G2.i[Treat==ToC])+r2[2]*(mG3-G3.i[Treat==ToC])+r2[3]*(mG4-G4.i[Treat==ToC]))/m
  w2<-w2/sum(w2)
  ########################################MR with 1011
  LG_1011<-function(r){
    
    -sum( log( 1+r[1]*(mG1-G1.i[Treat==ToC])+r[2]*(mG2-G2.i[Treat==ToC])+r[3]*(mG3-G3.i[Treat==ToC]) ) ) 
  }
  #opt_1011<-optim(c(0,0,0), LG_1011,control=list("maxit"=100,"trace"=FALSE), method="BFGS")
  constr3<-cbind(mG1-G1.i[Treat==ToC],mG2-G2.i[Treat==ToC],mG3-G3.i[Treat==ToC])
  opt_1011<-constrOptim(c(0,0,0), LG_1011,grad=NULL,ui=constr3,ci=rep(-1,m),control=list("maxit"=100,"trace"=FALSE))
  r3<-opt_1011$par
  
  w3<-1/(1+r3[1]*(mG1-G1.i[Treat==ToC])+r3[2]*(mG2-G2.i[Treat==ToC])+r3[3]*(mG3-G3.i[Treat==ToC]))/m
  w3<-w3/sum(w3)
  ########################################MR with 0111
  LG_0111<-function(r){
    
    -sum( log( 1+r[1]*(mG1-G1.i[Treat==ToC])+r[2]*(mG2-G2.i[Treat==ToC])+r[3]*(mG4-G4.i[Treat==ToC]) ) ) 
  }
  #opt_0111<-optim(c(0,0,0), LG_0111,control=list("maxit"=100,"trace"=FALSE), method="BFGS")
  constr4<-cbind(mG1-G1.i[Treat==ToC],mG2-G2.i[Treat==ToC],mG4-G4.i[Treat==ToC])
  opt_0111<-constrOptim(c(0,0,0), LG_0111,grad=NULL,ui=constr4,ci=rep(-1,m),control=list("maxit"=100,"trace"=FALSE))
  r4<-opt_0111$par
  
  w4<-1/(1+r4[1]*(mG1-G1.i[Treat==ToC])+r4[2]*(mG2-G2.i[Treat==ToC])+r4[3]*(mG4-G4.i[Treat==ToC]))/m
  w4<-w4/sum(w4)
  ########################################MR with 1111
  LG_1111<-function(r){
    
    -sum( log( 1+r[1]*(mG1-G1.i[Treat==ToC])+r[2]*(mG2-G2.i[Treat==ToC])+r[3]*(mG3-G3.i[Treat==ToC])+r[4]*(mG4-G4.i[Treat==ToC]) ) ) 
  }
  #opt_1111<-optim(c(0,0,0,0), LG_1111,control=list("maxit"=100,"trace"=FALSE), method="BFGS")
  constr5<-cbind(mG1-G1.i[Treat==ToC],mG2-G2.i[Treat==ToC],mG3-G3.i[Treat==ToC],mG4-G4.i[Treat==ToC])
  opt_1111<-constrOptim(c(0,0,0,0), LG_1111,grad=NULL,ui=constr5,ci=rep(-1,m),control=list("maxit"=100,"trace"=FALSE))
  r5<-opt_1111$par
  
  w5<-1/(1+r5[1]*(mG1-G1.i[Treat==ToC])+r5[2]*(mG2-G2.i[Treat==ToC])+r5[3]*(mG3-G3.i[Treat==ToC])+r5[4]*(mG4-G4.i[Treat==ToC]))/m
  w5<-w5/sum(w5)
  #Our method
  q_w0<-wtd.quantile(data$y[data$z==ToC],weights=w0,prob)
  q_w1<-wtd.quantile(data$y[data$z==ToC],weights=w1,prob)
  q_w2<-wtd.quantile(data$y[data$z==ToC],weights=w2,prob)
  q_w3<-wtd.quantile(data$y[data$z==ToC],weights=w3,prob)
  q_w4<-wtd.quantile(data$y[data$z==ToC],weights=w4,prob)
  q_w5<-wtd.quantile(data$y[data$z==ToC],weights=w5,prob)
  return(c(q_w0,q_w1,q_w2,q_w3,q_w4,q_w5,q_1,q_2))
  
}

MR_Boot<-function(data,prob=0.5,B=500){
  data<-as.data.frame(data)
  n<-dim(data)[1]
  ID<-matrix(sample(1:n,B*n,replace=TRUE),B,n)
  system.time(dd1<-apply(ID,1,QTE_MR_BOOT,odata=data,prob=prob,ToC=1))
  system.time(dd0<-apply(ID,1,QTE_MR_BOOT,odata=data,prob=prob,ToC=0))
  bsd<-apply(dd1-dd0,1,sd)
  return(bsd)
}





setwd("/u5/y63xie/MacProfile/project3/datasim5/bootstrap")
library(parallel)
cl<-makeCluster(detectCores())
clusterExport(cl,c("MR_Boot","QTE_MR_BOOT","constrOptim","inwei","glm","wtd.quantile","bsc","uniroot"))
set.seed(12345)
p=0.95
system.time(A1_QTE_sd<-parApply(cl,simdata1A,2,MR_Boot,prob=p))
apply(A1_QTE_sd,1,mean)
system.time(A2_QTE_sd<-parApply(cl,simdata2A,2,MR_Boot,prob=p))
apply(A2_QTE_sd,1,mean)
system.time(B1_QTE_sd<-parApply(cl,simdata1B,2,MR_Boot,prob=p))
apply(B1_QTE_sd,1,mean)
system.time(B2_QTE_sd<-parApply(cl,simdata2B,2,MR_Boot,prob=p))
apply(B2_QTE_sd,1,mean)
save(A1_QTE_sd,A2_QTE_sd,B1_QTE_sd,B2_QTE_sd,file="BSD_MR_0.95.Rdata")

load("/u5/y63xie/MacProfile/project3/datasim5/simdata1000_3C.Rdata")
load("/u5/y63xie/MacProfile/project3/datasim5/200/simdata200_3C.Rdata")
set.seed(12345)
p=0.5
system.time(C3_QTE_sd<-parApply(cl,simdata3C,2,MR_Boot,prob=p))
apply(C3_QTE_sd,1,mean)
save(C3_QTE_sd,file="/u5/y63xie/MacProfile/project3/datasim5/BSD_MR_C3.Rdata")


load("N:/MacProfile/project3/datasim5/200/QTE_MR3.Rdata")
load("/u5/y63xie/MacProfile/project3/datasim5/q_0.95/QTE_MR3.Rdata")
bias<-rbind(apply(A1_QTE,1,mean),apply(A2_QTE,1,mean),apply(B1_QTE,1,mean),apply(B2_QTE,1,mean))[,c(2:6)]
colnames(bias)<-c("MR_12A","MR_12B","MR_1AB","MR_2AB","MR_12AB")
rownames(bias)<-c("1A","2A","1B","2B")
xtable(bias,digits=4)

load("N:/MacProfile/project3/datasim5/200/BSD_MR.Rdata")
BSD<-rbind(apply(A1_QTE_sd,1,mean),apply(A2_QTE_sd,1,mean),apply(B1_QTE_sd,1,mean),apply(B2_QTE_sd,1,mean))[,c(2:6)]
colnames(BSD)<-c("MR_12A","MR_12B","MR_1AB","MR_2AB","MR_12AB")
rownames(BSD)<-c("1A","2A","1B","2B")
xtable(BSD,digits=4)

load("N:/MacProfile/project3/datasim3/true_q.Rdata")
q1<-(apply(Q1[1:4,],1,mean)-apply(Q1[5:8,],1,mean))[2]
q2<-(apply(Q2[1:4,],1,mean)-apply(Q2[5:8,],1,mean))[2]
q3<-(apply(Q3[1:4,],1,mean)-apply(Q3[5:8,],1,mean))[2]
q4<-(apply(Q4[1:4,],1,mean)-apply(Q4[5:8,],1,mean))[2]
CI1<-function(i){
  qte1<-A1_QTE[i+1,]
  qte2<-A2_QTE[i+1,]
  qte3<-B1_QTE[i+1,]
  qte4<-B2_QTE[i+1,]
  sd1<-A1_QTE_sd[i+1,]
  sd2<-A2_QTE_sd[i+1,]
  sd3<-B1_QTE_sd[i+1,]
  sd4<-B2_QTE_sd[i+1,]
  id1<-(q1<=qte1+1.96*sd1 & q1>=qte1-1.96*sd1)
  id2<-(q2<=qte2+1.96*sd2 & q2>=qte2-1.96*sd2)
  id3<-(q3<=qte3+1.96*sd3 & q3>=qte3-1.96*sd3)
  id4<-(q4<=qte4+1.96*sd4 & q4>=qte4-1.96*sd4)
  c(mean(id1),mean(id2),mean(id3),mean(id4))
}

CR<-sapply(1:5,CI1)
colnames(CR)<-c("MR_12A","MR_12B","MR_1AB","MR_2AB","MR_12AB")
rownames(CR)<-c("1A","2A","1B","2B")
xtable(CR,digits=4)
save(CR,file="N:/MacProfile/project3/datasim5/CR.Rdata")
