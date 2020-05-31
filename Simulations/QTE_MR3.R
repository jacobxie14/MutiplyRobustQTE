QTE_MR<-function(data,ToC=1,prob=0.5){
  #data<-Datasimulate(1000,setting="continuous",indep=TRUE,psm=1,outm="A")
  
  data<-as.data.frame(data)
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

bsc <- function (t, t0, itmax=100) {
  lo <- 1
  hi <- length(t)
  mid <- .5*(lo+hi)
  mid.lo <- floor(mid)
  mid.hi <- ceiling(mid)
  t.mid.lo <- t[mid.lo]
  t.mid.hi <- t[mid.hi]
  it <- 0
  while ((it<itmax)&((t.mid.lo>=t0)|(t.mid.hi<t0))) {
    if (t.mid.lo>=t0) hi <- mid.lo else lo <- mid.hi
    mid <- .5*(lo+hi)
    mid.lo <- floor(mid)
    mid.hi <- ceiling(mid)
    t.mid.lo <- t[mid.lo]
    t.mid.hi <- t[mid.hi]
    it <- it+1
  }
  mid.hi
}

# weighted quantiles
wtd.quantile <- function (y, weights, probs, normwt=TRUE) {
  od <- order(y)
  yo <- y[od]
  wo <- weights[od]
  if (normwt) wo <- wo/sum(wo)
  n <- length(wo)
  wos <- wo
  for (i in 2:n) wos[i] <- wos[i-1]+wo[i]
  np <- length(probs)
  res <- numeric(np)
  for (k in 1:np) {
    p <- probs[k]
    i.k <- bsc(wos,p)
    res[k] <- yo[i.k]
  }
  res
}

inwei<-function(ps,A,ATT=FALSE){
  w<-rep(NA,length(ps))
  if(!ATT){
    w[A==1]<-c(1/ps)[A==1]
    w[A==0]<-c(1/(1-ps))[A==0]
  }else{
    w[A==1]<-1
    w[A==0]<-c(ps/(1-ps))[A==0]
  }
  return(w)
}
setwd("/u5/y63xie/MacProfile/project3/datasim5/")
load("/u5/y63xie/MacProfile/project3/datasim5/simdata1000ind.Rdata")
setwd("/u5/y63xie/MacProfile/project3/datasim5/200/")
load("simdata200ind.Rdata")


library(parallel)
cl<-makeCluster(detectCores())
clusterExport(cl,c("QTE_MR","constrOptim","inwei","glm","wtd.quantile","bsc","uniroot"))
system.time(A1_QTE_q1<-parApply(cl,simdata1A,2,QTE_MR,prob=0.5,ToC=1))
system.time(A1_QTE_q0<-parApply(cl,simdata1A,2,QTE_MR,prob=0.5,ToC=0))
A1_QTE<-A1_QTE_q1-A1_QTE_q0
apply(A1_QTE,1,mean)
system.time(A2_QTE_q1<-parApply(cl,simdata2A,2,QTE_MR,prob=0.5,ToC=1))
system.time(A2_QTE_q0<-parApply(cl,simdata2A,2,QTE_MR,prob=0.5,ToC=0))
A2_QTE<-A2_QTE_q1-A2_QTE_q0
apply(A2_QTE,1,mean)
system.time(B1_QTE_q1<-parApply(cl,simdata1B,2,QTE_MR,prob=0.5,ToC=1))
system.time(B1_QTE_q0<-parApply(cl,simdata1B,2,QTE_MR,prob=0.5,ToC=0))
B1_QTE<-B1_QTE_q1-B1_QTE_q0
apply(B1_QTE,1,mean)
system.time(B2_QTE_q1<-parApply(cl,simdata2B,2,QTE_MR,prob=0.5,ToC=1))
system.time(B2_QTE_q0<-parApply(cl,simdata2B,2,QTE_MR,prob=0.5,ToC=0))
B2_QTE<-B2_QTE_q1-B2_QTE_q0
apply(B2_QTE,1,mean)

#######special 3C
system.time(C3_QTE_q1<-parApply(cl,simdata3C,2,QTE_MR,prob=0.5,ToC=1))
system.time(C3_QTE_q0<-parApply(cl,simdata3C,2,QTE_MR,prob=0.5,ToC=0))
C3_QTE<-C3_QTE_q1-C3_QTE_q0
apply(C3_QTE,1,mean)

apply(C3_QTE,1,sd)



save(A1_QTE,A1_QTE_q1,A1_QTE_q0,A2_QTE,A2_QTE_q1,A2_QTE_q0,B1_QTE,B1_QTE_q1,B1_QTE_q0,B2_QTE,B2_QTE_q1,B2_QTE_q0,
     file="/u5/y63xie/MacProfile/project3/datasim5/QTE_MR3.Rdata")
save(A1_QTE,A1_QTE_q1,A1_QTE_q0,A2_QTE,A2_QTE_q1,A2_QTE_q0,B1_QTE,B1_QTE_q1,B1_QTE_q0,B2_QTE,B2_QTE_q1,B2_QTE_q0,file="/u5/y63xie/MacProfile/project3/datasim5/data_R1&S5/QTE_MR3_size_5000.Rdata")
save(A1_QTE,A1_QTE_q1,A1_QTE_q0,A2_QTE,A2_QTE_q1,A2_QTE_q0,B1_QTE,B1_QTE_q1,B1_QTE_q0,B2_QTE,B2_QTE_q1,B2_QTE_q0,
     file="/u5/y63xie/MacProfile/project3/datasim5/200/QTE_MR3.Rdata")

setwd("/u5/y63xie/MacProfile/project3/datasim5/")
load("simdata5000ind.Rdata")



bias<-rbind(apply(A1_QTE,1,mean),apply(A2_QTE,1,mean),apply(B1_QTE,1,mean),apply(B2_QTE,1,mean))
colnames(bias)<-c("MR_AB","MR_12A","MR_12B","MR_1AB","MR_2AB","MR_12AB","OR_A","OR_B")
rownames(bias)<-c("1A","2A","1B","2B")
xtable(bias,digits=4)


STDV<-rbind(apply(A1_QTE,1,sd),apply(A2_QTE,1,sd),apply(B1_QTE,1,sd),apply(B2_QTE,1,sd))
colnames(STDV)<-c("MR_AB","MR_12A","MR_12B","MR_1AB","MR_2AB","MR_12AB","OR_A","OR_B")
rownames(STDV)<-c("1A","2A","1B","2B")
xtable(STDV,digits=4)

load("/Volumes/y63xie/MacProfile/project3/datasim3/true_q.Rdata")
b1<-apply(Q1,1,mean)[2]-apply(Q1,1,mean)[6]#1A
b2<-apply(Q2,1,mean)[2]-apply(Q2,1,mean)[6]#2A
b3<-apply(Q3,1,mean)[2]-apply(Q3,1,mean)[6]#1B
b4<-apply(Q4,1,mean)[2]-apply(Q4,1,mean)[6]#2B
q<-rbind(rep(b1,8),rep(b2,8),rep(b3,8),rep(b4,8))
MSE<-(bias-q)^2+STDV^2
xtable(sqrt(MSE),digits=4)


q1<-(apply(Q1[1:4,],1,mean)-apply(Q1[5:8,],1,mean))[2]
q2<-(apply(Q2[1:4,],1,mean)-apply(Q2[5:8,],1,mean))[2]
q3<-(apply(Q3[1:4,],1,mean)-apply(Q3[5:8,],1,mean))[2]
q4<-(apply(Q4[1:4,],1,mean)-apply(Q4[5:8,],1,mean))[2]
CI1<-function(i){
  qte1<-A1_QTE[i,]
  qte2<-A2_QTE[i,]
  qte3<-B1_QTE[i,]
  qte4<-B2_QTE[i,]
  sd1<-STDV[1,i]
  sd2<-STDV[2,i]
  sd3<-STDV[3,i]
  sd4<-STDV[4,i]
  id1<-(q1<=qte1+1.96*sd1 & q1>=qte1-1.96*sd1)
  id2<-(q2<=qte2+1.96*sd2 & q2>=qte2-1.96*sd2)
  id3<-(q3<=qte3+1.96*sd3 & q3>=qte3-1.96*sd3)
  id4<-(q4<=qte4+1.96*sd4 & q4>=qte4-1.96*sd4)
  c(mean(id1),mean(id2),mean(id3),mean(id4))
}

CR<-sapply(1:8,CI1)
colnames(CR)<-c("MR_AB","MR_12A","MR_12B","MR_1AB","MR_2AB","MR_12AB","OR_A","OR_B")
rownames(CR)<-c("1A","2A","1B","2B")
xtable(CR,digits=4)