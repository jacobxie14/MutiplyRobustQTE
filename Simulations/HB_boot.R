pe.hdr_boot <- function (ID,odata,prob,scenario,K=5, normalize=FALSE) {
  data<-odata[ID,]
  t<-data$z
  y1<-data$y[t==1];y0<-data$y[t==0]
  y<-data$y
  n <- length(data$y)
  if(scenario=="1A"){
    X.ps <-cbind(data$w1,data$w2,data$w4,data$w5,data$w7,data$w8)
    X.or <-as.matrix(data[,1:6])
  }else if(scenario=="2A"){
    X.ps <-cbind(-exp(data$w1/3),data$w2,data$w1-data$w4+3,data$w1*data$w5/10+0.5,data$w7,-exp(data$w8/2))
    X.or <-as.matrix(data[,1:6])
  }else if(scenario=="1B"){
    X.ps <-cbind(data$w1,data$w2,data$w4,data$w5,data$w7,data$w8)
    X.or <-cbind(data$w1,exp(data$w2/2),data$w3,exp(data$w4/3),data$w5*data$w6,data$w6,data$w1^2)
  }else if(scenario=="2B"){
    X.ps <-cbind(-exp(data$w1/3),data$w2,data$w1-data$w4+3,data$w1*data$w5/10+0.5,data$w7,-exp(data$w8/2))
    X.or <-cbind(data$w1,exp(data$w2/2),data$w3,exp(data$w4/3),data$w5*data$w6,data$w6,data$w1^2)
  }
  
  Model.x <- cbind(rep(1,n),X.or)
  MX0 <- Model.x[t==0,]
  MX1 <- Model.x[t==1,]
  #OR model
  fm0 <- lm(y0~0+MX0)
  bet0 <- fm0$coefficients
  fm1 <- lm(y1~0+MX1)
  bet1<- fm1$coefficients
  
  W.bet0 <- Model.x%*%bet0
  W.bet1 <- Model.x%*%bet1

  #ps model
  fm3 <- glm(t~X.ps-1,family=binomial)
  t.hat <- fm3$fitted.values
  t.hat.rnk <- rank(t.hat,ties.method="average")
  str <- factor(ceiling(K*t.hat.rnk/n))
  fm4 <- glm(t~str,family=binomial)
  t.hat<-fm4$fitted.values
  w0 <- ifelse(t==0,1/(1-t.hat),0)
  w1 <- ifelse(t==1,1/(t.hat),0)
  
  if (normalize) {
    w0 <- w0/sum(w0)
    w1 <- w1/sum(w1)
  }
  f0 <- function(q) {
    G <- pnorm(q-W.bet0)
    I <- as.numeric(y<=q)
    abs(mean((G-prob)+w0*(I-G)))
  }
  f1 <- function(q) {
    G <- pnorm(q-W.bet1)
    I <- as.numeric(y<=q)
    abs(mean((G-prob)+w1*(I-G)))
  }
  q0 <- quantile(y[t==0],probs=prob)
  q0 <- optim(q0,f0,lower=min(y[t==0]),upper=max(y[t==0]),method="Brent")$par
  q1 <- quantile(y[t==1],probs=prob)
  q1 <- optim(q1,f1,lower=min(y[t==1]),upper=max(y[t==1]),method="Brent")$par
  c(q0,q1,q1-q0)
}



HB_boot<-function(data,prob=0.5,scenario,B=500,normalize=FALSE){
  data<-as.data.frame(data)
  n<-dim(data)[1]
  ID<-matrix(sample(1:n,B*n,replace=TRUE),B,n)
  x<-matrix(0,3,dim(ID)[1])
  for(i in 1:dim(ID)[1]){
    x[,i]<-pe.hdr_boot(ID[i,],odata=data,prob=prob,scenario)
  }
  
  #system.time(dd<-apply(ID,1,pe.hdr_boot,odata=data,prob,scenario,normalize))
  bsd<-apply(x,1,sd)
  return(bsd)
}
library(parallel)
cl<-makeCluster(detectCores())
clusterExport(cl,c("pe.hdr_boot","optim","glm","HB_boot"))
p=0.25
set.seed(123456)
#data=1A
system.time(A1_QTE1<-parApply(cl,simdata1A,2,HB_boot,prob=p,scenario="1A",normalize=FALSE))
system.time(A1_QTE2<-parApply(cl,simdata1A,2,HB_boot,prob=p,scenario="2A",normalize=FALSE))
system.time(A1_QTE3<-parApply(cl,simdata1A,2,HB_boot,prob=p,scenario="1B",normalize=FALSE))
system.time(A1_QTE4<-parApply(cl,simdata1A,2,HB_boot,prob=p,scenario="2B",normalize=FALSE))#INCORRECT
apply(A1_QTE1,1,mean)
apply(A1_QTE2,1,mean)
apply(A1_QTE3,1,mean)
apply(A1_QTE4,1,mean)
#data=2A
system.time(A2_QTE1<-parApply(cl,simdata2A,2,HB_boot,prob=p,scenario="1A",normalize=FALSE))
system.time(A2_QTE2<-parApply(cl,simdata2A,2,HB_boot,prob=p,scenario="2A",normalize=FALSE))
system.time(A2_QTE3<-parApply(cl,simdata2A,2,HB_boot,prob=p,scenario="1B",normalize=FALSE))#INCORRECT
system.time(A2_QTE4<-parApply(cl,simdata2A,2,HB_boot,prob=p,scenario="2B",normalize=FALSE))
apply(A2_QTE1,1,mean)
apply(A2_QTE2,1,mean)
apply(A2_QTE3,1,mean)
apply(A2_QTE4,1,mean)
#data=2B
system.time(B2_QTE1<-parApply(cl,simdata2B,2,HB_boot,prob=p,scenario="1A",normalize=FALSE))#INCORRECT
system.time(B2_QTE2<-parApply(cl,simdata2B,2,HB_boot,prob=p,scenario="2A",normalize=FALSE))
system.time(B2_QTE3<-parApply(cl,simdata2B,2,HB_boot,prob=p,scenario="1B",normalize=FALSE))
system.time(B2_QTE4<-parApply(cl,simdata2B,2,HB_boot,prob=p,scenario="2B",normalize=FALSE))
apply(B2_QTE1,1,mean)
apply(B2_QTE2,1,mean)
apply(B2_QTE3,1,mean)
apply(B2_QTE4,1,mean)
#data=1B
system.time(B1_QTE1<-parApply(cl,simdata1B,2,HB_boot,prob=p,scenario="1A",normalize=FALSE))
system.time(B1_QTE2<-parApply(cl,simdata1B,2,HB_boot,prob=p,scenario="2A",normalize=FALSE))#INCORRECT
system.time(B1_QTE3<-parApply(cl,simdata1B,2,HB_boot,prob=p,scenario="1B",normalize=FALSE))
system.time(B1_QTE4<-parApply(cl,simdata1B,2,HB_boot,prob=p,scenario="2B",normalize=FALSE))
apply(B1_QTE1,1,mean)
apply(B1_QTE2,1,mean)
apply(B1_QTE3,1,mean)
apply(B1_QTE4,1,mean)

load("/u5/y63xie/MacProfile/project3/datasim5/simdata1000_3C.Rdata")
load("/u5/y63xie/MacProfile/project3/datasim5/200/simdata200_3C.Rdata")
set.seed(123456)
p=0.25
#data=3C
system.time(C3_QTE1_sd<-parApply(cl,simdata3C,2,HB_boot,prob=p,scenario="1A",normalize=FALSE))
system.time(C3_QTE2_sd<-parApply(cl,simdata3C,2,HB_boot,prob=p,scenario="2A",normalize=FALSE))
system.time(C3_QTE3_sd<-parApply(cl,simdata3C,2,HB_boot,prob=p,scenario="1B",normalize=FALSE))
system.time(C3_QTE4_sd<-parApply(cl,simdata3C,2,HB_boot,prob=p,scenario="2B",normalize=FALSE))
apply(C3_QTE1_sd,1,mean)
apply(C3_QTE2_sd,1,mean)
apply(C3_QTE3_sd,1,mean)
apply(C3_QTE4_sd,1,mean)
save(C3_QTE1_sd,C3_QTE2_sd,C3_QTE3_sd,C3_QTE4_sd,file="BSD_HB_C3_0.25.Rdata")


save(A1_QTE1,A1_QTE2,A1_QTE3,A1_QTE4,B1_QTE1,B1_QTE2,B1_QTE3,B1_QTE4,
     A2_QTE1,A2_QTE2,A2_QTE3,A2_QTE4,B2_QTE1,B2_QTE2,B2_QTE3,B2_QTE4,
     file="BSD_HB_0.25.Rdata")
save(B1_QTE1,B1_QTE2,B1_QTE3,B1_QTE4,file="/u5/y63xie/MacProfile/project3/datasim5/bootstrap/BSD_HB_0.95.Rdata")


load("N:/MacProfile/project3/datasim5/bootstrap/BSD_HB.Rdata")
load("BSD_HB_0.25.Rdata")
load("BSD_HB_0.75.Rdata")
load("BSD_HB_0.95.Rdata")
load("BSD_HB.Rdata")
a1<-c(apply(A1_QTE1,1,mean)[3],apply(A1_QTE2,1,mean)[3],apply(A1_QTE3,1,mean)[3],apply(A1_QTE4,1,mean)[3])
b1<-c(apply(A2_QTE1,1,mean)[3],apply(A2_QTE2,1,mean)[3],apply(A2_QTE3,1,mean)[3],apply(A2_QTE4,1,mean)[3])
c1<-c(apply(B1_QTE1,1,mean)[3],apply(B1_QTE2,1,mean)[3],apply(B1_QTE3,1,mean)[3],apply(B1_QTE4,1,mean)[3])
d1<-c(apply(B2_QTE1,1,mean)[3],apply(B2_QTE2,1,mean)[3],apply(B2_QTE3,1,mean)[3],apply(B2_QTE4,1,mean)[3])
bsd<-rbind(a1,b1,c1,d1)
colnames(bsd)<-c("HB_1A","HB_2A","HB_1B","HB_2B")
rownames(bsd)<-c("1A","2A","1B","2B")
xtable(bsd,digits=4)

bsd1<-rbind(A1_QTE1[3,],A2_QTE1[3,],B1_QTE1[3,],B2_QTE1[3,],
            A1_QTE2[3,],A2_QTE2[3,],B1_QTE2[3,],B2_QTE2[3,],
            A1_QTE3[3,],A2_QTE3[3,],B1_QTE3[3,],B2_QTE3[3,],
            A1_QTE4[3,],A2_QTE4[3,],B1_QTE4[3,],B2_QTE4[3,])
load("N:/MacProfile/project3/datasim5/q_0.95/HB.Rdata")
load("N:/MacProfile/project3/datasim5/200/HB.Rdata")
a<-c(apply(A1_QTE1,1,mean)[3],apply(A1_QTE2,1,mean)[3],apply(A1_QTE3,1,mean)[3],apply(A1_QTE4,1,mean)[3])
b<-c(apply(A2_QTE1,1,mean)[3],apply(A2_QTE2,1,mean)[3],apply(A2_QTE3,1,mean)[3],apply(A2_QTE4,1,mean)[3])
c<-c(apply(B1_QTE1,1,mean)[3],apply(B1_QTE2,1,mean)[3],apply(B1_QTE3,1,mean)[3],apply(B1_QTE4,1,mean)[3])
d<-c(apply(B2_QTE1,1,mean)[3],apply(B2_QTE2,1,mean)[3],apply(B2_QTE3,1,mean)[3],apply(B2_QTE4,1,mean)[3])
bias<-rbind(a,b,c,d)
colnames(bias)<-c("HB_1A","HB_2A","HB_1B","HB_2B")
rownames(bias)<-c("1A","2A","1B","2B")
xtable(bias,digits=4)

load("N:/MacProfile/project3/datasim3/true_q.Rdata")
i=2
q1<-apply(Q1,1,mean)[i]-apply(Q1,1,mean)[i+4]#1A
q2<-apply(Q2,1,mean)[i]-apply(Q2,1,mean)[i+4]#2A
q3<-apply(Q3,1,mean)[i]-apply(Q3,1,mean)[i+4]#1B
q4<-apply(Q4,1,mean)[i]-apply(Q4,1,mean)[i+4]#2B

CI2<-function(a=A1_QTE1[3,],b=A2_QTE1[3,],c=B1_QTE1[3,],d=B2_QTE1[3,],i){
  sd1<-bsd1[i+1,]
  sd2<-bsd1[i+2,]
  sd3<-bsd1[i+3,]
  sd4<-bsd1[i+4,]
  id1<-(q1<=a+1.96*sd1 & q1>=a-1.96*sd1)
  id2<-(q2<=b+1.96*sd2 & q2>=b-1.96*sd2)
  id3<-(q3<=c+1.96*sd3 & q3>=c-1.96*sd3)
  id4<-(q4<=d+1.96*sd4 & q4>=d-1.96*sd4)
  c(mean(id1),mean(id2),mean(id3),mean(id4))
}
cr1<-CI2(a=A1_QTE1[3,],b=A2_QTE1[3,],c=B1_QTE1[3,],d=B2_QTE1[3,],0)
cr2<-CI2(a=A1_QTE2[3,],b=A2_QTE2[3,],c=B1_QTE2[3,],d=B2_QTE2[3,],4)
cr3<-CI2(a=A1_QTE3[3,],b=A2_QTE3[3,],c=B1_QTE3[3,],d=B2_QTE3[3,],8)
cr4<-CI2(a=A1_QTE4[3,],b=A2_QTE4[3,],c=B1_QTE4[3,],d=B2_QTE4[3,],12)
cr<-cbind(cr1,cr2,cr3,cr4)
colnames(cr)<-c("HB_1A","HB_2A","HB_1B","HB_2B")
rownames(cr)<-c("1A","2A","1B","2B")
xtable(cr,digits=4)

