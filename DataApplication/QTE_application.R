bweight<-read.csv("bweight.csv",header=TRUE)

bweight2 = bweight[,names(bweight)!=c("CigsPerDay")]
data=bweight2
hist(bweight$Weight,breaks=100,prob=T)
mean(bweight$MomSmoke)

plot(bweight$MomAge,bweight$Weight)
plot(bweight$MomWtGain,bweight$Weight)
sw<-(bweight$Weight-mean(bweight$Weight))/sd(bweight$Weight)
plot(bweight$MomAge,sw)
plot(bweight$MomWtGain,sw)


a<-data$Weight[data$MomSmoke==1]
b<-data$Weight[data$MomSmoke==0]
fa<-ecdf(a)
fb<-ecdf(b)
ya<-fa(sort(a))
yb<-fb(sort(b))
plot(sort(a),ya,type="l",lty=2,lwd=2,xlab="Weight",ylab="Probability")
lines(sort(b),yb,type="l",lty=3,lwd=2,col="red")
legend("topleft",legend=c("Smoked group","Non-smoked group"),lty=c(2,3),lwd=c(2,2),col=c("black","red"),cex=0.8)

QTE_MR<-function(data,ToC=1,prob=0.5){
  
  data<-as.data.frame(data) #bweight2
  n<-dim(data)[1];n1<-sum(data$MomSmoke)
  n0<-n-n1
  if(ToC==1){
    m=n1
  }else{
    m=n0
  }
  Treat<-data$MomSmoke
  #PS models
  #logit link
  ps1<-glm(MomSmoke~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel),
           data=data,family=binomial)
  ps2<-glm(MomSmoke~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
                    + I(MomAge^2), data=data,family=binomial)
  ps3<-glm(MomSmoke~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
             + I(MomAge^2)+ Black:MomAge+ Married:MomAge+ Boy:MomAge+ 
             as.factor(MomEdLevel):MomAge, data=data,family=binomial)

  #probit link
  ps4<-glm(MomSmoke~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel),
           data=data,family=binomial(link = "probit"))
  ps5<-glm(MomSmoke~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
             + I(MomAge^2), data=data,family=binomial(link = "probit"))
  ps6<-glm(MomSmoke~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
             + I(MomAge^2)+ Black:MomAge+ Married:MomAge+ Boy:MomAge+ 
             as.factor(MomEdLevel):MomAge, data=data,family=binomial(link = "probit"))
  

  #cloglog link
  ps7<-glm(MomSmoke~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel),
           data=data,family=binomial(link = "cloglog"))
  ps8<-glm(MomSmoke~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
             + I(MomAge^2), data=data,family=binomial(link = "cloglog"))
  ps9<-glm(MomSmoke~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
             + I(MomAge^2)+ Black:MomAge+ Married:MomAge+ Boy:MomAge+ 
             as.factor(MomEdLevel):MomAge, data=data,family=binomial(link = "cloglog"))
  
  #Outcome models
  lr1<-lm(Weight~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel),
                 data=data[data$MomSmoke==ToC,])
  lr2<-lm(Weight~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
                + I(MomAge^2),data=data[data$MomSmoke==ToC,])
  lr3<-lm(Weight~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
            + I(MomAge^2)+ Black:MomAge+ Married:MomAge+ Boy:MomAge+ 
            as.factor(MomEdLevel):MomAge, data=data[data$MomSmoke==ToC,])
  
  beta1<-lr1$coefficients
  beta2<-lr2$coefficients
  beta3<-lr3$coefficients
  
  sig1<-sqrt(sum(lr1$residuals^2)/lr1$df.residual)
  sig2<-sqrt(sum(lr2$residuals^2)/lr2$df.residual)
  sig3<-sqrt(sum(lr3$residuals^2)/lr3$df.residual)
  
  model.X1<-model.matrix(~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel),data=data)
  model.X2<-model.matrix(~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
                           + I(MomAge^2),data=data)
  model.X3<-model.matrix(~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
                           + I(MomAge^2)+ Black:MomAge+ Married:MomAge+ Boy:MomAge+ 
                           as.factor(MomEdLevel):MomAge,data=data)
  
  #outcome model1
  ocm1<-function(x){
    G<-pnorm((x-as.numeric(model.X1%*%beta1))/sig1)
    mean(G)-prob
  }
  #q0<-quantile(data$y[data$z==0],0.5)
  lb <- min(as.numeric(model.X1%*%beta1))+sig1*qnorm(prob)
  ub <- max(as.numeric(model.X1%*%beta1))+sig1*qnorm(prob)
  q_1<-uniroot(ocm1,interval=c(lb,ub))$root
  #uniroot(ocm1,interval=c(-10,15))$root
  #outcome model 2
  ocm2<-function(x){
    G<-pnorm((x-as.numeric(model.X2%*%beta2))/sig2)
    mean(G)-prob
  }
  #q0<-quantile(data$y[data$z==0],0.5)
  lb <- min(as.numeric(model.X2%*%beta2))+sig2*qnorm(prob)
  ub <- max(as.numeric(model.X2%*%beta2))+sig2*qnorm(prob)
  q_2<-uniroot(ocm2,interval=c(lb,ub))$root
  #outcome model3
  ocm3<-function(x){
    G<-pnorm((x-as.numeric(model.X3%*%beta3))/sig3)
    mean(G)-prob
  }
  #q0<-quantile(data$y[data$z==0],0.5)
  lb <- min(as.numeric(model.X3%*%beta3))+sig3*qnorm(prob)
  ub <- max(as.numeric(model.X3%*%beta3))+sig3*qnorm(prob)
  q_3<-uniroot(ocm3,interval=c(lb,ub))$root
  #uniroot(ocm3,interval=c(-10,15))$root
  

  
  G1.i<-ps1$fitted.values
  mG1<-sum(G1.i)/n
  G2.i<-ps2$fitted.values
  mG2<-sum(G2.i)/n
  G3.i<-ps3$fitted.values
  mG3<-sum(G3.i)/n
  
  G4.i<-ps4$fitted.values
  mG4<-sum(G4.i)/n
  G5.i<-ps5$fitted.values
  mG5<-sum(G5.i)/n
  G6.i<-ps6$fitted.values
  mG6<-sum(G6.i)/n
  
  G7.i<-ps7$fitted.values
  mG7<-sum(G7.i)/n
  G8.i<-ps8$fitted.values
  mG8<-sum(G8.i)/n
  G9.i<-ps9$fitted.values
  mG9<-sum(G9.i)/n

  
  G10.i<-pnorm((q_1-as.numeric(model.X1%*%beta1))/sig1)
  mG10<-sum(G10.i)/n
  G11.i<-pnorm((q_2-as.numeric(model.X2%*%beta2))/sig2)
  mG11<-sum(G11.i)/n
  G12.i<-pnorm((q_3-as.numeric(model.X3%*%beta3))/sig3)
  mG12<-sum(G12.i)/n
  
  ########################################MR with all models
  LG_all<-function(r){
    
    -sum( log( 1+r[1]*(mG1-G1.i[Treat==ToC])+r[2]*(mG2-G2.i[Treat==ToC])+r[3]*(mG3-G3.i[Treat==ToC])
                +r[4]*(mG4-G4.i[Treat==ToC])+r[5]*(mG5-G5.i[Treat==ToC])+r[6]*(mG6-G6.i[Treat==ToC])
                +r[7]*(mG7-G7.i[Treat==ToC])+r[8]*(mG8-G8.i[Treat==ToC])+r[9]*(mG9-G9.i[Treat==ToC])
                +r[10]*(mG10-G10.i[Treat==ToC])+r[11]*(mG11-G11.i[Treat==ToC])+r[12]*(mG12-G12.i[Treat==ToC])) ) 
  }
  #opt_1111<-optim(c(0,0,0,0), LG_1111,control=list("maxit"=100,"trace"=FALSE), method="BFGS")
  constr<-cbind(mG1-G1.i[Treat==ToC],mG2-G2.i[Treat==ToC],mG3-G3.i[Treat==ToC],mG4-G4.i[Treat==ToC],
                 mG5-G5.i[Treat==ToC],mG6-G6.i[Treat==ToC],mG7-G7.i[Treat==ToC],mG8-G8.i[Treat==ToC],
                 mG9-G9.i[Treat==ToC],mG10-G10.i[Treat==ToC],mG11-G11.i[Treat==ToC],mG12-G12.i[Treat==ToC])
  opt_all<-constrOptim(rep(0,12), LG_all,grad=NULL,ui=constr,ci=rep(-1,m),control=list("maxit"=100,"trace"=FALSE))
  r<-opt_all$par
  
  w<-1/(1+r[1]*(mG1-G1.i[Treat==ToC])+r[2]*(mG2-G2.i[Treat==ToC])+r[3]*(mG3-G3.i[Treat==ToC])
         +r[4]*(mG4-G4.i[Treat==ToC])+r[5]*(mG5-G5.i[Treat==ToC])+r[6]*(mG6-G6.i[Treat==ToC])
         +r[7]*(mG7-G7.i[Treat==ToC])+r[8]*(mG8-G8.i[Treat==ToC])+r[9]*(mG9-G9.i[Treat==ToC])
         +r[10]*(mG10-G10.i[Treat==ToC])+r[11]*(mG11-G11.i[Treat==ToC])+r[12]*(mG12-G12.i[Treat==ToC]))/m
  w<-w/sum(w)
  #Our method
  q_w<-wtd.quantile(data$Weight[data$MomSmoke==ToC],weights=w,prob)

  return(q_w)
  
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

p = seq(0.05,0.95,by=0.05)
qte<-numeric(length(p))
for (i in 1:length(p)){
 q1<-QTE_MR(data=bweight2,ToC=1,prob=p[i])
 q0<-QTE_MR(data=bweight2,ToC=0,prob=p[i])
  qte[i] = q1-q0
    
}


n<-dim(bweight2)[1]
QTE_boot<-function(id,x,prob=0.5,ToC){
  x<-x[id,]
  rownames(x)=NULL
  out<-QTE_MR(data=x,ToC=ToC,prob=prob)
  return(out)
}

set.seed(827)
ID<-sample(1:n,500*n,replace=TRUE)
ID<-matrix(ID,500,n)


library(parallel)
cl<-makeCluster(detectCores())
clusterExport(cl,c("QTE_MR","constrOptim","inwei","glm","wtd.quantile","bsc","uniroot","QTE_boot"))
system.time(QTE_q1<-parApply(cl,ID,1,QTE_boot,bweight2,prob=0.5,ToC=1))
system.time(QTE_q0<-parApply(cl,ID,1,QTE_boot,bweight2,prob=0.5,ToC=0))
QTE<-QTE_q1-QTE_q0
save(QTE,QTE_q1,QTE_q0,file="QTE_0.5.Rdata")

q<-matrix(0,19,500)
QTE_q1<-q
QTE_q0<-q
for (i in 1:length(p)){
 QTE_q1[i,]<-parApply(cl,ID,1,QTE_boot,bweight2,prob=p[i],ToC=1)
 QTE_q0[i,]<-parApply(cl,ID,1,QTE_boot,bweight2,prob=p[i],ToC=0)
 q[i,]<-QTE_q1[i,]-QTE_q0[i,]
 if(i%%5==0){print(i)}
}
save(q,QTE_q1,QTE_q0,file="QTE_allp.Rdata")


system.time(QTE_q1<-parApply(cl,ID,1,QTE_boot,bweight2,prob=0.5,ToC=1))
system.time(QTE_q0<-parApply(cl,ID,1,QTE_boot,bweight2,prob=0.5,ToC=0))
QTE<-QTE_q1-QTE_q0
save(QTE,QTE_q1,QTE_q0,file="/u5/y63xie/MacProfile/project3/data_application/QTE_0.5.Rdata")




