pe.dr <- function (data,prob,normalize=FALSE, links="logit") {
  data<-as.data.frame(data)
  t<-data$MomSmoke
  y1<-data$Weight[data$MomSmoke==1];y0<-data$Weight[data$MomSmoke==0]
  y<-data$Weight
  n <- length(data$Weight)
  X.ps <-model.matrix(~Black+ Married+ Boy+ MomAge+ as.factor(MomEdLevel)+
                        + I(MomAge^2)+ Black:MomAge+ Married:MomAge+ Boy:MomAge+ 
                        as.factor(MomEdLevel):MomAge,data=data)

  X.or <-X.ps
  
  MX0 <- X.or[t==0,]
  MX1 <- X.or[t==1,]
  #OR model
  fm0  <- lm(y0~0+MX0)
  bet0 <- fm0$coefficients
  sig0 <- sqrt(sum(fm0$residuals^2)/fm0$df.residual)
  fm1  <- lm(y1~0+MX1)
  bet1 <- fm1$coefficients
  sig1 <- sqrt(sum(fm1$residuals^2)/fm1$df.residual)
  
  W.bet0 <- X.or%*%bet0
  W.bet1 <- X.or%*%bet1
  #ps model
  psm1 <- glm(t~X.ps-1,family=binomial(link = links))  
  t.hat <- psm1$fitted.values
  w0 <- ifelse(t==0,1/(1-t.hat),0)
  w1 <- ifelse(t==1,1/(t.hat),0)
  
  if (normalize) {
    w0 <- w0/sum(w0)
    w1 <- w1/sum(w1)
  }
  f0 <- function(q) {
    G <- pnorm((q-W.bet0)/sig0)
    I <- as.numeric(y<=q)
    abs(mean((G-prob)+w0*(I-G)))
  }
  f1 <- function(q) {
    G <- pnorm((q-W.bet1)/sig1)
    I <- as.numeric(y<=q)
    abs(mean((G-prob)+w1*(I-G)))
  }
  q0 <- quantile(y[t==0],probs=prob)
  q0 <- optim(q0,f0,lower=min(y[t==0]),upper=max(y[t==0]),method="Brent")$par
  q1 <- quantile(y[t==1],probs=prob)
  q1 <- optim(q1,f1,lower=min(y[t==1]),upper=max(y[t==1]),method="Brent")$par
  c(q0,q1,q1-q0)
}
prob=0.5
q_logit = pe.dr(bweight2,prob,normalize=FALSE, links="logit") 
q_probit = pe.dr(bweight2,prob,normalize=FALSE, links="probit") 
q_cloglog = pe.dr(bweight2,prob,normalize=FALSE, links="cloglog") 
q_logit
q_probit 
q_cloglog


 QTE_boot2<-function(id,x,prob=0.5,links){
   x<-x[id,]
   rownames(x)=NULL
   out<-pe.dr(data=x,prob=prob,links=links)
   return(out)
 }

clusterExport(cl,c("pe.dr","optim","glm","QTE_boot2"))
p=0.5
system.time(DR1<-parApply(cl,ID,1,QTE_boot2,bweight2,prob=p,links="logit"))
system.time(DR2<-parApply(cl,ID,1,QTE_boot2,bweight2,prob=p,links="probit"))
system.time(DR3<-parApply(cl,ID,1,QTE_boot2,bweight2,prob=p,links="cloglog"))
save(DR1,DR2,DR3,file="DR_0.5.Rdata")

p = seq(0.05,0.95,by=0.05)
dr1<-numeric(length(p))
dr2<-dr1
dr3<-dr1
for (i in 1:length(p)){
  q1 = pe.dr(bweight2,p[i],normalize=FALSE, links="logit") [3]
  q2 = pe.dr(bweight2,p[i],normalize=FALSE, links="probit") [3]
  q3 = pe.dr(bweight2,p[i],normalize=FALSE, links="cloglog") [3]
  dr1[i] = q1
  dr2[i] = q2
  dr3[i] = q3
  
} 


load("QTE_allp.Rdata")
load("other_quantiles.Rdata")

op <- par(cex = 0.8)
y.low <- qte - 1.96*apply(q,1,sd)
y.high <- qte + 1.96*apply(q,1,sd)

p = seq(0.05,0.95,by=0.05)

plot(p,qte,type = 'n', ylim = c(-350, -100),
     xlab = 'Probability', ylab = 'Quantile Treatment Effect')
lines(p, y.low, col = 'grey')
lines(p, y.high, col = 'grey')

polygon(c(p, rev(p)), c(y.high, rev(y.low)),
        col = rgb(211, 211, 211, maxColorValue=255), border = NA)


lines(p,qte,lwd=2)

lines(x=p,y=dr1,col='red',lty=2,lwd=2)
lines(x=p,y=dr2,col='blue',lty=3,lwd=2)
lines(x=p,y=dr3,col='green',lty=4,lwd=2)
legend(0.4,-250,legend=c("MR","DR_logit",'DR_probit',"DR_cloglog"),col=c("black","red","blue","green"),
        lty=c(1,2,3,4),bty = 'n',y.intersp = 0.6)
