#generate data without correlation (combind VS no combind)
#generate data with new correlation (combind VS no combind)


Datasimulate<-function(n,setting="mixed",indep=FALSE,outm,psm){
  if(indep==TRUE){
    if(setting=="continuous"){
      w1<-rnorm(n,0,1);w2<-rnorm(n,0,1);w3<-rnorm(n,0,1);w4<-rnorm(n,0,1)
      w5<-rnorm(n,0,1);w6<-rnorm(n,0,1);w7<-rnorm(n,0,1);w8<-rnorm(n,0,1)
      w9<-rnorm(n,0,1)
    }else{
      ind1<-rbinom(n,1,0.5)
      ind3<-rbinom(n,1,0.5)
      ind4<-rbinom(n,1,0.5)
      ind5<-rbinom(n,1,0.5)
      ind8<-rbinom(n,1,0.5)
      w1<-ind1*rnorm(n,1,1)+(1-ind1)*rnorm(n,-1,1)
      w3<-ind3*rnorm(n,1,1)+(1-ind3)*rnorm(n,-1,1)
      w4<-ind4*rnorm(n,1,1)+(1-ind4)*rnorm(n,-1,1)
      w5<-ind5*rnorm(n,1,1)+(1-ind5)*rnorm(n,-1,1)
      w8<-ind8*rnorm(n,1,1)+(1-ind8)*rnorm(n,-1,1)
      
      w2<-rbinom(n,1,0.5);w6<-rbinom(n,1,0.5);w7<-rbinom(n,1,0.5)
      w9<-rbinom(n,1,0.5)
      w1<-rnorm(n,0,1);w3<-rnorm(n,0,1);w4<-rnorm(n,0,1);w5<-rnorm(n,0,1)
      w8<-rnorm(n,0,1);w2<-rbinom(n,1,0.5);w6<-rbinom(n,1,0.5);w7<-rbinom(n,1,0.5)
      w9<-rbinom(n,1,0.5)
    }
  }else{
    w1=rnorm(n, mean=0, sd=1)
    w2=rnorm(n, mean=0, sd=1)
    w4=F.sample.cor(w1, -0.33)
    w5=F.sample.cor(w2, -0.5)
    
    w3=rnorm(n, mean=0, sd=1)
    w6=F.sample.cor(w3, -0.4)
    
    w7=rnorm(n, mean=0, sd=1)
    w8=F.sample.cor(w7, 0.5)
    w9=rnorm(n, mean=0, sd=2)
    if(setting=="mixed"){
      w2=ifelse(w2>mean(w2),1,0)  # dichotomizing covariates
      w6=ifelse(w6>mean(w6),1,0)
      w7=ifelse(w7>mean(w7),1,0)
      w9=ifelse(w9>mean(w9),1,0)
    }
  }
  
  if(psm==1){
    fx<-log(2)*w1+log(1.4)*w2+log(2)*w4+log(1.4)*w5+log(2)*w7+log(1.4)*w8
    
  }else if(psm==2){
    x1<- -exp(w1/3);x4<-w1-w4+3;x5<-w1*w5/10+0.5;x8<--exp(w8/2)
    fx<-log(2)*x1+log(1.4)*w2+log(2)*x4+log(1.4)*x5+log(2)*w7+log(1.4)*x8
    
  }else if(psm==3){
    x1<-exp(w1);x4<- -2*w4/(1+exp(w1))+0.5;x5<- -w4*w5/2-2;x8<-abs(w8)-1
    fx<-log(2)*x1+log(1.4)*w2^2+log(2)*x4+log(1.4)*x5+log(2)*w7+log(1.4)*x8
    
    }
  
  ps<-1/(1+exp(fx))
  z<-rbinom(length(ps),1,ps)
  #hist(ps)
  #hist(z)
  #sum(z==1)
  if(outm=="A"){
    y1<- -2.4+1.68*w1+1.68*w2+1.68*w3+3.47*w4+3.47*w5+3.47*w6+5*w1+rnorm(n)
    y0<- -2.4+1.68*w1+1.68*w2+1.68*w3+3.47*w4+3.47*w5+3.47*w6+rnorm(n)
    
  }else if(outm=="B"){
    y1<--2.4+1.68*w1+1.68*exp(w2/2)+1.68*w3+3.47*exp(w4/3)+3.47*w5*w6+3.47*w6+1.68*w1^2+5*w1+rnorm(n)
    y0<--2.4+1.68*w1+1.68*exp(w2/2)+1.68*w3+3.47*exp(w4/3)+3.47*w5*w6+3.47*w6+1.68*w1^2+rnorm(n)
  }else if(outm=="C"){
    y1<--2.4+1.68*w1+1.68*w2^2+1.68*exp(w3/3)+3.47*abs(w4)+3.47*w5+3.47*w6^2+1.68*w1^3/3+5*w1+rnorm(n)
    y0<--2.4+1.68*w1+1.68*w2^2+1.68*exp(w3/3)+3.47*abs(w4)+3.47*w5+3.47*w6^2+1.68*w1^3/3+rnorm(n)
    }
  
  y<-y1*z+y0*(1-z)
  data<-data.frame(w1,w2,w3,w4,w5,w6,w7,w8,w9,z,y1,y0,y)
  return(data)
  
}

F.sample.cor=function(x, rho) {
  y <- (rho * (x - mean(x)))/sqrt(var(x)) + sqrt(1 - rho^2) * rnorm(length(x))
  #cat("Sample corr = ", cor(x, y), "\n")
  return(y)
}


################################################################################################################################
#generate data with new correlation (combind VS no combind)
set.seed(1234)
simdata1A<-replicate(1000,Datasimulate(200,setting="continuous",indep=TRUE,psm=1,outm="A"))
simdata1B<-replicate(1000,Datasimulate(1000,setting="continuous",indep=TRUE,psm=1,outm="B"))
simdata2A<-replicate(1000,Datasimulate(1000,setting="continuous",indep=TRUE,psm=2,outm="A"))
simdata2B<-replicate(1000,Datasimulate(1000,setting="continuous",indep=TRUE,psm=2,outm="B"))
save(simdata1A,simdata1B,simdata2A,simdata2B,file="/u5/y63xie/MacProfile/project3/datasim5/simdata1000ind.Rdata")
save(simdata1A,simdata1B,simdata2A,simdata2B,file="N:/MacProfile/project3/datasim5/simdata1000ind.Rdata")

set.seed(1234)
simdata3C<-replicate(1000,Datasimulate(200,setting="continuous",indep=TRUE,psm=3,outm="C"))
save(simdata3C,file="/u5/y63xie/MacProfile/project3/datasim5/200/simdata200_3C.Rdata")


data<-Datasimulate(1000000,setting="continuous",indep=TRUE,psm=3,outm="C")
#estimated ATE bias
mean(data$y[data$z==1])-mean(data$y[data$z==0])
#true ATE bias
mean(data$y1)-mean(data$y0)
#ture qte
quantile(data$y1,0.25)-quantile(data$y0,0.25)
quantile(data$y[data$z==1],0.25)-quantile(data$y[data$z==0],0.25)

quantile(data$y1,0.5)-quantile(data$y0,0.5)
quantile(data$y[data$z==1],0.5)-quantile(data$y[data$z==0],0.5)

quantile(data$y1,0.75)-quantile(data$y0,0.75)
quantile(data$y[data$z==1],0.75)-quantile(data$y[data$z==0],0.75)

quantile(data$y1,0.95)-quantile(data$y0,0.95)
quantile(data$y[data$z==1],0.95)-quantile(data$y[data$z==0],0.95)


quantile(data$y1,0.001)-quantile(data$y0,0.001)


y<-data$y
z<-data$z
a<-y[z==1]
b<-y[z==0]
hist(a,breaks=50)
hist(b,breaks=50,add=TRUE,col="red")
fa<-ecdf(a)
fb<-ecdf(b)
ks.test(a,b)
length(a)
length(b)
plot(fa)
plot(fb,add=TRUE)


###TRUE ecdf
a<-data$y1
b<-data$y0
quantile(a1,0.5)-quantile(b1,0.5)
hist(a,breaks=50)
hist(b,breaks=50,add=TRUE,col="red")
fa<-ecdf(a)
fb<-ecdf(b)
ks.test(a,b)
length(a)
length(b)
plot(fa,col="red")
plot(fb,add=TRUE)






set.seed(123)
true_quantile<-function(n,psm,outm){
  datax<-Datasimulate(n,setting="continuous",indep=TRUE,psm=psm,outm=outm)
  #datax<-Datasimulate(1000000,setting="mixed",indep=TRUE,psm=2,outm="B")
  y00<-datax$y0
  y11<-datax$y1
  q1_true1<-quantile(y11,0.25)
  q1_true2<-quantile(y11,0.5)
  q1_true3<-quantile(y11,0.75)
  q1_true4<-quantile(y11,0.95)
  #3.218089
  q0_true1<-quantile(y00,0.25)
  q0_true2<-quantile(y00,0.5)
  q0_true3<-quantile(y00,0.75)
  q0_true4<-quantile(y00,0.95)
  c(q1_true1,q1_true2,q1_true3,q1_true4,q0_true1,q0_true2,q0_true3,q0_true4)
}

library(parallel)
cl<-makeCluster(detectCores())
clusterExport(cl,c("Datasimulate","true_quantile","quantile","F.sample.cor"))
n<-1000000
num<-rep(n,10000)
set.seed(1)
system.time(Q4<-parSapply(cl,num,true_quantile,psm=2,outm="B"))
apply(Q4,1,mean)
apply(Q4[1:4,],1,mean)-apply(Q4[5:8,],1,mean)
system.time(Q3<-parSapply(cl,num,true_quantile,psm=1,outm="B"))
apply(Q3,1,mean)
apply(Q3[1:4,],1,mean)-apply(Q3[5:8,],1,mean)
system.time(Q1<-parSapply(cl,num,true_quantile,psm=1,outm="A"))
apply(Q1,1,mean)
apply(Q1[1:4,],1,mean)-apply(Q1[5:8,],1,mean)
system.time(Q2<-parSapply(cl,num,true_quantile,psm=2,outm="A"))
apply(Q2,1,mean)
apply(Q2[1:4,],1,mean)-apply(Q2[5:8,],1,mean)

save(Q1,Q2,Q3,Q4,file="/u5/y63xie/MacProfile/project3/datasim3/true_q.Rdata")

set.seed(1)
system.time(Q5<-parSapply(cl,num,true_quantile,psm=3,outm="C"))
q5<-apply(Q5[1:4,],1,mean)-apply(Q5[5:8,],1,mean)
save(Q5,q5,file="/u5/y63xie/MacProfile/project3/datasim5/true_q5.Rdata")
