envelope.mkreg<-function(residual,n,beta,exvar,alpha,tau,link="logit",conf=0.95)
{
  source("sample.mkreg.R")
  source("mkreg-env.R")
  smk<-sr<-matrix(NA,nrow=n,ncol=100)
  for (i in 1:100){
    smk[,i]<-sample.mkreg(n,beta=beta,exvar=exvar,alpha=alpha,tau=tau,link=link)
    mkr<-mkreg.env(smk[,i],exvar=exvar,tau=tau,resid=1,graph=F,check=F,print=F,link=link)
    if(mkr$RMC==0){
      sr[,i]<-sort(mkr$residual)}else{
        while(mkr$RMC==1){smk[,i]<-sample.mkreg(n,beta=beta,exvar=exvar,alpha=alpha,tau=tau,link=link)
        mkr<-mkreg.env(smk[,i],exvar=exvar,tau=tau,resid=1,graph=F,check=F,print=F,link=link)}
        sr[,i]<-sort(mkr$residual)
      }
  }
  m <- apply(sr,1,mean)
  qlb<-function(x)
  {
    quantile(x,(1-conf)/2)
  }
  qub<-function(x)
  {
    quantile(x,1-(1-conf)/2)
  }
  
  lb<-apply(sr,1,qlb)
  ub<-apply(sr,1,qub)
  
  tb <- range(residual,lb,ub)
  par(pty="s")
  qqnorm(residual,main="",xlab="Normal quantile",ylab="Empirical quantile", ylim=tb, pch=16)
  par(new=T)
  qqnorm(lb,axes=F,main="",xlab="",ylab="",type="l",ylim=tb,lty=1)
  par(new=T)
  qqnorm(ub,axes=F,main="",xlab="",ylab="", type="l",ylim=tb,lty=1)
  par(new=T)
  qqnorm(m,axes=F,main="",xlab="", ylab="", type="l",ylim=tb,lty=2)
}