# tau: quantil, when set 0.5 is the median
# link: "logit", "probit" or "cloglog"

sample.mkreg <- function(n,exvar.beta=NA,exvar.nu=NA,exvar.rho=NA,beta=c(0.0),nu=c(0.0),rho=c(0.0),alpha=10,tau=0.5,link="logit")
{
  
  
  #n <- length(y) 
  
  Z<-rep(1,n)
  c=0
  exvar.beta=as.matrix(exvar.beta)
  k=if(is.matrix(exvar.beta)){ncol(exvar.beta)}else{1}
  X <- matrix(c(rep(1,n),exvar.beta), nrow=n, ncol=(k+1), byrow=F)
  
  
    ##funções de ligação
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  
  #gerar amostra modified kumaraswamy reparametrizada
  rimk<-function(n,mu,alpha,lambda)
  {
    u<- runif(n,0,1)
    mu[is.na(mu)]<-.Machine$double.eps
    mu[mu<.Machine$double.eps]<-.Machine$double.eps
    mu[mu>0.9999999]<-0.9999999
    critical<-(alpha-alpha/mu)
    critical[is.na(critical)]<- -.Machine$double.eps
    critical[is.nan(critical)]<- -36.04365
    critical[critical < (-36.04365)]<- -36.04365
    den.cr=log(1-exp(critical))
    den.cr[is.nan(den.cr)]<--36.04365
    beta0_cond<-log(1-tau)/den.cr
    y<-ifelse(u<=lambda,0,alpha/(alpha-log(1-(1-(u-lambda)/(1-lambda))^(1/beta0_cond))))    
    if(any(y <0)| any(y >1)| any(is.nan(y))){
      while(any(y <0)| any(y >1)| any(is.nan(y))){
      u<- runif(n,0,1)
      y<-ifelse(u<=lambda,0,alpha/(alpha-log(1-(1-(u-lambda)/(1-lambda))^(1/beta0_cond))))    
      }
    }
    return(y)
  }
  

  lambda<-linkinv(Z%*%as.matrix(nu))
  mu <- linkinv(X%*%as.matrix(beta))
  y <- rimk(n, mu, alpha,lambda)
  #print(y)
  return(y)
}
