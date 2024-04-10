# tau: quantil, when set 0.5 is the median
# link: "logit", "probit" or "cloglog"

sample.mkreg <- function(n,beta=c(0.0),exvar=NA,alpha=10,tau=0.5,link="logit")
{
  
  exvar=as.matrix(exvar)
  k=if(is.matrix(exvar)){ncol(exvar)}else{1}
  X <- matrix(c(rep(1,n),exvar), nrow=n, ncol=(k+1), byrow=F)
  
  
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
  rmkum_a<-function(n,mu,alpha)
  {
    u<- runif(n,min=.Machine$double.eps)
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
    y<- alpha/(alpha-log(1-(1-u)^(1/beta0_cond)))
    if(any(y <= .Machine$double.eps)| any(y >= 0.9999999)){
      while(any(y <= .Machine$double.eps)| any(y >= 0.9999999)){
      u<- runif(n,min=.Machine$double.eps)
      y<- alpha/(alpha-log(1-(1-u)^(1/beta0_cond)))
      }
    }
    return(y)
  }
  
  mu <- linkinv(X%*%as.matrix(beta))
  y    <- rmkum_a(n, mu, alpha)

  return(y)
}
