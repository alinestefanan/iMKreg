# y: time series
# S: seasonal periodicity
# exvar: covariate column matrix
# tau: quantil, when set 0.5 is the median
# resid: 1 = quantile residual, 2 = deviance residual
# link: "logit", "probit" or "cloglog"

mkreg.env <- function(y,exvar.beta=NA,tau=0.5,resid=1,graph=T,print=T,check=F,link="logit")
{
  n <- length(y) 
  exvar=as.matrix(exvar.beta)
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
  if (any(linktemp == "logit"))
  {  
    stats <- make.link(linktemp)
    dif2link = function(muhat) eval(D(expression(-1/((muhat-1)*muhat)),"muhat"),list(muhat=muhat))
  }else if (any(linktemp == "probit"))
  {  
    stats <- make.link(linktemp)
    library(VGAM)
    dif2link = function(muhat) probitlink(muhat,  deriv = 2)#segunda derivada
  }else if (any(linktemp == "cloglog"))
  {  
    stats <- make.link(linktemp)
    dif2link = function(muhat) eval(D(expression(1/(log(1-muhat)*(muhat-1))),"muhat"),list(muhat=muhat))#segunda derivada
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv,
                         mu.eta = stats$mu.eta,#derivada de mu em relação a eta
                         diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- link$diflink
 
  #função densidade de probabilidade modified kumaraswamy
  dmk <- Vectorize(function(y,alpha,beta,log = FALSE){
    critical.y=exp(alpha-alpha/y)
    critical.y[is.infinite(critical.y)]<- .Machine$double.xmax
    density<-(alpha*beta*critical.y*(1-critical.y)^(beta-1))/(y^2)
    density[is.na(density)] <- .Machine$double.eps
    density[is.nan(density)] <- .Machine$double.eps
    density[density<.Machine$double.eps] <- .Machine$double.eps
    density[density>.Machine$double.xmax] <- .Machine$double.xmax
    logden <- log(density)
    val <- ifelse(log, logden, exp(logden)) 
    return(val)
  }) #função Murilo
  
  #modified kumaraswamy cumulative distribution function 
  pmk <- Vectorize(function(q,alpha,beta,log.p = FALSE){
    cdf <-  1-(1-exp(alpha-alpha/q))^beta
    cdf[is.na(cdf)]<- .Machine$double.eps
    cdf[cdf<.Machine$double.eps]<-.Machine$double.eps
    cdf[cdf>0.9999999]<-0.9999999
    val <- ifelse(log.p, log(cdf), cdf)
    return(val)
  })#função Murilo
  
  dmk_alpha<-function(alpha){
    beta=log(1-tau)/log(1-exp(alpha-alpha/median(y)))
    if (is.na(beta)){beta=.Machine$double.eps}
    if(beta<.Machine$double.eps) beta<-.Machine$double.eps
    if (is.infinite(beta)){beta=.Machine$double.xmax}
    critical.y=exp(alpha-alpha/y)
    critical.y[is.infinite(critical.y)]<-.Machine$double.xmax
    density<-(alpha*beta*critical.y*(1-critical.y)^(beta-1))/(y^2)
    density[is.na(density)] <- .Machine$double.eps
    density[is.nan(density)] <- .Machine$double.eps
    density[density<.Machine$double.eps]<-.Machine$double.eps
    return(density)
  }
  
  loglik <- function(z) 
  {
    beta <- z[1:(k+1)]
    alpha <- z[(k+2)] # mk parameter      
    mu <- linkinv(X%*%as.matrix(beta))
    mu[is.na(mu)]<-.Machine$double.eps
    mu[mu<.Machine$double.eps]<-.Machine$double.eps
    mu[mu>0.9999999]<-0.9999999
    critical<-alpha-alpha/mu
    critical[is.na(critical)]<--.Machine$double.eps
    critical[is.nan(critical)]<--36.04365
    critical[critical< (-36.04365)]<--36.04365
    den.cr=log(1-exp(critical))
    den.cr[is.nan(den.cr)]<--36.04365
    critical.y<-exp(alpha-alpha/y)
    critical.y[is.infinite(critical.y)]<-.Machine$double.xmax
    critical.ly<-log(1-critical.y)
    critical.ly[is.nan(critical.ly)]<--36.04365
    critical.ly[critical.ly< (-36.04365)]<--36.04365#para exp dar .Machine$double.eps
    # l=log(alpha)+alpha-alpha/y+(log(1-tau)/den.cr -1)*critical.ly-2*log(y)+log(log(1-tau)/den.cr)
    # print("l");print(sum(l))
    ll <- dmk(y, alpha, log(1-tau)/den.cr, log = TRUE)#log-density modified kumaraswamy quantile re-parametrization
    return(sum(ll))
    # print("ll");print(sum(ll))
    # return(sum(l))
  }#fim loglik
  
  score <- function(z) 
  {
    beta <- z[1:(k+1)]
    alpha <- z[(k+2)] # mk parameter      
    mu <- linkinv(X%*%as.matrix(beta))
    # print("mu");print(mu)
    mu.eta <-  link$mu.eta
    # mu[is.na(mu)]<-.Machine$double.eps
    # mu[mu<.Machine$double.eps]<-.Machine$double.eps
    # mu[mu>0.9999999]<-0.9999999
    critical<-alpha-alpha/mu
    # print("critical");print(critical)
    # critical[is.na(critical)]<--.Machine$double.eps
    # critical[is.nan(critical)]<--36.04365
    # critical[critical< (-36.04365)]<--36.04365
    den.cr<-log(1-exp(critical))
    # print("den.cr");print(den.cr)
    den.cr[is.nan(den.cr)]<--36.04365
    critical.y<-exp(alpha-alpha/y)
    # print("critical.y");print(critical.y)
    # critical.y[is.infinite(critical.y)]<-.Machine$double.xmax
    critical.ly<-log(1-critical.y)
    # print("critical.ly");print(critical.ly)
    critical.ly[is.nan(critical.ly)]<--36.04365
    # critical.ly[critical.ly<(-36.04365)]<--36.04365
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU
    ###########################################################################################################
    num1<-exp(alpha)*alpha*(den.cr+log(1-tau)*critical.ly)
    den1<-(mu^2)*(exp(alpha/mu)-exp(alpha))*((den.cr)^2)
    mustar<-num1/den1
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU      
   
    mT <- diag(as.vector(mu.eta(X%*%as.matrix(beta))))
    Ubeta <-  t(X) %*% mT %*% as.matrix(mustar)
    
    ##############################################################################
    ##### START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO alpha
    num1<-(mu-1)*exp(alpha)*log(1-tau)*critical.ly
    den1<- mu*(exp(alpha/mu)-exp(alpha))*(den.cr)^2
    num2.part1<-(y-1)*critical.y
    num2.part2<-(log(1-tau)/den.cr-1)
    den2<-y-y*critical.y
    den3<- mu*(exp(alpha/mu)-exp(alpha))*den.cr
    Ualpha<-sum(1+1/alpha+num1/den1-num2.part1*num2.part2/den2+((mu-1)*exp(alpha))/den3-1/y)
    ##### END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO alpha
    ##############################################################################
  
    return(c(Ubeta,Ualpha))
  }#end score
  
  # initial values for estimation
  Ynew = linkfun(y) 
  ajuste = lm.fit(X, Ynew) 
  mqo = c(ajuste$coef)
  mqo[is.na(mqo)]<-0
  
  library(GenSA)
  on.dmk.alpha<-function(alpha){-sum(log(dmk_alpha(alpha)))}
  
  gen.semchute <- GenSA(lower = c(.Machine$double.eps), 
                        upper = c(100),
                        fn = on.dmk.alpha, control=list(max.time=2))
  alpha<-gen.semchute$par
  reg <- c(mqo, alpha) # initializing the parameter values
  #reg=c(0,0,alpha)
 # reg <- c(mqo,0)
  z <- c()
  opt.error<- tryCatch(optim(reg, loglik, score,
                             method = "BFGS",
                             control = list(fnscale = -1)), error = function(e) return("error"))
  if(opt.error[1] == "error")
  {z$RMC=1
  warning("optim error")
  return(z)
  }
  opt <- optim(reg, loglik, score, 
               method = "BFGS", hessian = T, 
               control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
  if (opt$conv != 0)
  {
    warning("FUNCTION DID NOT CONVERGE!")
    z$RMC=1
    return(z)
  }
 
  z$conv <- opt$conv
  coef <- (opt$par)[1:(k+2)]
  z$coeff <- coef
  z$loglik <- opt$value
  beta <-coef[1:(k+1)]
  alpha <-coef[(k+2)] # mk parameter
  
  z$beta <- beta
  z$alpha <- alpha
  
  z$RMC=0
   
  muhat <- linkinv(X%*%as.matrix(beta))
  
  z$fitted <- ts(muhat,start=start(y),frequency=frequency(y))
  
  z$serie <- y
  #quantile residuals
  
  
  z$fitted[is.na(z$fitted)]<-.Machine$double.eps
  critical<-(alpha-alpha/z$fitted)
  # critical[is.na(critical)]<--.Machine$double.eps
  # critical[is.nan(critical)]<--36.04365
  # critical[critical< (-36.04365)]<--36.04365
  den.cr=log(1-exp(critical))
  den.cr[is.nan(den.cr)]<--36.04365
  z$resid1 <- as.vector(qnorm(pmk(y[1:n],alpha, log(1-tau)/den.cr[1:n],log.p = FALSE ) ))
  #deviance residuals
  l_tilde <- (dmk(y[1:n], alpha, log(1-tau)/log(1-exp(alpha-alpha/y[1:n])), log = TRUE))#y[1:n] where was mu
  l_hat <- (dmk(y[1:n], alpha, log(1-tau)/den.cr[1:n], log = TRUE))#z$fitted[1:n] where was mu
  #print(sum(l_hat));print(opt$value)
  dt <- (l_tilde-l_hat)
  dt[which(dt<0)]<-0
  
  r2a<-sign(y[1:n]-z$fitted[1:n])
  r2b<-sqrt(2*(dt))
  z$resid2<-r2a*r2b
  
  z$deviance <- 2*sum(dt)
  z$dof.dev=n-k#desconsidera intercepto do eta e alpha da distribuição
  z$p_deviance <- 1 - pchisq(z$deviance, z$dof.dev)
 
  mresult<-matrix(round(c(z$loglik,z$aic,z$bic,z$deviance),4),nrow=4,ncol=1)
  rownames(mresult)<-c("Log-likelihood","AIC","BIC","Deviance")
  colnames(mresult)<-c("")
  z$mresult<-mresult


  if(resid==1) {
    residual <- z$resid1
  }
  
  if(resid==2) {
    residual <- z$resid2
  }
  
  z$residual<-residual
 
  return(z)
}#fim estimação
