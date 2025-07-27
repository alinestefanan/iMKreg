# y: time series
# S: seasonal periodicity
# exvar.beta,exvar.nu: covariate column matrix
# tau: quantil, when set 0.5 is the median
# link: "logit", "probit" or "cloglog"

imkreg0Z1 <- function(y,exvar.beta=NA,exvar.nu=NA,exvar.rho=NA,tau=0.5,graph=T,print=T,check=F,link="logit")
{
  n <- length(y) 
  
  exvar.nu=as.matrix(exvar.nu)
  c=if(is.matrix(exvar.nu)){ncol(exvar.nu)}else{1}
  Z<-matrix(c(rep(1,n),exvar.nu), nrow=n, ncol=(c+1), byrow=F)
  
  m=0
  A<-matrix(rep(1,n), nrow=n, ncol=1, byrow=F)
 
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
  
  #função densidade de probabilidade inflated modified kumaraswamy
  dmk <- Vectorize(function(y,alpha,beta,lambda,p,log = FALSE){
    critical.y=exp(alpha-alpha/y)
    critical.y[is.infinite(critical.y)]<- .Machine$double.xmax
    density<-ifelse(y!=0 & y!=1,(1-lambda)*(alpha*beta*critical.y*(1-critical.y)^(beta-1))/(y^2),NA)
    density<-ifelse(y==0,lambda*(1-p),density)
    density<-ifelse(y==1,lambda*p,density)
    density[is.na(density)] <- .Machine$double.eps
    density[is.nan(density)] <- .Machine$double.eps
    density[density<.Machine$double.eps] <- .Machine$double.eps
    density[density>.Machine$double.xmax] <- .Machine$double.xmax
    logden <- log(density)
    val <- ifelse(log, logden, exp(logden)) 
    return(val)
  }) 
  
  #modified kumaraswamy cumulative inflated distribution function 
  pmk <- Vectorize(function(q,alpha,beta,lambda,p,log.p = FALSE){
    cdf <-  ifelse(y!=0 & y!=1,lambda*(1-p)+(1-lambda)*(1-(1-exp(alpha-alpha/q))^beta),NA)
    cdf <-  ifelse(y==0,lambda*(1-p),cdf)
    cdf <-  ifelse(y==1,lambda*(1-p)+lambda*p+(1-lambda)*(1-(1-exp(alpha-alpha/q))^beta),cdf)
    cdf[is.na(cdf)]<- .Machine$double.eps
    cdf[cdf<.Machine$double.eps]<-.Machine$double.eps
    cdf[cdf>0.9999999]<-0.9999999
    val <- ifelse(log.p, log(cdf), cdf)
    return(val)
  })
  
  dmk_alpha<-function(alpha){
    y1<-y[which(y!=0 & y!=1)]
    beta=log(1-tau)/log(1-exp(alpha-alpha/median(y1)))
    if (is.na(beta)){beta=.Machine$double.eps}
    if(beta<.Machine$double.eps) beta<-.Machine$double.eps
    if (is.infinite(beta)){beta=.Machine$double.xmax}
    lambda=length(y[which(y==0 | y==1)])/length(y)
    p=length(y[which(y==1)])/length(y)
    critical.y=exp(alpha-alpha/y)
    critical.y[is.infinite(critical.y)]<-.Machine$double.xmax
    density<-ifelse(y!=0 & y!=1,(1-lambda)*(alpha*beta*critical.y*(1-critical.y)^(beta-1))/(y^2),NA)
    density<-ifelse(y==0,lambda*(1-p),density)
    density<-ifelse(y==1,lambda*p,density)
    density[is.na(density)] <- .Machine$double.eps
    density[is.nan(density)] <- .Machine$double.eps
    density[density<.Machine$double.eps]<-.Machine$double.eps
    return(density)
  }
  
  # função quantílica inflated modified kumaraswamy reparametrizada
  rimk<-function(u,mu,alpha,lambda,p)
  {
    mu[is.na(mu)] <- .Machine$double.eps
    mu[mu<.Machine$double.eps] <- .Machine$double.eps
    mu[mu>0.9999999]<-0.9999999
    critical<-(alpha-alpha/mu)
    critical[is.na(critical)]<--.Machine$double.eps
    critical[is.nan(critical)]<--36.04365
    critical[critical< (-36.04365)]<--36.04365
    den.cr=log(1-exp(critical))
    den.cr[is.nan(den.cr)]<--36.04365
    beta0_cond<-log(1-tau)/den.cr   
    y<-ifelse(u<=lambda*(1-p),0,NA)
    y<-ifelse(u>=(1-lambda*p),1,y)
    y<-ifelse(u>lambda*(1-p) & u<(1-lambda*p),alpha/(alpha-log(1-(1-(u-lambda*(1-p))/(1-lambda))^(1/beta0_cond))),y)
    
    return(y)
  }
  
  loglik <- function(z) 
  {
    beta <- z[1:(k+1)]
    alpha <- z[(k+2)] #mk parameter
    nu<-z[(k+3):(k+3+c)]
    rho<-z[(k+4+c):(k+4+c+m)]
    lambda<-linkinv(Z%*%as.matrix(nu))
    p<-linkinv(A%*%as.matrix(rho))
    mu <- linkinv(X%*%as.matrix(beta))
    mu[is.na(mu)] <- .Machine$double.eps
    mu[mu<.Machine$double.eps] <- .Machine$double.eps
    mu[mu>0.9999999]<-0.9999999
    critical<-(alpha-alpha/mu)
    critical[is.na(critical)]<--.Machine$double.eps
    critical[is.nan(critical)]<--36.04365
    critical[critical< (-36.04365)]<--36.04365
    den.cr=log(1-exp(critical))
    den.cr[is.nan(den.cr)]<--36.04365
    critical.y=exp(alpha-alpha/y)
    critical.y[is.infinite(critical.y)]<-.Machine$double.xmax
    critical.ly<-log(1-critical.y)
    critical.ly[is.nan(critical.ly)]<--36.04365
    critical.ly[critical.ly< (-36.04365)]<--36.04365
    l=ifelse(y==0,log(lambda)+log(1-p),NA)
    l=ifelse(y==1,log(lambda)+log(p),l)
    l=ifelse(y!=0 & y!=1,log(1-lambda)+log(alpha)+alpha-alpha/y+(log(1-tau)/den.cr -1)*critical.ly-2*log(y)+log(log(1-tau)/den.cr),l)
    # print("l");print(l)
    # print("log(alpha)+alpha-alpha/y");print(log(alpha)+alpha-alpha/y)
    # print("(log(1-tau)/den.cr[i] -1)*critical.ly[i]");print((log(1-tau)/den.cr[i] -1)*critical.ly[i])
    # print("-2*log(y[i])");print(-2*log(y[i]))
    # print('log(log(1-tau)/den.cr[i])');print(log(log(1-tau)/den.cr[i]))
    # print("l");print(sum(l))
    # ll <- dmk(y, alpha, log(1-tau)/den.cr, lambda, p,log = TRUE)#log-density modified kumaraswamy quantile re-parametrization
    # print("ll");print(ll)
    # print("sum(ll)");print(sum(ll))
    return(sum(l))
  }#fim loglik
  
  score <- function(z) 
  {
    beta <- z[1:(k+1)]
    alpha <- z[(k+2)] # mk parameter
    nu<-z[(k+3):(k+3+c)]
    rho<-z[(k+4+c):(k+4+c+m)]
    lambda<-linkinv(Z%*%as.matrix(nu))
    p<-linkinv(A%*%as.matrix(rho))
    mu <- linkinv(X%*%as.matrix(beta))
    mu[is.na(mu)] <- .Machine$double.eps
    mu[mu<.Machine$double.eps] <- .Machine$double.eps
    mu[mu>0.9999999]<-0.9999999
    mu.eta <-  link$mu.eta
    critical<-(alpha-alpha/mu)
    critical[is.na(critical)]<--.Machine$double.eps
    critical[is.nan(critical)]<--36.04365
    critical[critical< (-36.04365)]<--36.04365
    den.cr=log(1-exp(critical))
    den.cr[is.nan(den.cr)]<--36.04365
    critical.y=exp(alpha-alpha/y)
    critical.y[is.infinite(critical.y)]<-.Machine$double.xmax
    critical.ly<-log(1-critical.y)
    critical.ly[is.nan(critical.ly)]<--36.04365
    critical.ly[critical.ly< (-36.04365)]<--36.04365
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU
    ###########################################################################################################
    num1<-exp(alpha)*alpha*(den.cr+log(1-tau)*critical.ly)
    den1<-(mu^2)*(exp(alpha/mu)-exp(alpha))*((den.cr)^2)
    mustar<-ifelse(y==0 | y==1,0,num1/den1)
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU 
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA
    ###########################################################################################################
    lambdastar<-ifelse(y==0,1/lambda,NA)
    lambdastar<-ifelse(y==1,1/lambda,lambdastar)
    lambdastar<-ifelse(y!=0 & y!=1,-1/(1-lambda),lambdastar)
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA 
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO p
    ###########################################################################################################
    pstar<-ifelse(y==1,1/p,NA)
    pstar<-ifelse(y==0,-1/(1-p),pstar)
    pstar<-ifelse(y!=0 & y!=1,0,pstar)
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO p
    
    mT <- diag(as.vector(mu.eta(X%*%as.matrix(beta))))
    lT <- diag(as.vector(mu.eta(Z%*%as.matrix(nu))))
    pT <- diag(as.vector(mu.eta(A%*%as.matrix(rho))))
    Ubeta <-  t(X) %*% mT %*% as.matrix(mustar)
    Unu <-  t(Z) %*% lT %*% as.matrix(lambdastar)
    Urho <-  t(A) %*% pT %*% as.matrix(pstar)
    
    ##############################################################################
    ##### START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO alpha
    num1<-(mu-1)*exp(alpha)*log(1-tau)*critical.ly
    den1<- mu*(exp(alpha/mu)-exp(alpha))*(den.cr)^2
    num2.part1<-(y-1)*critical.y
    num2.part2<-(log(1-tau)/den.cr-1)
    den2<-y-y*critical.y
    den3<- mu*(exp(alpha/mu)-exp(alpha))*den.cr
    Ua<-ifelse(y==0 | y==1,0,1+1/alpha+num1/den1-num2.part1*num2.part2/den2+((mu-1)*exp(alpha))/den3-1/y)
    Ualpha<-sum(Ua)
    ##### END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO alpha
    ##############################################################################
    
    return(c(Ubeta,Ualpha,Unu,Urho))
  }#end score
  
  # initial values for estimation
  Ynew = linkfun(y[y!=0 & y!=1])
  ajuste = lm.fit(X[y!=0 & y!=1,], Ynew)
  
  mqo = c(ajuste$coef)
  mqo[is.na(mqo)]<-0
  
  library(GenSA)
  on.dmk.alpha<-function(alpha){-sum(log(dmk_alpha(alpha)))}
  
  gen.semchute <- GenSA(lower = c(.Machine$double.eps),
                        upper = c(100),
                        fn = on.dmk.alpha, control=list(max.time=2))
  alpha<-gen.semchute$par
  reg <- c(mqo, alpha,length(y[y==0|y==1])/n,rep(0,c),length(y[y==1])/length(y[y==0|y==1]),rep(0,m)) # initializing the parameter values
  #reg <- c(mqo, 0,length(y[y==0])/length(y),rep(0,c)) # initializing the parameter values
  
  # reg <- c(0,rep(0,k), alpha,length(y[which(y==0)])/length(y),rep(0,c)) # initializing the parameter values
  # reg <- c(0,rep(0,k), 0,length(y[which(y==0)])/length(y),rep(0,c)) # initializing the parameter values
  
  #reg=c(0,0,alpha)
  # reg <- c(mqo,0)
  # print(reg)
  z <- c()
  # opt.error<- tryCatch(optim(reg, loglik, score,
  #                            method = "BFGS",
  #                            control = list(fnscale = -1)), error = function(e) return("error"))
  # if(opt.error[1] == "error")
  # {z$RMC=1
  # warning("optim error")
  # return(z)
  # }
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
  coef <- (opt$par)[1:(k+c+4+m)]
  z$coeff <- coef
  z$loglik <- opt$value
  beta <-coef[1:(k+1)] 
  alpha <-coef[(k+2)] # mk parameter
  nu<-coef[(k+3):(k+3+c)]
  rho<-coef[(k+4+c):(k+4+c+m)]
  z$alpha <- alpha
  z$beta <- beta
  z$nu<-nu
  z$rho<-rho
  z$RMC=0
  lambdahat<-linkinv(Z%*%as.matrix(nu))
  phat<-linkinv(A%*%as.matrix(rho))
  # print("z$coeff");print(z$coeff)
  # print("rho");print(rho)
  # print("phat");print(phat)
  muhat <- linkinv(X%*%as.matrix(beta))
  z$fitted<-ts(rimk(u=rep(0.5,n),mu=muhat,alpha=alpha,lambda=lambdahat,p=phat),start=start(y),frequency=frequency(y)) 
  # print(z$fitted)
  # plot(z$fitted,type="l")
  
  obs.inf<-function(y,muhat)
  {
    muhat[is.na(muhat)]<-.Machine$double.eps
    muhat[muhat<.Machine$double.eps]<-.Machine$double.eps
    muhat[muhat>0.9999999]<-0.9999999
    critical<-(alpha-alpha/muhat)
    critical[is.na(critical)]<--.Machine$double.eps
    critical[is.nan(critical)]<--36.04365
    critical[critical< (-36.04365)]<--36.04365
    den.cr=log(1-exp(critical))
    den.cr[is.nan(den.cr)]<--36.04365
    e.ay=exp(alpha/y)
    e.ay[is.infinite(e.ay)]<-.Machine$double.xmax
    e.am=exp(alpha/muhat)
    e.am[is.infinite(e.am)]<-.Machine$double.xmax
    critical.y=exp(alpha-alpha/y)
    critical.y[is.infinite(critical.y)]<-.Machine$double.xmax
    critical.ly<-log(1-critical.y)
    critical.ly[is.nan(critical.ly)]<--36.04365
    critical.ly[critical.ly< (-36.04365)]<--36.04365
    
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU
    ###########################################################################################################
    numerador<- ( exp(alpha)*alpha* ( den.cr*( (2*muhat*exp(alpha)+ e.am*(alpha-2*muhat))*log(1-tau)*critical.ly+exp(alpha)*alpha ) + (2*muhat*exp(alpha) +  e.am*(alpha-2*muhat))*den.cr*den.cr + 2*exp(alpha)*alpha*log(1-tau)*critical.ly) ) 
    denominador <-( (muhat^4)*((exp(alpha)- e.am)^2) *den.cr*den.cr*den.cr)
    muhatstar.sec<-ifelse(y==0 |y==1 ,0,numerador/denominador)
    # print("muhatstar.sec");print(muhatstar.sec)
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECTO TO MU   
    
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA
    ###########################################################################################################
    lambdahatstar.sec<-ifelse(y==0,-1/(lambdahat^2),NA)
    lambdahatstar.sec<-ifelse(y==1,-1/(lambdahat^2),lambdahatstar.sec)
    lambdahatstar.sec<-ifelse(y!=0 & y!=1,-1/(lambdahat-1)^2,lambdahatstar.sec)
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA
    
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO p
    ###########################################################################################################
    phatstar.sec<-ifelse(y==1,-1/(phat^2),NA)
    phatstar.sec<-ifelse(y==0,-1/(phat-1)^2,phatstar.sec)
    phatstar.sec<-ifelse(y!=0 & y!=1,0,phatstar.sec)
    # print("phatstar.sec");print(phatstar.sec)
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO p 
    
    deta.dbetabeta<- array(0,dim=c((k+1),(k+1),n))
    
    for(i in 1:n)
    {
      for(b in 1:(k+1))
      {
        for(a in 1:(k+1))
        {
          deta.dbetabeta[a,b,i] <- 0 
        }
      }
    }
    
    deta.dnunu<- array(0,dim=c((c+1),(c+1),n))
    
    for(i in 1:n)
    {
      for(b in 1:(c+1))
      {
        for(a in 1:(c+1))
        {
          deta.dnunu[a,b,i] <- 0 
        }
      }
    }
    
    deta.drhorho<- array(0,dim=c((m+1),(m+1),n))
    
    for(i in 1:n)
    {
      for(b in 1:(m+1))
      {
        for(a in 1:(m+1))
        {
          deta.drhorho[a,b,i] <- 0 
        }
      }
    }
    
    deta.dnurho<- array(0,dim=c((c+1),(m+1),n))
    
    for(i in 1:n)
    {
      for(b in 1:(m+1))
      {
        for(a in 1:(c+1))
        {
          deta.dnurho[a,b,i] <- 0 
        }
      }
    }
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO ALPHA 
    ###########################################################################################################
    num1<-2*((muhat-1)^2)*exp(2*alpha)*log(1-tau)*critical.ly
    den1<-(muhat^2)*((exp(alpha)-e.am)^2)*((den.cr)^3)
    num2<-((muhat-1)^2)*exp(alpha)*log(1-tau)*critical.ly
    den2<-(muhat^2)*(e.am-exp(alpha))*((den.cr)^2)
    num3<-((muhat-1)^2)*exp(2*alpha)*log(1-tau)*critical.ly
    den3<-(muhat^2)*((exp(alpha)-e.am)^2)*((den.cr)^2)
    num4<-((muhat-1)^2)*exp(2*alpha)
    den4<-den3
    num5<-((muhat-1)^2)*exp(alpha)
    den5<-(muhat^2)*(e.am-exp(alpha))*den.cr
    num6<-num4
    den6<-(muhat^2)*((exp(alpha)-e.am)^2)*den.cr
    num7<-((y-1)^2)*critical.y*(log(1-tau)/den.cr-1)
    den7<-(y^2)*(1-critical.y)
    num8<-((y-1)^2)*exp(2*alpha-2*alpha/y)*(log(1-tau)/den.cr-1)
    den8<-(y^2)*((critical.y-1)^2)
    num9<-2*(muhat-1)*exp(2*alpha)*(y-1)*log(1-tau)
    den9<-muhat*y*(exp(alpha)-e.am)*(e.ay-exp(alpha))*((den.cr)^2)
    Uaa<-ifelse(y==0 | y==1,0,num1/den1+num2/den2+num3/den3+num4/den4+num5/den5+num6/den6-num7/den7-num8/den8+num9/den9-1/(alpha^2))
    
    # print(Uaa)
    Ualphaalpha<-sum(Uaa)
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO ALPHA
    
    ####START DERIVATIVE FROM [DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU] IN RESPECT TO ALPHA
    ###########################################################################################################
    e.a=exp(alpha)
    numerador<-(e.a*(den.cr*(log(1-tau)*(muhat*(-e.a)*alpha*(y-1)*(e.a-e.am)-y*(muhat*e.a-e.am*(muhat*alpha+muhat-alpha))*(e.a-e.ay)*critical.ly)+(muhat-1)*e.a*alpha*y*(e.a-e.ay))+2*(muhat-1)*e.a*alpha*y*(e.a-e.ay)*log(1-tau)*critical.ly+y*(muhat*e.a-e.am*(muhat*alpha+muhat-alpha))*(-(e.a-e.ay))*((den.cr)^2)))
    denominador<-(muhat^3)*y*((e.a-e.am)^2)*(e.a-e.ay)*((den.cr)^3)
    Umualpha<-ifelse(y==0| y==1,0,numerador/denominador)
    #print(Umualpha)
    ########################################################################################################### 
    ####END DERIVATIVE FROM [DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU] IN RESPECT TO ALPHA
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU
    ###########################################################################################################
    num1<-exp(alpha)*alpha*(den.cr+log(1-tau)*critical.ly)
    den1<-(muhat^2)*(exp(alpha/muhat)-exp(alpha))*((den.cr)^2)
    mustar<-ifelse(y==0| y==1,0,num1/den1)
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU 
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA
    ###########################################################################################################
    lambdastar<-ifelse(y==0,1/lambdahat,NA)
    lambdastar<-ifelse(y==1,1/lambdahat,lambdastar)
    lambdastar<-ifelse(y!=0 & y!=1,-1/(1-lambdahat),lambdastar)
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA 
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO p
    ###########################################################################################################
    pstar<-ifelse(y==1,1/phat,NA)
    pstar<-ifelse(y==0,-1/(1-phat),pstar)
    pstar<-ifelse(y!=0 & y!=1,0,pstar)
    # print("pstar");print(pstar)
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO p    
    mT <- diag(as.vector(mu.eta(X%*%as.matrix(beta)))) 
    num.mT2<-dif2link(muhat)
    mT2<-diag(as.vector(-num.mT2*(mu.eta(X%*%as.matrix(beta))^3)))
    mV <- diag(as.vector(muhatstar.sec))
    mV0 <- diag(as.vector(mustar))
    mualpha<-diag(as.vector(Umualpha))
    vI <- matrix(rep(1,n),ncol=1)
    lT <- diag(as.vector(mu.eta(Z%*%as.matrix(nu))))
    num.lT2<-dif2link(lambdahat)
    lT2<-diag(as.vector(-num.lT2*(mu.eta(Z%*%as.matrix(nu))^3)))
    lV <- diag(as.vector(lambdahatstar.sec))
    lV0 <- diag(as.vector(lambdastar))
    pT <- diag(as.vector(mu.eta(A%*%as.matrix(rho))))
    num.pT2<-dif2link(phat)
    pT2<-diag(as.vector(-num.pT2*(mu.eta(A%*%as.matrix(rho))^3)))
    pV <- diag(as.vector(phatstar.sec))
    pV0 <- diag(as.vector(pstar))
    
    # print("pT");print(summary(as.vector(pT)))
    # print("num.pT2");print(summary(as.vector(num.pT2)))
    # print("pT2");print(summary(as.vector(pT2)))
    # print("pV");print(summary(as.vector(pV)))
    # print("pV0");print(summary(as.vector(pV0)))
    KBB=matrix(rep(NA,(k+1)*(k+1)),ncol=(k+1))
    if(length(KBB)==1){
      KBB <- -(t(X)%*%mV%*%(mT^2)%*%X + t(X)%*%mV0%*%mT2%*%X + t(vI)%*%mV0%*%mT%*%deta.dbetabeta)
    }else{
      for(j in 1:(k+1)){
        for(i in 1:(k+1)){
          KBB[i,j] <- -(t(as.vector(X[,i]))%*%mV%*%(mT^2)%*%as.vector(X[,j]) + t(as.vector(X[,i]))%*%mV0%*%mT2%*%as.vector(X[,j]) + t(vI)%*%mV0%*%mT%*%as.matrix(deta.dbetabeta[i,j,]))
        }
      }
    }
    KBalpha <- -t(X)%*% mualpha %*% mT %*% vI
    KBnu=matrix(rep(0,(k+1)*(c+1)),ncol=(c+1))
    KBrho=matrix(rep(0,(k+1)*(m+1)),ncol=(m+1))
    
    KalphaB <- t(KBalpha)
    Kalphaalpha <- -Ualphaalpha
    Kalphanu<-t(as.matrix(rep(0,c+1)))
    Kalpharho<-t(as.matrix(rep(0,m+1)))
    KnuB<-t(KBnu)
    Knualpha<-t(Kalphanu)
    Knunu<-matrix(rep(NA,(c+1)*(c+1)),ncol=(c+1))
    if(length(Knunu)==1){
      Knunu <- -(t(Z)%*%lV%*%(lT^2)%*%Z + t(Z)%*%lV0%*%lT2%*%Z + t(vI)%*%lV0%*%lT%*%deta.dnunu)
    }else{
      for(j in 1:(c+1)){
        for(i in 1:(c+1)){
          Knunu[i,j] <- -(t(as.vector(Z[,i]))%*%lV%*%(lT^2)%*%as.vector(Z[,j]) + t(as.vector(Z[,i]))%*%lV0%*%lT2%*%as.vector(Z[,j]) + t(vI)%*%lV0%*%lT%*%as.matrix(deta.dnunu[i,j,]))
        }
      }
    }
    Knurho<-matrix(rep(0,(c+1)*(m+1)),ncol=(m+1))
    
    KrhoB<-t(KBrho)
    Krhoalpha<-t(Kalpharho)
    Krhonu<-t(Knurho)
    Krhorho<-matrix(rep(NA,(m+1)*(m+1)),ncol=(m+1))
    if(length(Krhorho)==1){
      Krhorho <- -(t(A)%*%pV%*%(pT^2)%*%A + t(A)%*%pV0%*%pT2%*%A + t(vI)%*%pV0%*%pT%*%deta.drhorho)
    }else{
      for(j in 1:(m+1)){
        for(i in 1:(m+1)){
          Krhorho[i,j] <- -(t(as.vector(A[,i]))%*%pV%*%(pT^2)%*%as.vector(A[,j]) + t(as.vector(A[,i]))%*%pV0%*%pT2%*%as.vector(A[,j]) + t(vI)%*%pV0%*%pT%*%as.matrix(deta.drhorho[i,j,]))
        }
      }
    }
    # print(KBB);print(KBalpha);print(KBnu)
    # print(KalphaB);print(Kalphaalpha);print(Kalphanu)
    # print(KnuB);print(Knualpha);print(Knunu)
    K <- rbind(
      cbind(KBB,KBalpha,KBnu,KBrho),
      cbind(KalphaB,Kalphaalpha,Kalphanu,Kalpharho),
      cbind(KnuB,Knualpha,Knunu,Knurho),
      cbind(KrhoB,Krhoalpha,Krhonu,Krhorho)
    )
    return(K)
  }
  K<-obs.inf(y,muhat)
  #print(K)
  Ksolve<- tryCatch(solve(K), error = function(e) return("error"))
  if(Ksolve[1] == "error")
  {z$RMC=1#used at Monte-Carlo simulation for discard from the sample
  warning("Analytic Observed Information Matrix is not positive semi-definite")
  return(z)#if Analytic Observed Information Matrix is not positive semi-definite, do not calculate
  }else{sol=try(solve(K))}
  
  v<-diag(sol)#Variância assintótica dos esimadores
  #print(v)
  for (i in 1:length(v))
  {
    if(is.na(v[i]) | is.nan(v[i]) | v[i]<0 )  {
      z$RMC=1
      warning("Analytic Observed Information Matrix is not positive semi-definite")
      return(z)#if Analytic Observed Information Matrix is not positive semi-definite, do not calculate
    }
  }
  
  z$zstat<-z$coeff/sqrt(v)
  #print(z$zstat)
  resp<-rep(0,length(z$zstat))
  for (i in 1:length(resp)){
    if(abs(z$zstat[i])>qnorm(0.975))
    {
      resp[i] <- "H0 rejected"
    } else {resp[i] <- "H0 not rejected"}
  }
  LI<-z$coeff-qnorm(0.975)*sqrt(v)
  LS<-z$coeff+qnorm(0.975)*sqrt(v)
  z$LI<-LI
  z$LS<-LS
  z$pvalues<-(1-pnorm(abs(z$zstat)))*2
  first_col<-c(0:k,"alpha",0:c,0:m)
  result <- matrix(c(first_col,round(c(z$coeff,z$zstat,LI,LS,z$pvalues),4),resp), nrow=length(z$coeff), ncol=7, byrow=F)
  colnames(result) <- c("Parameter","MLE","Wald's Statistic","Lower bound","Upper bound","p-value","Wald'S Test result")
  rownames(result)<-c(rep("beta",(k+1)),"",rep("nu",(c+1)),rep("rho",(m+1)))
  z$coef.result<-result
  z$aic <- -2*(z$loglik)+2*(length(opt$par)) 
  z$bic <- -2*(z$loglik)+(length(opt$par))*log(n)
  result2<-matrix(round(c(z$loglik,z$aic,z$bic),4),nrow=3,ncol=1)
  rownames(result2)<-c("Log-likelihood","AIC","BIC")
  #print(z$coef.result)
  ###START error metrics
  mae<-sum(abs(y[1:n]-z$fitted[1:n]))/n
  
  sq<-rep(NA,n)
  for(i in 1:n)
  {
    sq[i]<-(y[i]-z$fitted[i])^2
  }  
  
  mse<-sum(sq[1:n])/n
  
  rmse<-sqrt(mse)
  
  mape<-sum(abs((y[1:n]-z$fitted[1:n])/y[1:n])*100)/n
  
  MdRAE<-median(abs(y[2:n]-z$fitted[2:n])/abs(y[2:n]-y[1:(n-1)]))
  
  MAEnaive<-sum(abs(y[2:n]-y[1:(n-1)]))/(n-1)
  MASE<-mae/MAEnaive #Its value greater than one (1) indicates the algorithm is performing poorly compared to the naïve forecast.
  
  MAEnaive.star<-sum(abs(y[2:n]-y[1:(n-1)]))/n
  MASE.star<-mae/MAEnaive.star
  
  #Mean directional accuracy
  sign.y<-sign(y[2:n]-y[1:(n-1)])
  sign.f<-sign(z$fitted[2:n]-y[1:(n-1)])
  MDA.cont<-0
  for (i in 1:(n-1)){  
    if(sign.y[i]==sign.f[i]){MDA.cont<-MDA.cont+1}  
  }
  MDA<-MDA.cont/n
  
  MASE=MASE.star
  
  z$accuracyfitted<-accuracy<-matrix(round(c(mae,mse,rmse,mape,MdRAE,MASE,MDA),4), nrow=1, ncol=7, byrow=T)
  colnames(z$accuracyfitted) <-colnames(accuracy) <- c("MAE","MSE","RMSE","MAPE","MdRAE", "MASE","MDA")
  rownames(accuracy) <- c("")
  # print(z$accuracyfitted)
  ytofit<-y
  
  
  ###########################
  
  z$serie <- y

  muhat[is.na(muhat)]<-.Machine$double.eps
  muhat[muhat<.Machine$double.eps]<-.Machine$double.eps
  muhat[muhat>0.9999999]<-0.9999999
  critical<-(alpha-alpha/muhat)
  critical[is.na(critical)]<--.Machine$double.eps
  critical[is.nan(critical)]<--36.04365
  critical[critical< (-36.04365)]<--36.04365
  den.cr=log(1-exp(critical))
  den.cr[is.nan(den.cr)]<--36.04365
  
  ########################################################################
  ######## randomized quantile residuals with uniform distribution  ######
  ########################################################################
  
  ui<-rep(NA,n)
  for(i in 1:n)
  {
    if(y[i]==0) ui[i] <- runif(1,0,lambdahat[i]*(1-phat[i]))
    if(y[i]==1) ui[i] <- runif(1,(1-lambdahat[i]*phat[i]),1)
    if(y[i]!=0 & y[i]!=1) ui[i] <- pmk(y[i],alpha=z$alpha,beta=log(1-tau)/den.cr[i], lambda=lambdahat[i],p=phat[i],log.p = FALSE)
  }
  z$residual <- residual <- qnorm(ui)

  mresult<-matrix(round(c(z$loglik,z$aic,z$bic),4),nrow=3,ncol=1)
  rownames(mresult)<-c("Log-likelihood","AIC","BIC")
  colnames(mresult)<-c("")
  z$mresult<-mresult
  
 
  
  ###################################################
  ######### GRAPHICS ################################
  
  if(graph==T)
  {
    t<-seq(-5,n+6,by=1)
    w1<-5
    h1<-4
    pdf("cor_y_covX.pdf",width=5, height=4)
    {
      if (k<=3){
        par(mfrow=c(1,k))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))}
      # for (i in k){
      # plot(as.vector(X[,1+i]),as.vector(y),main=" ",xlab=paste("Cov ",i),ylab="y")
      # }
      if(k==1){
        plot(as.vector(X[,2]),as.vector(y),main=" ",xlab=paste("Cov X1"),ylab="y")#,legend("topleft",paste("cor=",cor(X[,2],y))))#,#pch=vpch,
        text(x = median(as.vector(X[,2])), y = median(as.vector(y)), label = paste("cor=",round(cor(X[,2],y),4)), cex = 1)
      }
      if (k==2){
        plot(as.vector(X[,2]),as.vector(y),main=" ",xlab=paste("Cov X1"),ylab="y")#,legend("topleft",paste("cor=",cor(X[,2],y))))
        text(x = median(as.vector(X[,2])), y = median(as.vector(y)), label = paste("cor=",round(cor(X[,2],y),4)), cex = 1)
        
        plot(as.vector(X[,3]),as.vector(y),main=" ",xlab=paste("Cov X2"),ylab="y")
        text(x = median(as.vector(X[,3])), y = median(as.vector(y)), label = paste("cor=",round(cor(X[,3],y),4)), cex = 1)
      }
      if (k==3){
        plot(as.vector(X[,2]),as.vector(y),main=" ",xlab=paste("Cov X1"),ylab="y")#,legend("topleft",paste("cor=",cor(X[,2],y)),pt.bg="white", lty=c(1,2), bty="n"))
        text(x = median(as.vector(X[,2])), y = median(as.vector(y)), label = paste("cor=",round(cor(X[,2],y),4)), cex = 1)
        
        plot(as.vector(X[,3]),as.vector(y),main=" ",xlab=paste("Cov X2"),ylab="y")
        text(x = median(as.vector(X[,3])), y = median(as.vector(y)), label = paste("cor=",round(cor(X[,3],y),4)), cex = 1)
        
        plot(as.vector(X[,4]),as.vector(y),main=" ",xlab=paste("Cov X3"),ylab="y")
        text(x = median(as.vector(X[,4])), y = median(as.vector(y)), label = paste("cor=",round(cor(X[,4],y),4)), cex = 1)
      }
    }
    dev.off()
    if(k>1){
      pdf("cor_covX_covX.pdf",width=5, height=4)
      {
        if (k==2){
          par(mfrow=c(1,1))
          par(mar=c(2.8, 2.7, 1, 1))
          par(mgp=c(1.7, 0.45, 0))
          
          plot(as.vector(X[,3]),as.vector(X[,2]),main=" ",xlab=paste("Cov X2"),ylab="Cov X1")
          text(x = median(as.vector(X[,3])), y = median(as.vector(X[,2])), label = paste("cor=",round(cor(X[,3],X[,2]),4)), cex = 1)
        }
        if (k==3){
          par(mfrow=c(1,3))
          par(mar=c(2.8, 2.7, 1, 1))
          par(mgp=c(1.7, 0.45, 0))
          
          plot(as.vector(X[,3]),as.vector(X[,2]),main=" ",xlab=paste("Cov X2"),ylab="Cov X1")
          text(x = median(as.vector(X[,3])), y = median(as.vector(X[,2])), label = paste("cor=",round(cor(X[,3],X[,2]),4)), cex = 1)
          
          plot(as.vector(X[,4]),as.vector(X[,2]),main=" ",xlab=paste("Cov X3"),ylab="Cov X1")
          text(x = median(as.vector(X[,4])), y = median(as.vector(X[,2])), label = paste("cor=",round(cor(X[,4],X[,2]),4)), cex = 1)
          
          plot(as.vector(X[,4]),as.vector(X[,3]),main=" ",xlab=paste("Cov X3"),ylab="Cov X2")
          text(x = median(as.vector(X[,4])), y = median(as.vector(X[,3])), label = paste("cor=",round(cor(X[,4],X[,3]),4)), cex = 1)
          
        }
      }
      dev.off()
      
      pdf("cor_matrix_covX.pdf",width=5, height=4)
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))
        library(corrplot)
        corrplot(cor(exvar.beta),method='number')              # visualize the multicollinearity
      }
      dev.off()
    }
    pdf("cor_y_covZ.pdf",width=5, height=4)
    {
      if (c<=3){
        par(mfrow=c(1,c))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))}
      # for (i in k){
      # plot(as.vector(X[,1+i]),as.vector(y),main=" ",xlab=paste("Cov ",i),ylab="y")
      # }
      if(c==1){
        plot(as.vector(Z[,2]),as.vector(y),main=" ",xlab=paste("Cov Z1"),ylab="y")#,legend("topleft",paste("cor=",cor(X[,2],y))))#,#pch=vpch,
        text(x = median(as.vector(Z[,2])), y = median(as.vector(y)), label = paste("cor=",round(cor(Z[,2],y),4)), cex = 1)
      }
      if (c==2){
        plot(as.vector(Z[,2]),as.vector(y),main=" ",xlab=paste("Cov Z1"),ylab="y")#,legend("topleft",paste("cor=",cor(X[,2],y))))
        text(x = median(as.vector(Z[,2])), y = median(as.vector(y)), label = paste("cor=",round(cor(Z[,2],y),4)), cex = 1)
        
        plot(as.vector(Z[,3]),as.vector(y),main=" ",xlab=paste("Cov Z2"),ylab="y")
        text(x = median(as.vector(Z[,3])), y = median(as.vector(y)), label = paste("cor=",round(cor(Z[,3],y),4)), cex = 1)
      }
      if (c==3){
        plot(as.vector(Z[,2]),as.vector(y),main=" ",xlab=paste("Cov Z1"),ylab="y")#,legend("topleft",paste("cor=",cor(X[,2],y)),pt.bg="white", lty=c(1,2), bty="n"))
        text(x = median(as.vector(Z[,2])), y = median(as.vector(y)), label = paste("cor=",round(cor(Z[,2],y),4)), cex = 1)
        
        plot(as.vector(Z[,3]),as.vector(y),main=" ",xlab=paste("Cov Z2"),ylab="y")
        text(x = median(as.vector(Z[,3])), y = median(as.vector(y)), label = paste("cor=",round(cor(Z[,3],y),4)), cex = 1)
        
        plot(as.vector(Z[,4]),as.vector(y),main=" ",xlab=paste("Cov Z3"),ylab="y")
        text(x = median(as.vector(Z[,4])), y = median(as.vector(y)), label = paste("cor=",round(cor(Z[,4],y),4)), cex = 1)
      }
    }
    dev.off()
    if(c>1){
      pdf("cor_covZ_covZ.pdf",width=5, height=4)
      {
        if (c==2){
          par(mfrow=c(1,1))
          par(mar=c(2.8, 2.7, 1, 1))
          par(mgp=c(1.7, 0.45, 0))
          
          plot(as.vector(Z[,3]),as.vector(Z[,2]),main=" ",xlab=paste("Cov Z2"),ylab="Cov Z1")
          text(x = median(as.vector(Z[,3])), y = median(as.vector(Z[,2])), label = paste("cor=",round(cor(Z[,3],Z[,2]),4)), cex = 1)
        }
        if (c==3){
          par(mfrow=c(1,3))
          par(mar=c(2.8, 2.7, 1, 1))
          par(mgp=c(1.7, 0.45, 0))
          
          plot(as.vector(Z[,3]),as.vector(Z[,2]),main=" ",xlab=paste("Cov Z2"),ylab="Cov Z1")
          text(x = median(as.vector(Z[,3])), y = median(as.vector(Z[,2])), label = paste("cor=",round(cor(Z[,3],Z[,2]),4)), cex = 1)
          
          plot(as.vector(Z[,4]),as.vector(Z[,2]),main=" ",xlab=paste("Cov Z3"),ylab="Cov Z1")
          text(x = median(as.vector(Z[,4])), y = median(as.vector(Z[,2])), label = paste("cor=",round(cor(Z[,4],Z[,2]),4)), cex = 1)
          
          plot(as.vector(Z[,4]),as.vector(Z[,3]),main=" ",xlab=paste("Cov Z3"),ylab="Cov Z2")
          text(x = median(as.vector(Z[,4])), y = median(as.vector(Z[,3])), label = paste("cor=",round(cor(Z[,4],Z[,3]),4)), cex = 1)
          
        }
      }
      dev.off()
      
      pdf("cor_matrix_covZ.pdf",width=5, height=4)
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))
        library(corrplot)
        corrplot(cor(exvar.nu),method='number')              # visualize the multicollinearity
      }
      dev.off()
    }

    pdf("resid_v_fitted.pdf",width=5, height=4)
    {
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
      plot(as.vector(z$fitted[1:n]),as.vector(residual), main=" ", pch = "+",
           xlab="Fitted values",ylab="Residuals",ylim=c(-4,4))
      lines(t,rep(-3,n+12),lty=2,col=1)
      lines(t,rep(3,n+12),lty=2,col=1)
      lines(t,rep(-2,n+12),lty=3,col=1)
      lines(t,rep(2,n+12),lty=3,col=1)
    }
    dev.off()
    pdf("obs_v_fit.pdf",width=5, height=4)### abre no navegador google chrome só
    {
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
      plot(as.vector(z$fitted), as.vector(ytofit), main=" ", pch = "+",
           xlab="Fitted values",ylab="Observed data",
           xlim=c(0.95*min(y),max(y)*1.05),
           ylim=c(0.95*min(y),max(y)*1.05))
      lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    }
    dev.off()
    pdf("resid_density.pdf",width=5, height=4)
    {
      par(mfrow=c(1,1))
      par(mar=c(1.5, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
      densidade<-density(residual)
      plot(densidade,ylab="Density",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
      lines(densidade$x,dnorm(densidade$x),lty=2)
      legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
             pt.bg="white", lty=c(1,2), bty="n")
    }
    dev.off()
    
    pdf("qq_plot.pdf",width=5, height=4)
    {  
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
      qqnorm(residual, pch = "+",
             xlim=c(0.95*min(residual),max(residual)*1.05),
             ylim=c(0.95*min(residual),max(residual)*1.05),
             main="",xlab="Normal quantile",ylab="Empirical quantile")
      lines(c(-10,10),c(-10,10),lty=2)
    }
    dev.off()
    # pdf("envelope_plot.pdf",width=5, height=4)
    # { 
    #   par(mfrow=c(1,1))
    #   par(mar=c(2.8, 2.7, 1, 1)) 
    #   par(mgp=c(1.7, 0.45, 0))
    #   library(hnp)
    #   hnp(residual,half=F, sim = 1000, conf = 0.95, ylab="Empirical quantile", xlab="Normal quantile")
    # }
    # dev.off()
    pdf("adjusted.pdf",width=5, height=4)
    {
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
      par(mgp=c(1.7, 0.45, 0))
      plot(ytofit,type="l",ylab="Serie",xlab="Time")
      lines(z$fitted,col="blue",lty=2)
      legend("bottomleft",c("Observed data","Fitted values"),#pch=vpch,
             pt.bg="white", lty=c(1,2), bty="n",col=c(1,"blue"))
    }
    dev.off()
    
  }#END GRAPHICS
  
  
  ########################################################################
  ########################   residual analysis   ########################
  ########################################################################
  loglik_null <- function(z)
  {
    beta <- z[1]
    alpha <- z[2]
    nu<-z[3]
    rho<-z[4]
    lambda<-linkinv(Z[,1]*nu)
    p<-linkinv(A[,1]*rho)
    mu <- linkinv(X[,1]*beta)
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
    l=ifelse(y==0,log(lambda)+log(1-p),NA)
    l=ifelse(y==1,log(lambda)+log(p),l)
    l=ifelse(y!=0 & y!=1,log(1-lambda)+log(alpha)+alpha-alpha/y+(log(1-tau)/den.cr -1)*critical.ly-2*log(y)+log(log(1-tau)/den.cr),l)
    return(sum(l))
  }
  
  # loglik_null <- function(z)
  # {
  #   beta <- z[1]
  #   alpha <- z[2]
  #   lambda_null<-length(which(y==0 | y==1))/n
  #  # p_null<-length(which(y==1))/length(which(y==0 | y==1))
  #   mu <- linkinv(X[,1]*beta)
  #   mu[is.na(mu)]<-.Machine$double.eps
  #   mu[mu<.Machine$double.eps]<-.Machine$double.eps
  #   mu[mu>0.9999999]<-0.9999999
  #   critical<-alpha-alpha/mu
  #   critical[is.na(critical)]<--.Machine$double.eps
  #   critical[is.nan(critical)]<--36.04365
  #   critical[critical< (-36.04365)]<--36.04365
  #   den.cr=log(1-exp(critical))
  #   den.cr[is.nan(den.cr)]<--36.04365
  #   critical.y<-exp(alpha-alpha/y)
  #   critical.y[is.infinite(critical.y)]<-.Machine$double.xmax
  #   critical.ly<-log(1-critical.y)
  #   critical.ly[is.nan(critical.ly)]<--36.04365
  #   critical.ly[critical.ly< (-36.04365)]<--36.04365#para exp dar .Machine$double.eps
  #   l=ifelse(y!=0 & y!=1,log(1-lambda_null)+log(alpha)+alpha-alpha/y+(log(1-tau)/den.cr -1)*critical.ly-2*log(y)+log(log(1-tau)/den.cr),0)
  #   return(sum(l))
  # }
  # 
  # ini_null<- c(linkfun(mean(y[y!=0&y!=1])),reg[k+2])
  
  
  ini_null<- c(linkfun(mean(y[y!=0&y!=1])),reg[k+2],length(y[y==0|y==1])/n,length(y[y==1])/length(y[y==0|y==1]))
  # print(ini_null)
  #opti.error<- tryCatch(optim(ini_null, loglik_null,method = "BFGS", control = list(fnscale = -1)), error = function(e) return("error"))
  # if(opti.error[1] == "error")
  # {z$r2 <-NA
  # }else{
  opt_null <- optim(ini_null, loglik_null,method = "BFGS", control = list(fnscale = -1)) # , maxit = 500, reltol = 1e-9))
  #opt_null<-optim(alpha,on.dmk.alpha,method = "BFGS")
  r2 <- 1-exp((-2/n)*(opt$value-opt_null$value))
  z$r2 <- r2
  #}
  # print(z$r2)
  #null hypothesis: normality
  library(nortest)
  andersondarling=ad.test(residual)
  z$andersondarling<-andersondarling$statistic
  z$p_andersondarling<-andersondarling$p.value
  
  if(n<=5000){
    shapiro=shapiro.test(residual)
    z$shapiro=shapiro$statistic
    z$p_shapiro=shapiro$p.value
  }
  # print(z$p_andersondarling)
  
  #non autocorrelated
  ljungbox<- Box.test(residual, lag = 10, type = "Ljung-Box", fitdf = (k+c))
  z$ljungbox<-ljungbox$statistic
  z$p_ljungbox<-ljungbox$p.value
  
  dw<-function(res){
    alternative = "two.sided"
    dw <- sum(diff(res)^2)/sum(res^2)
    Q1 <- chol2inv(qr.R(qr(cbind(exvar.beta,exvar.nu))))
    if (n<100) {
      A <- diag(c(1, rep(2, n - 2), 1))
      A[abs(row(A) - col(A)) == 1] <- -1
      MA <- diag(rep(1, n)) - cbind(exvar.beta,exvar.nu) %*% Q1 %*% t(cbind(exvar.beta,exvar.nu))
      MA <- MA %*% A
      ev <- eigen(MA)$values[1:(n - (c+k))]
      if (any(Im(ev) > 1e-10)) 
        warning("imaginary parts of eigenvalues discarded")
      ev <- Re(ev)
      ev <- ev[ev > 1e-10]
      pdw <- function(dw) .Fortran("pan", as.double(c(dw, ev)), as.integer(length(ev)), as.double(0), 
                                   as.integer(15), x = double(1), PACKAGE = "lmtest")$x
      pval <- switch(alternative, two.sided = (2 * min(pdw(dw), 1 - pdw(dw))), less = (1 - pdw(dw)), greater = pdw(dw))
      if (is.na(pval) || ((pval > 1) | (pval < 0))) {
        warning("exact p value cannot be computed (not in [0,1]), approximate p value will be used")
      }
    }else{
      if (n < max(5, (k+c))) {
        warning("not enough observations for computing an approximate p value, set to 1")
        pval <- 1
      }
      else {
        AX <- matrix(as.vector(filter(cbind(exvar.beta,exvar.nu), c(-1, 2, -1))), ncol = (c+k))
        AX[1, ] <- cbind(exvar.beta,exvar.nu)[1, ] - cbind(exvar.beta,exvar.nu)[2, ]
        AX[n, ] <- cbind(exvar.beta,exvar.nu)[n, ] - cbind(exvar.beta,exvar.nu)[(n - 1), ]
        XAXQ <- t(cbind(exvar.beta,exvar.nu)) %*% AX %*% Q1
        P <- 2 * (n - 1) - sum(diag(XAXQ))
        Q <- 2 * (3 * n - 4) - 2 * sum(diag(crossprod(AX) %*% Q1)) + sum(diag(XAXQ %*% XAXQ))
        dmean <- P/(n - (c+k))
        dvar <- 2/((n - (c+k)) * (n - (c+k) + 2)) * (Q - P * dmean)
        pval <- switch(alternative, two.sided = (2 * pnorm(abs(dw - dmean), sd = sqrt(dvar), lower.tail = FALSE)), less = pnorm(dw, mean = dmean, sd = sqrt(dvar), lower.tail = FALSE), greater = pnorm(dw, mean = dmean, sd = sqrt(dvar)))
      }
    }
    alternative <- switch(alternative, two.sided = "true autocorrelation is not 0", less = "true autocorrelation is less than 0", greater = "true autocorrelation is greater than 0")
    names(dw) <- "DW"
    RVAL <- list(statistic = dw, method = "Durbin-Watson test", 
                 alternative = alternative, p.value = pval)
    class(RVAL) <- "htest"
    return(RVAL)
  }
  durbin<-dw(residual)
  z$durbin<-durbin$statistic
  z$p_durbin<-durbin$p.value
  #null hypothesis: non-multicollinearity
  z$cor.beta <- cor(exvar.beta)                                         # independent variables correlation matrix 
  z$cor.nu <- cor(exvar.nu)
  #null hypothesis: non-heteroscedasticity (constant variance)
   SBP<-function(resi){#adapted to fits a linear regression model to the residuals of the mkreg model 
    sigma2 <- sum(resi^2)/length(resi)
    w <- resi^2 - sigma2
    aux <- lm.fit(exvar.beta, w)
    bp <- n * sum(aux$fitted.values^2)/sum(w^2)
    method <- "studentized Breusch-Pagan test"
    names(bp) <- "BP"
    df <- c(df = aux$rank - 1)
    RVAL <- list(statistic = bp, parameter = df, method = method, 
                 p.value = pchisq(bp, df, lower.tail = FALSE))
    class(RVAL) <- "htest"
    return(RVAL)
  }
  breusch=SBP(residual)
  z$breusch=breusch$statistic
  z$p_breusch=breusch$p.value
  
  rownames(z$accuracyfitted) <-rownames(accuracy) <- c("Accuracy fitted")
  if(n<=5000){
    diagnostic<-matrix(round(c(z$andersondarling,z$shapiro,z$durbin,z$ljungbox,z$breusch,
                               z$p_andersondarling,z$p_shapiro,z$p_durbin,z$p_ljungbox,z$p_breusch
    ),4), nrow=2, ncol=5, byrow=T)
    colnames(diagnostic) <- c("Anderson-Darling test","Shapiro-Wilk test","Durbin-Watson test","Ljung-Box test","studentized Breusch-Pagan test")
  }else{
    diagnostic<-matrix(round(c(z$andersondarling,z$durbin,z$ljungbox,z$breusch,
                               z$p_andersondarling,z$p_durbin,z$p_ljungbox,z$p_breusch
    ),4), nrow=2, ncol=4, byrow=T)
    colnames(diagnostic) <- c("Anderson-Darling test","Durbin-Watson test","Ljung-Box test","studentized Breusch-Pagan test")
    
  }
  rownames(diagnostic) <- c("Statistic","P-value")
  z$diagnostic=diagnostic
  
  if(print==T){
    print("MKreg",quote=F)
    print(z$coef.result,quote=F)
    message("")
    print(c("Log-likelihood =",round(z$loglik,4)),quote=F)
    print(c("aic =",round(z$aic,4),"bic =",round(z$bic,4)),quote=F)
    print(c("R-squared=",round(z$r2,4)),quote=F)
    message("")  
    print("Randomized quantile residuals:",quote=F)
    print(summary(z$residual))
    message("")
    print(z$diagnostic)
    message("")
    print(z$accuracyfitted)
    message("")
  }
  
  if(check==TRUE){
    opt2 <- optim(reg, loglik, method = "BFGS", hessian = T, control = list(fnscale = -1))    
    library(rootSolve)
    print("mle estimação com gradiente analítico")
    print(opt$par)
    print("mle estimação com gradiente numérico")
    print(opt2$par)
    print("verificando derivadas")
    print(gradient(score,opt2$par))
    print(hessian(loglik,opt2$par))
    print("rbind(score(opt$par),gradient(loglik,opt$par))")
    print(rbind(score(opt$par),gradient(loglik,opt$par)))
    print("rbind(score(opt2$par),gradient(loglik,opt2$par))")
    print(rbind(score(opt2$par),gradient(loglik,opt2$par)))
    print("-hessiana = Matriz de informação observada condicional")
    print(round(K,4))
    print("hessiana numerica")
    print(round(-opt$hessian,4))
    print("comparando meu cálculo com hessiana da estimação numérica")
    print(round((K+opt2$hessian),2))
    print("comparando meu cálculo com hessiana numérica da estimação analítica")
    print(round((K+opt$hessian),2))
    print("soma diferença hessiana otimização numérica")
    print(round(sum(abs(K+opt2$hessian)),2))
    print("soma diferença hessiana numérica otimização analítica")
    print(round(sum(abs(K+opt$hessian)),2))
  }
  
  return(z)
}#fim estimação
