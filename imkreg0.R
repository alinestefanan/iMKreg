# y: time series
# S: seasonal periodicity
# exvar.beta,exvar.nu: covariate column matrix
# tau: quantil, when set 0.5 is the median
# resid: 1 = quantile residual, 2 = deviance residual
# link: "logit", "probit" or "cloglog"

mkreg <- function(y,exvar.beta=NA,exvar.nu=NA,tau=0.5,resid=1,graph=T,print=T,check=F,link="logit")
{
  n <- length(y) 

  exvar.nu=as.matrix(exvar.nu)
  c=if(is.matrix(exvar.nu)){ncol(exvar.nu)}else{1}
  Z<-matrix(c(rep(1,n),exvar.nu), nrow=n, ncol=(c+1), byrow=F)
  
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
  dmk <- Vectorize(function(y,alpha,beta,lambda0,log = FALSE){
    critical.y=exp(alpha-alpha/y)
    density=c()
    for (i in 1:length(y)){
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
      if (y[i]==0){
        density[i]=lambda0[i]
      }else{
        density[i]<-(1-lambda0[i])*(alpha*beta[i]*critical.y[i]*(1-critical.y[i])^(beta[i]-1))/(y[i]^2)
        
      }
    }
    density[is.na(density)] <- .Machine$double.eps
    density[is.nan(density)] <- .Machine$double.eps
    
    for (i in 1:length(density)){
      if(density[i]<.Machine$double.eps) density[i]<-.Machine$double.eps
      if(density[i]>.Machine$double.xmax) density[i]<-.Machine$double.xmax
    }
   logden <- log(density)
    val <- ifelse(log, logden, exp(logden)) 
    return(val)
  }) 
  
  #modified kumaraswamy cumulative inflated distribution function 
  pmk <- Vectorize(function(q,alpha,beta,lambda0,log.p = FALSE){
    cdf=rep(0,length(q))
    for(i in length(q)){
      if (y[i]==0){
        cdf[i] <- lambda0[i]
      }else{
      cdf[i] <-  lambda0[i]+(1-lambda0[i])*(1-(1-exp(alpha-alpha/q[i]))^beta[i])
      if (is.na(cdf[i])) cdf[i]<-.Machine$double.eps
      if (cdf[i]<.Machine$double.eps) cdf[i]<-.Machine$double.eps
      if (cdf[i]>0.9999999) cdf[i]<-0.9999999#1-.Machine$double.eps
    }}
    val <- ifelse(log.p, log(cdf), cdf)
    return(val)
  })
  
  dmk_alpha<-function(alpha){
    beta=log(1-tau)/log(1-exp(alpha-alpha/median(y)))
    lambda0=length(y[which(y==0)])/length(y)
    if (is.na(beta)){beta=.Machine$double.eps}
    if(beta<.Machine$double.eps) beta<-.Machine$double.eps
    if (is.infinite(beta)){beta=.Machine$double.xmax}
    critical.y=exp(alpha-alpha/y)
    density=c()
    for (i in 1:length(y)){
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
      if (y[i]==0){
        density[i]=lambda0
      }else{
        density[i]<-(1-lambda0)*(alpha*beta*critical.y[i]*(1-critical.y[i])^(beta-1))   /(y[i]^2)
    }}
    
    density[is.na(density)] <- .Machine$double.eps
    density[is.nan(density)] <- .Machine$double.eps
    for (i in 1:length(density)){
      if(density[i]<.Machine$double.eps) density[i]<-.Machine$double.eps
    }
    return(density)
  }
 
  # função quantílica inflated modified kumaraswamy reparametrizada
  rimk<-function(u,mu,alpha,lambda0)
  {
    beta0_cond<-den.cr<-critical<-c()
    for (i in 1:length(u)){
      if (u[i]<=lambda0[i]){y[i]=0}else{
    if(is.na(mu[i])) mu[i]<-.Machine$double.eps
    if(mu[i]<.Machine$double.eps) mu[i]<-.Machine$double.eps
    if(mu[i]>0.9999999) mu[i]<-0.9999999#1-.Machine$double.eps
    critical[i]<-(alpha-alpha/mu[i])
    if(is.na(critical[i])) critical[i]<- -.Machine$double.eps
    if(is.nan(critical[i])){critical[i]=709.7827}
    if(critical[i] < -36.04365) critical[i] <- -36.04365
    den.cr[i]=log(1-exp(critical[i]))
    if(is.nan(den.cr[i])){den.cr[i]=-36.04365}
    beta0_cond[i]<-log(1-tau)/den.cr[i]

    y[i]<- alpha/(alpha-log(1-(1-(u[i]-lambda0[i])/(1-lambda0[i]))^(1/beta0_cond[i])))
}}
   # print("y");print(y)
    return(y)
  }
  
  loglik <- function(z) 
  {
    beta0 <- z[1]
    beta <- z[2:(k+1)]
    alpha <- z[(k+2)] # alpha parameter
    nu0 <- z[(k+3)]
    nu<-z[(k+4):(k+3+c)]
    #print(z)
    error<-rep(0,n) # E(error)=0 
    eta0<-eta<-rep(NA,n)
   
    critical<-c()
    den.cr<-c()
    for(i in 1:n)
    {
      eta0[i] <- Z[i,1]%*%as.matrix(nu0)+ Z[i,2:ncol(Z)]%*%as.matrix(nu)  
      eta[i] <- X[i,1]*beta0 + X[i,2:ncol(X)]%*%as.matrix(beta)  
    }
    lambda0<-linkinv(eta0)
    mu <- linkinv(eta)
    critical<-c()
    den.cr<-c()
    for(i in 1:n)
    {
      if(is.na(mu[i])) mu[i]<-.Machine$double.eps
      if(mu[i]<.Machine$double.eps) mu[i]<-.Machine$double.eps
      if(mu[i]>(0.9999999)) mu[i]<-0.9999999#1-.Machine$double.eps

      critical[i]<-(alpha-alpha/mu[i])
      if(is.na(critical[i])) critical[i]<- -.Machine$double.eps
      if(is.nan(critical[i])){critical[i]=-36.04365}
      if(critical[i] < -36.04365) critical[i] <- -36.04365
      
      den.cr[i]=log(1-exp(critical[i]))
      if(is.nan(den.cr[i])){den.cr[i]=-36.04365}
    } 
    l=critical.y=critical.ly=c()
    for (i in 1:length(y)){
      if(y[i]==0){l[i]=log(lambda0)}else{
      critical.y[i]=exp(alpha-alpha/y[i])
      critical.ly[i]<-log(1-critical.y[i])
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
      if(is.nan(critical.ly[i])){critical.ly[i]=-36.04365}
      if(critical.ly[i] < -36.04365) critical.ly[i] <- -36.04365 #para exp dar .Machine$double.eps
      
      l[i]=log(1-lambda0)+log(alpha)+alpha-alpha/y[i]+(log(1-tau)/den.cr[i] -1)*critical.ly[i]-2*log(y[i])+log(log(1-tau)/den.cr[i])
    }}
    # print("l");print(l)
    # print("log(alpha)+alpha-alpha/y[i]");print(log(alpha)+alpha-alpha/y[i])
    # print("(log(1-tau)/den.cr[i] -1)*critical.ly[i]");print((log(1-tau)/den.cr[i] -1)*critical.ly[i])
    # print("-2*log(y[i])");print(-2*log(y[i]))
    # print('log(log(1-tau)/den.cr[i])');print(log(log(1-tau)/den.cr[i]))
    # print("l");print(sum(l))
    # ll <- dmk(y, alpha, log(1-tau)/den.cr, lambda0,log = TRUE)#log-density modified kumaraswamy quantile re-parametrization
    # print("ll");print(sum(ll))
    return(sum(l))
  }#fim loglik
  
  score <- function(z) 
  {
    beta0 <- z[1]
    beta <- z[2:(k+1)]
    alpha <- z[(k+2)] # alpha parameter
    nu0 <- z[(k+3)]
    nu<-z[(k+4):(k+3+c)]
    error<-rep(0,n) # E(error)=0 
    eta0<-eta<-rep(NA,n)
    
    critical<-c()
    den.cr<-c()
    for(i in 1:n)
    {
      eta0[i] <- Z[i,1]%*%as.matrix(nu0)+ Z[i,2:ncol(Z)]%*%as.matrix(nu)
      eta[i] <- X[i,1]*beta0 + X[i,2:ncol(X)]%*%as.matrix(beta) 
    }
    lambda0<-linkinv(eta0)
    mu <- linkinv(eta)
    mu.eta <-  link$mu.eta
    
    critical<-c()
    den.cr<-c()
    for(i in 1:n)
    {
      if(is.na(mu[i]))mu[i]<-.Machine$double.eps
      if(mu[i]<.Machine$double.eps) mu[i]<-.Machine$double.eps
      if(mu[i]>(0.9999999)) mu[i]<-0.9999999#1-.Machine$double.eps
      critical[i]<-(alpha-alpha/mu[i])
      if(is.na(critical[i])) critical[i]<- -.Machine$double.eps
      if(is.nan(critical[i])){critical[i]=-36.04365}
      if(critical[i] < -36.04365) critical[i] <- -36.04365
      den.cr[i]=log(1-exp(critical[i]))
      if(is.nan(den.cr[i])){den.cr[i]=-36.04365}
    }
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU
    ###########################################################################################################
    mustar<-den1<-num1<-critical.y<-critical.ly<-c()
    for (i in 1:length(y)){
      if(y[i]==0){mustar[i]=0}else{
        critical.y[i]=exp(alpha-alpha/y[i])
        critical.ly[i]<-log(1-critical.y[i])
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
      if(is.nan(critical.ly[i])){critical.ly[i]=-36.04365}
      if(critical.ly[i] < -36.04365) critical.ly[i] <- -36.04365 #para exp dar .Machine$double.eps
        num1[i]<-exp(alpha)*alpha*(den.cr[i]+log(1-tau)*critical.ly[i])
        den1[i]<-(mu[i]^2)*(exp(alpha/mu[i])-exp(alpha))*((den.cr[i])^2)
        mustar[i]<-num1[i]/den1[i]
        }}

    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU 
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA0
    ###########################################################################################################
    lambda0star<-c()
    for (i in 1:length(y)){
      if(y[i]==0){
    lambda0star[i]<-1/lambda0[i]
      }else{
        lambda0star[i]<-1/(lambda0[i]-1)
      }
    }
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA0 
    
    mT <- diag(mu.eta(eta))
    l0T <- diag(mu.eta(eta0))
    
    deta.dnu0<-deta.dbeta0 <- matrix(rep(NA,n),ncol=1)#intercept
    for(i in 1:n)
    {
      deta.dnu0[i,1]<-deta.dbeta0[i,1] <- X[i,1] 
    }
    
    deta.dbeta <- matrix(rep(NA,n*k),ncol=k)#covariates
    for(i in 1:n)
    {
      for(j in 1:k)
      {
        deta.dbeta[i,j] <- X[i,1+j]
      }
    }
    deta.dnu <- matrix(rep(NA,n*c),ncol=c)#covariates
    for(i in 1:n)
    {
      for(j in 1:c)
      {
        deta.dnu[i,j] <- Z[i,1+j]
      }
    }
    
    yl0star <- matrix((lambda0star),ncol=1)
    ymstar <- matrix((mustar),ncol=1)
 
    Ubeta0 <-  t(deta.dbeta0) %*% mT %*% ymstar
    Ubeta <-  t(deta.dbeta) %*% mT %*% ymstar
    Unu0 <-  t(deta.dnu0) %*% l0T %*% yl0star
    Unu <-  t(deta.dnu) %*% l0T %*% yl0star
    ##############################################################################
    ##### START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO alpha
    Ua<-den3<-den2<-num2.part2<-num2.part1<-den1<-num1<-c()
    for (i in 1:length(y)){
      if(y[i]==0){Ua[i]<-0}else{
    num1[i]<-(mu[i]-1)*exp(alpha)*log(1-tau)*critical.ly[i]
    den1[i]<- mu[i]*(exp(alpha/mu[i])-exp(alpha))*(den.cr[i])^2
    num2.part1[i]<-(y[i]-1)*critical.y[i]
    num2.part2[i]<-(log(1-tau)/den.cr[i]-1)
    den2[i]<-y[i]-y[i]*critical.y[i]
    den3[i]<- mu[i]*(exp(alpha/mu[i])-exp(alpha))*den.cr[i]
    Ua[i]<-num1[i]/den1[i]-num2.part1[i]*num2.part2[i]/den2[i]+((mu[i]-1)*exp(alpha))/den3[i]-1/y[i]
      }}
      Ualpha<-sum(1+1/alpha+Ua)
      
    ##### END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO alpha
    ##############################################################################
  
    return(c(Ubeta0,Ubeta,Ualpha,Unu0,Unu))
  }#end score
  
  # initial values for estimation
  index<-which(y!=0)
  x <- cbind(X[index,])#intercepto, covariáveis

  Ynew = linkfun(y[index])
  ajuste = lm.fit(x, Ynew)

  mqo = c(ajuste$coef)
  for(i in 1:length(mqo))
  {
    if (is.na(mqo[i])){
      mqo[i] <- 0
    }
  }
  library(GenSA)
  on.dmk.alpha<-function(alpha){-sum(log(dmk_alpha(alpha)))}

  gen.semchute <- GenSA(lower = c(.Machine$double.eps),
                        upper = c(100),
                        fn = on.dmk.alpha, control=list(max.time=2))
  alpha<-gen.semchute$par
  reg <- c(mqo, alpha,length(y[which(y==0)])/length(y),rep(0,c)) # initializing the parameter values
  
 # reg <- c(0,rep(0,k), alpha,length(y[which(y==0)])/length(y),rep(0,c)) # initializing the parameter values
  # reg <- c(0,rep(0,k), 0,length(y[which(y==0)])/length(y),rep(0,c)) # initializing the parameter values
  
  #reg=c(0,0,alpha)
 # reg <- c(mqo,0)
 # print(reg)
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
  #opt2 <- optim(reg, loglik, method = "BFGS", hessian = T, control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
  #opt<-opt2
  if (opt$conv != 0)
  {
    warning("FUNCTION DID NOT CONVERGE!")
    z$RMC=1
    return(z)
  }
 
  #library(rootSolve)
  #print("verificando derivadas")
  #print(gradient(score,opt2$par))
  #print(hessian(loglik,opt2$par))
  
  # print(rbind(score(opt$par),gradient(loglik,opt$par)))
  # print(rbind(score(opt2$par),gradient(loglik,opt2$par)))
  
  z$conv <- opt$conv
  coef <- (opt$par)[1:(k+c+3)]
  z$coeff <- coef
  z$loglik <- opt$value
  beta0 <-coef[1]#intercept
  beta <-coef[2:(k+1)] #covariates
  alpha <-coef[(k+2)] # alpha parameter
  nu0<-coef[(k+3)]
  nu<-coef[(k+4):(k+3+c)]
  z$beta0 <- beta0
  z$beta <- beta
  z$alpha <- alpha
  z$nu0<-nu0
  z$nu<-nu
  z$RMC=0
   
  eta0hat <-etahat <- rep(NA,n)
  critical<-c()
  den.cr<-c()
  
  for(i in 1:n)
  {
    eta0hat[i] <- Z[i,1]%*%as.matrix(nu0)+ Z[i,2:ncol(Z)]%*%as.matrix(nu)
    etahat[i] <- X[i,1]*beta0 + X[i,2:ncol(X)]%*%as.matrix(beta) 
  }
  lambda0hat<-linkinv(eta0hat)
  muhat <- linkinv(etahat)
# print("muhat");print(muhat)
# print("alpha");print(alpha)
# print("lambda0hat");print(lambda0hat)
# print("n");print(n)
#print(rimk(u=0.5,mu=muhat,alpha=alpha,lambda0=lambda0hat))
  z$fitted<-ts(rimk(u=rep(0.5,n),mu=muhat,alpha=alpha,lambda0=lambda0hat),start=start(y),frequency=frequency(y)) 
  # print(z$fitted)
  # plot(z$fitted,type="l")
  z$etahat <- etahat
  #continuar aqui 
  obs.inf<-function(y,muhat)
  {
    critical<-c()
    den.cr<-c()
    for(i in 1:n)
    {
      if(is.na(muhat[i])) muhat[i]<-.Machine$double.eps
      if(muhat[i]<.Machine$double.eps) muhat[i]<-.Machine$double.eps
      if(muhat[i]>0.9999999) muhat[i]<-0.9999999#1-.Machine$double.eps
      critical[i]<-(alpha-alpha/muhat[i])
      if(is.na(critical[i])) critical[i]<- -.Machine$double.eps
      if(is.nan(critical[i])){critical[i]=-36.04365}
      if(critical[i] < -36.04365) critical[i] <- -36.04365
      den.cr[i]=log(1-exp(critical[i]))
      if(is.nan(den.cr[i])){den.cr[i]=-36.04365}
    }
    
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU
    ###########################################################################################################
    muhatstar.sec<-denominador<- numerador<-c()
    e.ay=e.am=critical.y=critical.ly<-c()
    #print("muhatstar.sec");print(muhatstar.sec)
    for (i in 1:length(y)){
      if(y[i]==0){muhatstar.sec[i] <-0}else{
      e.ay[i]=exp(alpha/y[i])
      e.am[i]=exp(alpha/muhat[i])
      critical.y[i]=exp(alpha-alpha/y[i])
      critical.ly[i]<-log(1-critical.y[i])
      if(is.infinite(e.ay[i])){e.ay[i]=.Machine$double.xmax}
      if(is.infinite(e.am[i])){e.am[i]=.Machine$double.xmax}
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
      if(is.nan(critical.ly[i])){critical.ly[i]=-36.04365}
      if(critical.ly[i] < -36.04365) critical.ly[i] <- -36.04365 #para exp dar .Machine$double.eps
      numerador[i]<- ( exp(alpha)*alpha* ( den.cr[i]*( (2*muhat[i]*exp(alpha)+ e.am[i]*(alpha-2*muhat[i]))*log(1-tau)*critical.ly[i]+exp(alpha)*alpha ) + (2*muhat[i]*exp(alpha) +  e.am[i]*(alpha-2*muhat[i]))*den.cr[i]*den.cr[i] + 2*exp(alpha)*alpha*log(1-tau)*critical.ly[i] ) ) 
      denominador[i] <-( (muhat[i]^4)*((exp(alpha)- e.am[i])^2) *den.cr[i]*den.cr[i]*den.cr[i] )
      muhatstar.sec[i] <-  numerador[i]/denominador[i]
      }}
   # print("muhatstar.sec");print(muhatstar.sec)
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECTO TO MU   
    
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA0
    ###########################################################################################################
    lambda0hatstar.sec<-c()
    for (i in 1:length(y)){
      if(y[i]==0){
        lambda0hatstar.sec[i]<--(1/lambda0hat[i])^2
      }else{
        lambda0hatstar.sec[i]<--(1/(lambda0hat-1))^2
      }
    }
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA0 
    
    deta.dnu0<-deta.dbeta0 <- matrix(rep(NA,n),ncol=1)#intercept
    for(i in 1:n)
    {
      deta.dnu0[i,1]<-deta.dbeta0[i,1] <- X[i,1] 
    }
    
    deta.dbeta <- matrix(rep(NA,n*k),ncol=k)#covariates
    for(i in 1:n)
    {
      for(j in 1:k)
      {
        deta.dbeta[i,j] <- X[i,1+j] 
      }
    }
    deta.dnu <- matrix(rep(NA,n*c),ncol=c)#covariates
    for(i in 1:n)
    {
      for(j in 1:c)
      {
        deta.dnu[i,j] <- Z[i,1+j]
      }
    }
    
    deta.dnu0nu0<-deta.dbeta0beta0<- matrix(0, ncol=1,nrow=n)
    deta.dbeta0beta<- matrix(0, ncol=k,nrow=n)
    deta.dbetabeta<- array(0,dim=c(k,k,n))
    deta.dnu0nu<-matrix(0, ncol=c,nrow=n)
    deta.dnunu<-  array(0,dim=c(c,c,n))
    for(i in 1:n)
    {
      deta.dnu0nu0[i,]<-deta.dbeta0beta0[i,]<- 0 
      deta.dbeta0beta[i,]<-rep(0,k) 
      deta.dnu0nu[i,]<-rep(0,c)
      for(b in 1:k)
      {
        for(a in 1:k)
        {
          deta.dbetabeta[a,b,i] <- 0 
        }
      }
      for(b in 1:c)
      {
        for(a in 1:c)
        {
          deta.dnunu[a,b,i] <- 0 
        }
      }
    }
    
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO ALPHA (CONFERIDA, IGUAL AO ÚLTIMO TERMO DA HESSIANA se usar mu ao invés de muhat)
    ###########################################################################################################
    Uaa<-num1<-den1<-num2<-den2<-num3<-den3<-num4<-den4<-num5<-den5<-num6<-den6<-num7<-den7<-num8<-den8<-num9<-den9<-c()
    for (i in 1:length(y)){
    if(y[i]==0){Uaa[i]<-0}else{
    num1[i]<-2*((muhat[i]-1)^2)*exp(2*alpha)*log(1-tau)*critical.ly[i]
    den1[i]<-(muhat[i]^2)*((exp(alpha)-e.am[i])^2)*((den.cr[i])^3)
    num2[i]<-((muhat[i]-1)^2)*exp(alpha)*log(1-tau)*critical.ly[i]
    den2[i]<-(muhat[i]^2)*(e.am[i]-exp(alpha))*((den.cr[i])^2)
    num3[i]<-((muhat[i]-1)^2)*exp(2*alpha)*log(1-tau)*critical.ly[i]
    den3[i]<-(muhat[i]^2)*((exp(alpha)-e.am[i])^2)*((den.cr[i])^2)
    num4[i]<-((muhat[i]-1)^2)*exp(2*alpha)
    den4[i]<-den3[i]
    num5[i]<-((muhat[i]-1)^2)*exp(alpha)
    den5[i]<-(muhat[i]^2)*(e.am[i]-exp(alpha))*den.cr[i]
    num6[i]<-num4[i]
    den6[i]<-(muhat[i]^2)*((exp(alpha)-e.am[i])^2)*den.cr[i]
    num7[i]<-((y[i]-1)^2)*critical.y[i]*(log(1-tau)/den.cr[i]-1)
    den7[i]<-(y[i]^2)*(1-critical.y[i])
    num8[i]<-((y[i]-1)^2)*exp(2*alpha-2*alpha/y[i])*(log(1-tau)/den.cr[i]-1)
    den8[i]<-(y[i]^2)*((critical.y[i]-1)^2)
    num9[i]<-2*(muhat[i]-1)*exp(2*alpha)*(y[i]-1)*log(1-tau)
    den9[i]<-muhat[i]*y[i]*(exp(alpha)-e.am[i])*(e.ay[i]-exp(alpha))*((den.cr[i])^2)
    Uaa[i]<-sum(num1[i]/den1[i]+num2[i]/den2[i]+num3[i]/den3[i]+num4[i]/den4[i]+num5[i]/den5[i]+num6[i]/den6[i]-num7[i]/den7[i]-num8[i]/den8[i]+num9[i]/den9[i])
    }}
    Ualphaalpha<-sum(Uaa-1/(alpha^2))
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO ALPHA
    
    ####START DERIVATIVE FROM [DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU] IN RESPECT TO ALPHA
    ###########################################################################################################
    e.a=exp(alpha)
    Umualpha<-numerador<-denominador<-c()
    for (i in 1:length(y)){
      if(y[i]==0){Umualpha[i]<-0}else{
    numerador[i]<-(e.a*(den.cr[i]*(log(1-tau)*(muhat[i]*(-e.a)*alpha*(y[i]-1)*(e.a-e.am[i])-y[i]*(muhat[i]*e.a-e.am[i]*(muhat[i]*alpha+muhat[i]-alpha))*(e.a-e.ay[i])*critical.ly[i])+(muhat[i]-1)*e.a*alpha*y[i]*(e.a-e.ay[i]))+2*(muhat[i]-1)*e.a*alpha*y[i]*(e.a-e.ay[i])*log(1-tau)*critical.ly[i]+y[i]*(muhat[i]*e.a-e.am[i]*(muhat[i]*alpha+muhat[i]-alpha))*(-(e.a-e.ay[i]))*((den.cr[i])^2)))
    denominador[i]<-(muhat[i]^3)*y[i]*((e.a-e.am[i])^2)*(e.a-e.ay[i])*((den.cr[i])^3)
    Umualpha[i]<-numerador[i]/denominador[i]
     }}
    
    #print(Umualpha)
    ########################################################################################################### 
    ####END DERIVATIVE FROM [DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU] IN RESPECT TO ALPHA
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU
    ###########################################################################################################
    
    mustar<-num1<-den1<-c()
    for (i in 1:length(y)){
      if(y[i]==0){mustar[i]<-0}else{
        critical.y[i]=exp(alpha-alpha/y[i])
        critical.ly[i]<-log(1-critical.y[i])
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
      if(is.nan(critical.ly[i])){critical.ly[i]=-36.04365}
      if(critical.ly[i] < -36.04365) critical.ly[i] <- -36.04365 #para exp dar .Machine$double.eps
      num1[i]<-exp(alpha)*alpha*(den.cr[i]+log(1-tau)*critical.ly[i])
      den1[i]<-(muhat[i]^2)*(exp(alpha/muhat[i])-exp(alpha))*((den.cr[i])^2)
      mustar[i]<-num1[i]/den1[i]
    }}
   
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU 
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA0
    ###########################################################################################################
    lambda0star<-c()
    for (i in 1:length(y)){
      if(y[i]==0){
        lambda0star[i]<-1/lambda0hat[i]
      }else{
        lambda0star[i]<-1/(lambda0hat[i]-1)
      }
    }
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO LAMBDA0 
    mT <- diag(mu.eta(etahat))   
    num.mT2=dif2link(muhat)
    mT2=diag(-num.mT2*(mu.eta(etahat)^3))
    mV <- diag(muhatstar.sec)
    mV0 <- diag(mustar)
    
    l0T <- diag(mu.eta(eta0hat))
    num.l0T2=dif2link(lambda0hat)
    l0T2=diag(-num.l0T2*(mu.eta(eta0hat)^3))
    l0V <- diag(lambda0hatstar.sec)
    l0V0 <- diag(lambda0star)
    
    mualpha<-diag(Umualpha)
    vI <- matrix(rep(1,n),ncol=1)
    
    KB0B0 <- -(t(deta.dbeta0)%*%mV%*%(mT^2)%*%deta.dbeta0 + t(deta.dbeta0)%*%mV0%*%mT2%*%deta.dbeta0 + t(vI)%*%mV0%*%mT%*%deta.dbeta0beta0)
    KB0B <- -(t(deta.dbeta0)%*%mV%*%(mT^2)%*%deta.dbeta + t(deta.dbeta0)%*%mV0%*%mT2%*%deta.dbeta + t(vI)%*%mV0%*%mT%*%deta.dbeta0beta)
    KB0alpha <- -t(deta.dbeta0)%*% mualpha %*% mT %*% vI
    KB0nu0<-0
    KB0nu<-0
    KBB0<-t(KB0B)
    KBB=matrix(rep(NA,k*k),ncol=k)
    if(length(KBB)==1){
      KBB <- -(t(deta.dbeta)%*%mV%*%(mT^2)%*%deta.dbeta + t(deta.dbeta)%*%mV0%*%mT2%*%deta.dbeta + t(vI)%*%mV0%*%mT%*%deta.dbetabeta)
    }else{
      for(j in 1:k){
        for(i in 1:k){
          KBB[i,j] <- -(t(deta.dbeta[,i])%*%mV%*%(mT^2)%*%deta.dbeta[,j] + t(deta.dbeta[,i])%*%mV0%*%mT2%*%deta.dbeta[,j] + t(vI)%*%mV0%*%mT%*%as.matrix(deta.dbetabeta[i,j,]))
        }
      }
    }
    KBalpha <- -t(deta.dbeta)%*% mualpha %*% mT %*% vI
    KBnu0<-0
    KBnu<-0
    KalphaB0 <- t(KB0alpha)
    KalphaB <- t(KBalpha)
    Kalphaalpha <- -Ualphaalpha
    Kalphanu0<-0
    Kalphanu<-0
    Knu0B0<-t(KB0nu0)
    Knu0B<-t(KBnu0)
    Knu0alpha<-t(Kalphanu0)
    Knu0nu0<--(t(deta.dnu0)%*%l0V%*%(l0T^2)%*%deta.dnu0 + t(deta.dnu0)%*%l0V0%*%l0T2%*%deta.dnu0 + t(vI)%*%l0V0%*%l0T%*%deta.dnu0nu0)
    Knu0nu<--(t(deta.dnu0)%*%l0V%*%(l0T^2)%*%deta.dnu + t(deta.dnu0)%*%l0V0%*%l0T2%*%deta.dnu + t(vI)%*%l0V0%*%l0T%*%deta.dnu0nu)
    KnuB0<-t(KB0nu)
    KnuB<-t(KBnu)
    Knualpha<-t(Kalphanu)
    Knunu0<-t(Knu0nu)
    Knunu<-matrix(rep(NA,c*c),ncol=c)
    if(length(Knunu)==1){
      Knunu <- -(t(deta.dnu)%*%l0V%*%(l0T^2)%*%deta.dnu + t(deta.dnu)%*%l0V0%*%l0T2%*%deta.dnu + t(vI)%*%l0V0%*%l0T%*%deta.dnunu)
    }else{
      for(j in 1:c){
        for(i in 1:c){
          Knunu[i,j] <- -(t(deta.dnu[,i])%*%l0V%*%(l0T^2)%*%deta.dnu[,j] + t(deta.dnu[,i])%*%l0V0%*%l0T2%*%deta.dnu[,j] + t(vI)%*%l0V0%*%l0T%*%as.matrix(deta.dnunu[i,j,]))
        }
      }
    }
    
    K <- rbind(
      cbind(KB0B0,KB0B,KB0alpha,KB0nu0,KB0nu),
      cbind(KBB0,KBB,KBalpha,KBnu0,KBnu),
      cbind(KalphaB0,KalphaB,Kalphaalpha,Kalphanu0,Kalphanu),
      cbind(Knu0B0,Knu0B,Knu0alpha,Knu0nu0,Knu0nu),
      cbind(KnuB0,KnuB,Knualpha,Knunu0,Knunu)
    )
    return(K)
  }
  K<-obs.inf(y,muhat)
  # print("-hessiana = Matriz de informação observada condicional")
  # print(round(K,4))
  # print(opt$par)
  # print("hessiana numerica")
  # print(round(-opt$hessian,4))
  # print("comparando meu cálculo com hessiana da estimação numérica")
  # print(round((K+opt2$hessian),2))
  # print("comparando meu cálculo com hessiana numérica da estimação analítica")
  # print(round((K+opt$hessian),2))
  # print("soma diferença hessiana otimização numérica")
  # print(round(sum(abs(K+opt2$hessian)),2))
  # print("soma diferença hessiana numérica otimização analítica")
  # print(round(sum(abs(K+opt$hessian)),2))
  # print(solve(K))
  Ksolve<- tryCatch(solve(K), error = function(e) return("error"))
  if(Ksolve[1] == "error")
  {z$RMC=1#used at Monte-Carlo simulation for discard from the sample
  warning("Analytic Observed Information Matrix is not positive semi-definite")
  return(z)#if Analytic Observed Information Matrix is not positive semi-definite, do not calculate
  }else{sol=try(solve(K))}

  v<-diag(sol)#Variância assintótica dos esimadores
  # print(v)
  for (i in 1:length(v))
  {
    if(is.na(v[i]) | is.nan(v[i]) | v[i]<0 )  {
      z$RMC=1
      warning("Analytic Observed Information Matrix is not positive semi-definite")
      return(z)#if Analytic Observed Information Matrix is not positive semi-definite, do not calculate
    }
  }
  #return(z) # # # # # # # # # # # # # # tirar depois, abaixo nao precisa na simulacao # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  z$zstat<-z$coeff/sqrt(v)
  resp<-rep(0,length(z$zstat))
  for (i in 1:length(resp)){
    if(abs(z$zstat[i])>qnorm(0.975))
    {
      resp[i] <- "H0 rejected"
    } else {resp[i] <- "H0 not rejected"}
  }
  LI<-z$coeff-qnorm(0.975)*sqrt(v)
  LS<-z$coeff+qnorm(0.975)*sqrt(v)
  z$pvalues<-(1-pnorm(abs(z$zstat)))*2
  first_col<-c("intercept", 1:k,"alpha estimator","nu0",1:c)
  result <- matrix(c(first_col,round(c(z$coeff,z$zstat,LI,LS,z$pvalues),4),resp), nrow=length(z$coeff), ncol=7, byrow=F)
  colnames(result) <- c("Estimator","MLE","Wald's Statistic","Lower bound","Upper bound","p-value","Wald'S Test result")
  rownames(result)<-c("", rep("cov X",k),"","", rep("cov Z",c))
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
  
  ytofit<-ts(c(y[1:n]),start=start(y),frequency=frequency(y))
  
  
  ###########################
  
  z$serie <- y
  #quantile residuals
  
  critical<-c()
  den.cr<-c()
  for(i in 1:n)
  {
    if(is.na(muhat[i])) muhat[i]<-.Machine$double.eps
    if(muhat[i]<.Machine$double.eps) muhat[i]<-.Machine$double.eps
    if(muhat[i]>0.9999999) muhat[i]<-0.9999999#1-.Machine$double.eps
    critical[i]<-(alpha-alpha/muhat[i])
    if(is.na(critical[i])) critical[i]<- -.Machine$double.eps
    if(is.nan(critical[i])){critical[i]=-36.04365}
    if(critical[i] < -36.04365) critical[i] <- -36.04365
    den.cr[i]=log(1-exp(critical[i]))
    if(is.nan(den.cr[i])){den.cr[i]=-36.04365}
   }
  z$resid1 <- as.vector(qnorm(pmk(y[1:n],alpha, log(1-tau)/den.cr[1:n],lambda0=lambda0hat,log.p = FALSE ) ))
  #deviance residuals
  l_tilde <- (dmk(y[1:n], alpha, log(1-tau)/log(1-exp(alpha-alpha/y[1:n])),lambda0=lambda0hat, log = TRUE))#y[1:n] where was mu
  l_hat <- (dmk(y[1:n], alpha, log(1-tau)/den.cr[1:n],lambda0=lambda0hat, log = TRUE))#muhat where was mu
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
 # print(z$mresult)
  

  if(resid==1) {
    residual <- z$resid1
  }
  
  if(resid==2) {
    residual <- z$resid2
  }
  
  z$residual<-residual
 # print(z$residual)
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
    pdf("resid_v_ind.pdf",width=5, height=4)
    {
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1))
      par(mgp=c(1.7, 0.45, 0))
      plot(residual,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
      lines(t,rep(-3,n+12)#length(residual))
            ,lty=2,col=1)
      lines(t,rep(3,n+12)#length(residual))
            ,lty=2,col=1)
      lines(t,rep(-2,n+12)#length(residual))
            ,lty=3,col=1)
      lines(t,rep(2,n+12)#length(residual))
            ,lty=3,col=1)
    }
    dev.off()
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
    pdf("envelope_plot.pdf",width=5, height=4)
    { 
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
    library(hnp)
    hnp(residual,half=F, sim = 1000, conf = 0.95, ylab="Empirical quantile", xlab="Normal quantile")
    }
    dev.off()
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
  #Multiple R-squared
  #Adjusted R-squared
  
  # loglik_null <- function(z) 
  # {
  #   beta0 <- z[1]
  #   beta <- z[2:(k+1)]
  #   alpha <- z[(k+2)]
  #   
  #   l=0
  #     return(l)
  # }
  # ini_null<- c(mean(z$beta0),mean(z$beta),0)
  # 
  # opt_null <- optim(ini_null, loglik_null,method = "BFGS", control = list(fnscale = -1)) # , maxit = 500, reltol = 1e-9))
  # r2 <- 1-exp((-2/n)*(opt$value-opt_null$value))
  # z$r2 <- r
  r2 <- 1-exp((-2/n)*(opt$value-sum(l_tilde)))
  z$r2 <- r2
 # print(z$r2)
#print("aqui")
  #null hypothesis: normality
  library(nortest)
  andersondarling=ad.test(residual)
  z$andersondarling<-andersondarling$statistic
  z$p_andersondarling<-andersondarling$p.value
  
  shapiro=shapiro.test(residual)
  z$shapiro=shapiro$statistic
  z$p_shapiro=shapiro$p.value
 # print(z$p_andersondarling)
  
  #non autocorrelated
  ljungbox<- Box.test(residual, lag = 10, type = "Ljung-Box", fitdf = k)
  z$ljungbox<-ljungbox$statistic
    z$p_ljungbox<-ljungbox$p.value
    
    
    dw<-function(res){
      alternative = "two.sided"
      dw <- sum(diff(res)^2)/sum(res^2)
      Q1 <- chol2inv(qr.R(qr(exvar.beta)))
      if (n<100) {
        A <- diag(c(1, rep(2, n - 2), 1))
        A[abs(row(A) - col(A)) == 1] <- -1
        MA <- diag(rep(1, n)) - X %*% Q1 %*% t(exvar.beta)
        MA <- MA %*% A
        ev <- eigen(MA)$values[1:(n - k)]
        if (any(Im(ev) > tol)) 
          warning("imaginary parts of eigenvalues discarded")
        ev <- Re(ev)
        ev <- ev[ev > tol]
        pdw <- function(dw) .Fortran("pan", as.double(c(dw, 
                                                        ev)), as.integer(length(ev)), as.double(0), 
                                     as.integer(iterations), x = double(1), PACKAGE = "lmtest")$x
        pval <- switch(alternative, two.sided = (2 * min(pdw(dw), 
                                                         1 - pdw(dw))), less = (1 - pdw(dw)), greater = pdw(dw))
        if (is.na(pval) || ((pval > 1) | (pval < 0))) {
          warning("exact p value cannot be computed (not in [0,1]), approximate p value will be used")
          exact <- FALSE
        }
      }else{
          if (n < max(5, k)) {
            warning("not enough observations for computing an approximate p value, set to 1")
            pval <- 1
          }
          else {
            AX <- matrix(as.vector(filter(exvar.beta, c(-1, 2, -1))), 
                         ncol = k)
            AX[1, ] <- exvar.beta[1, ] - exvar.beta[2, ]
            AX[n, ] <- exvar.beta[n, ] - exvar.beta[(n - 1), ]
            XAXQ <- t(exvar.beta) %*% AX %*% Q1
            P <- 2 * (n - 1) - sum(diag(XAXQ))
            Q <- 2 * (3 * n - 4) - 2 * sum(diag(crossprod(AX) %*% 
                                                  Q1)) + sum(diag(XAXQ %*% XAXQ))
            dmean <- P/(n - k)
            dvar <- 2/((n - k) * (n - k + 2)) * (Q - P * 
                                                   dmean)
            pval <- switch(alternative, two.sided = (2 * 
                                                       pnorm(abs(dw - dmean), sd = sqrt(dvar), lower.tail = FALSE)), 
                           less = pnorm(dw, mean = dmean, sd = sqrt(dvar), 
                                        lower.tail = FALSE), greater = pnorm(dw, 
                                                                             mean = dmean, sd = sqrt(dvar)))
          }
      }
      alternative <- switch(alternative, two.sided = "true autocorrelation is not 0", 
                            less = "true autocorrelation is less than 0", greater = "true autocorrelation is greater than 0")
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
 # variance inflation factor (VIF) for each independent variable, and a VIF value greater than 1.5 indicates multicollinearity.
  #library(car) 
  #vif()
  
  #library(caTools)
  #library(quantmod)
  
 # Multicollinearity test can be checked by
  z$cor.beta <- cor(exvar.beta)                                         # independent variables correlation matrix 
  z$cor.nu <- cor(exvar.nu)
  #null hypothesis: non-heteroscedasticity (constant variance)
  #https://cran.r-project.org/web/packages/olsrr/vignettes/heteroskedasticity.html
  # Heteroscedasticity
  SBP<-function(resi){#adapted to fits a linear regression model to the residuals of the mkreg model 
    sigma2 <- sum(resi^2)/n
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

  diagnostic<-matrix(round(c(z$andersondarling,z$shapiro,z$durbin,z$ljungbox,z$breusch,
                             z$p_andersondarling,z$p_shapiro,z$p_durbin,z$p_ljungbox,z$p_breusch
  ),4), nrow=2, ncol=5, byrow=T)
  colnames(diagnostic) <- c("Anderson-Darling test","Shapiro-Wilk test","Durbin-Watson test","Ljung-Box test","studentized Breusch-Pagan test")
  rownames(diagnostic) <- c("Statistic","P-value")
  z$diagnostic=diagnostic
  
  if(print==T){
    print("MKreg",quote=F)
  print(z$coef.result,quote=F)
  message("")
  print(c("Log-likelihood =",round(z$loglik,4)),quote=F)
  print(c("aic =",round(z$aic,4),"bic =",round(z$bic,4)),quote=F)
  print(c("Deviance =",round(z$deviance,4)," DF:",z$dof.dev),quote=F)
  print(c("R-squared=",round(z$r2,4)),quote=F)
  message("")  
  if(resid==1) {
    print("Quantile residuals:",quote=F)
  }
  
  if(resid==2) {
    print("Deviance residuals:",quote=F)
  }
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
