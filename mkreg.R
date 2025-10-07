# exvar.beta: covariate column matrix
# tau: quantil, when set 0.5 is the median
# resid: 1 = quantile residual, 2 = deviance residual
# link: "logit", "probit" or "cloglog"

mkreg01 <- function(y,exvar.beta=NA,exvar.nu=NA,exvar.rho=NA,tau=0.5,resid=1,print=T,check=F,link="logit")
{
  n <- length(y) 
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
    l=log(alpha)+alpha-alpha/y+(log(1-tau)/den.cr -1)*critical.ly-2*log(y)+log(log(1-tau)/den.cr)
    # print("alpha");print(alpha)
    # print("l");print(sum(l))
    # ll <- dmk(y, alpha, log(1-tau)/den.cr, log = TRUE)#log-density modified kumaraswamy quantile re-parametrization
    # print("ll");print(sum(ll))
    return(sum(l))
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
  
  obs.inf<-function(y,muhat)
  {
    # muhat[is.na(muhat)]<-.Machine$double.eps
    # muhat[muhat<.Machine$double.eps]<-.Machine$double.eps
    # muhat[muhat>0.9999999]<-0.9999999
    critical<-(alpha-alpha/muhat)
    # critical[is.na(critical)]<--.Machine$double.eps
    # critical[is.nan(critical)]<--36.04365
    # critical[critical< (-36.04365)]<--36.04365
    den.cr=log(1-exp(critical))
    den.cr[is.nan(den.cr)]<--36.04365
    e.ay<-exp(alpha/y)
    # e.ay[is.infinite(e.ay)]<-.Machine$double.xmax
    e.am<-exp(alpha/muhat)
    # e.am[is.infinite(e.am)]<-.Machine$double.xmax
    critical.y=exp(alpha-alpha/y)
    # critical.y[is.infinite(critical.y)]<-.Machine$double.xmax
    critical.ly<-log(1-critical.y)
    critical.ly[is.nan(critical.ly)]<--36.04365
    # critical.ly[critical.ly< (-36.04365)]<--36.04365
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU
    ###########################################################################################################
    numerador<- ( exp(alpha)*alpha* ( den.cr*( (2*muhat*exp(alpha)+ e.am*(alpha-2*muhat))*log(1-tau)*critical.ly+exp(alpha)*alpha ) + (2*muhat*exp(alpha) +  e.am*(alpha-2*muhat))*den.cr*den.cr + 2*exp(alpha)*alpha*log(1-tau)*critical.ly ) ) 
    denominador <-( (muhat^4)*((exp(alpha)- e.am)^2) *den.cr*den.cr*den.cr )
    muhatstar.sec <-  numerador/denominador
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECTO TO MU   
    
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
    
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO ALPHA (CONFERIDA, IGUAL AO ÚLTIMO TERMO DA HESSIANA se usar mu ao invés de muhat)
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
    
    Ualphaalpha<-sum(num1/den1+num2/den2+num3/den3+num4/den4+num5/den5+num6/den6-num7/den7-num8/den8+num9/den9-1/(alpha^2))
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO ALPHA
    
    ####START DERIVATIVE FROM [DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU] IN RESPECT TO ALPHA
    ###########################################################################################################
    e.a=exp(alpha)
    numerador<-(e.a*(den.cr*(log(1-tau)*(muhat*(-e.a)*alpha*(y-1)*(e.a-e.am)-y*(muhat*e.a-e.am*(muhat*alpha+muhat-alpha))*(e.a-e.ay)*critical.ly)+(muhat-1)*e.a*alpha*y*(e.a-e.ay))+2*(muhat-1)*e.a*alpha*y*(e.a-e.ay)*log(1-tau)*critical.ly+y*(muhat*e.a-e.am*(muhat*alpha+muhat-alpha))*(-(e.a-e.ay))*((den.cr)^2)))
    denominador<-(muhat^3)*y*((e.a-e.am)^2)*(e.a-e.ay)*((den.cr)^3)
    Umualpha<-numerador/denominador
    ########################################################################################################### 
    ####END DERIVATIVE FROM [DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU] IN RESPECT TO ALPHA
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU
    ###########################################################################################################
    num1<-exp(alpha)*alpha*(den.cr+log(1-tau)*critical.ly)
    den1<-(muhat^2)*(exp(alpha/muhat)-exp(alpha))*((den.cr)^2)
    mustar<-num1/den1
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU 
    
    mT <- diag(as.vector(mu.eta(X%*%as.matrix(beta)))) 
    num.mT2=dif2link(muhat)
    mT2=diag(as.vector(-num.mT2*(mu.eta(X%*%as.matrix(beta))^3)))
    mV <- diag(as.vector(muhatstar.sec))
    mV0 <- diag(as.vector(mustar))
    mualpha<-diag(as.vector(Umualpha))
    vI <- matrix(rep(1,n),ncol=1)
    
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
    
    
    KalphaB <- t(KBalpha)
    Kalphaalpha <- -Ualphaalpha
    
    K <- rbind(
      cbind(KBB,KBalpha),
      cbind(KalphaB,Kalphaalpha)
    )
    return(K)
  }
  K<-obs.inf(y,muhat)
  
  Ksolve<- tryCatch(solve(K), error = function(e) return("error"))
  if(Ksolve[1] == "error")
  {z$RMC=1#used at Monte-Carlo simulation for discard from the sample
  warning("Analytic Observed Information Matrix is not positive semi-definite")
  return(z)#if Analytic Observed Information Matrix is not positive semi-definite, do not calculate
  }else{sol=try(solve(K))}
  
  v<-diag(sol)#Variância assintótica dos esimadores
  
  for (i in 1:length(v))
  {
    if(is.na(v[i]) | is.nan(v[i]) | v[i]<0 )  {
      z$RMC=1
      warning("Analytic Observed Information Matrix is not positive semi-definite")
      return(z)#if Analytic Observed Information Matrix is not positive semi-definite, do not calculate
    }
  }
  
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
  z$LI<-LI
  z$LS<-LS
  z$pvalues<-(1-pnorm(abs(z$zstat)))*2
  first_col<-c(0:k,"alpha")
  result <- matrix(c(first_col,round(c(z$coeff,z$zstat,LI,LS,z$pvalues),4),resp), nrow=length(z$coeff), ncol=7, byrow=F)
  colnames(result) <- c("Parameter","MLE","Wald's Statistic","Lower bound","Upper bound","p-value","Wald'S Test result")
  rownames(result)<-c(rep("beta",(k+1)),"")
  z$coef.result<-result
  z$aic <- -2*(z$loglik)+2*(length(opt$par)) 
  z$bic <- -2*(z$loglik)+(length(opt$par))*log(n)
  result2<-matrix(round(c(z$loglik,z$aic,z$bic),4),nrow=3,ncol=1)
  rownames(result2)<-c("Log-likelihood","AIC","BIC")
  
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
   
  ########################################################################
  ########################   residual analysis   ########################
  ########################################################################
  
  #null hypothesis: normality
  library(nortest)
  andersondarling=ad.test(residual)
  z$andersondarling<-andersondarling$statistic
  z$p_andersondarling<-andersondarling$p.value
  
  if(length(residual)<=5000){
  shapiro=shapiro.test(residual)
  z$shapiro=shapiro$statistic
  z$p_shapiro=shapiro$p.value
  }else{
    shapiro=shapiro.test(residual[1:5000])
    z$shapiro=shapiro$statistic
    z$p_shapiro=shapiro$p.value
  }
  
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
      MA <- diag(rep(1, n)) - exvar.beta %*% Q1 %*% t(exvar.beta)
      MA <- MA %*% A
      ev <- eigen(MA)$values[1:(n - k)]
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
      if (n < max(5, k)) {
        warning("not enough observations for computing an approximate p value, set to 1")
        pval <- 1
      }
      else {
        AX <- matrix(as.vector(filter(exvar.beta, c(-1, 2, -1))), ncol = k)
        AX[1, ] <- exvar.beta[1, ] - exvar.beta[2, ]
        AX[n, ] <- exvar.beta[n, ] - exvar.beta[(n - 1), ]
        XAXQ <- t(exvar.beta) %*% AX %*% Q1
        P <- 2 * (n - 1) - sum(diag(XAXQ))
        Q <- 2 * (3 * n - 4) - 2 * sum(diag(crossprod(AX) %*% Q1)) + sum(diag(XAXQ %*% XAXQ))
        dmean <- P/(n - k)
        dvar <- 2/((n - k) * (n - k + 2)) * (Q - P * dmean)
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

  
  #  non-multicollinearity test can be checked by
  z$var <- cor(exvar.beta)                                         # independent variables correlation matrix 
  
  #null hypothesis: non-heteroscedasticity (constant variance)
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

