sample.mkreg <- function(n,exvar.beta=NA,exvar.nu=NA,exvar.rho=NA,beta=c(0.0),nu=c(0.0),rho=c(0.0),alpha=10,tau=0.5,link="logit")
{
  if(is.na(beta[1])==F & beta[1]!=0 & is.na(nu[1])==F & nu[1]!=0 & is.na(rho[1])==F & rho[1]!=0){
    source("sample.imkreg0Z1A.R")
    sample.mkreg(n=n,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,beta=beta,nu=nu,rho=rho,alpha=alpha,tau=tau,link=link)
  }else{
    if(is.na(beta[1])==F & beta[1]!=0 & is.na(nu[1])==F & nu[1]!=0 & any(is.na(rho[1])==T | rho[1]==0)){
      source("sample.imkregX0.R")
      sample.mkreg(n=n,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,beta=beta,nu=nu,rho=rho,alpha=alpha,tau=tau,link=link)
    }else{
    if(is.na(beta[1])==F & beta[1]!=0 & is.na(nu[1])==F & nu[1]!=0 & any(is.na(rho[1])==T | rho[1]==0)){
      parameters=c(beta,alpha,nu)
    }else{
      if(is.na(beta[1])==F & beta[1]!=0 & any(is.na(nu[1])==T | nu==0) & is.na(rho[1])==F & rho[1]!=0){
        parameters=c(beta,alpha,rho)
      }else{
        source("sample.mkreg.R")
        sample.mkreg(n=n,beta=beta,exvar.beta=exvar.beta,exvar.nu=NA,exvar.rho=NA,nu=NA, rho=NA,alpha=alpha,tau=tau,link=link)
        }
    }}}
}
  