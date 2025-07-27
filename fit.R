mkreg <- function(y,exvar.beta=NA,exvar.nu=NA,exvar.rho=NA,tau=0.5,graph=T,print=T,check=F,link="logit")
{
  library(lmtest)
  zero<-which(y==0)
  one<-which(y==1)
  zero.one<-c(zero,one)
  if(length(zero)!=0 & length(zero.one)==length(zero)){
      if(is.na(exvar.beta[1])==F  & is.na(exvar.nu[1])==F  & is.na(exvar.rho[1])==T){
        source("imkregX0Z.R")
        imkreg0Z(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
      }else if(is.na(exvar.beta[1])==F  & is.na(exvar.nu[1])==T & is.na(exvar.rho[1])==T ){
        source("imkregX0.R")
        imkreg0(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
      }else if(is.na(exvar.beta[1])==T  & is.na(exvar.nu[1])==T & is.na(exvar.rho[1])==T ){
        source("imkreg0.R")
        imkreg0(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
      }
  }else if(length(one)!=0 & length(zero.one)==length(one)){
    if(is.na(exvar.beta[1])==F  & is.na(exvar.nu[1])==T & is.na(exvar.rho[1])==F ){
      source("imkregX1A.R")
      imkreg1A(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
    }else if(is.na(exvar.beta[1])==F  & is.na(exvar.nu[1])==T & any(is.na(exvar.rho[1])==T )){
      source("imkregX1.R")
      imkreg1(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
    }else if(is.na(exvar.beta[1])==T  & is.na(exvar.nu[1])==T & any(is.na(exvar.rho[1])==T )){
      source("imkreg1.R")
      imkreg1(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
    }
  }else if(length(zero.one)!=0){
    if(is.na(exvar.beta[1])==T  & is.na(exvar.nu[1])==T  & is.na(exvar.rho[1])==T ){
      source("imkreg01.R")
      imkreg01(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
    }else if(is.na(exvar.beta[1])==F  & is.na(exvar.nu[1])==F  & is.na(exvar.rho[1])==F ){
      source("imkregX0Z1A.R")
      imkreg0Z1A(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
    }else if(is.na(exvar.beta[1])==F  & is.na(exvar.nu[1])==F  & is.na(exvar.rho[1])==T ){
      source("imkregX0Z1.R")
      imkreg0Z1(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
    }else if(is.na(exvar.beta[1])==F  & is.na(exvar.nu[1])==T & is.na(exvar.rho[1])==F ){
      source("imkregX01A.R")
      imkreg01A(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
    }else if(is.na(exvar.beta[1])==F  & is.na(exvar.nu[1])==T  & is.na(exvar.rho[1])==T ){
      source("imkregX01.R")
      imkreg01(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
    }else if(is.na(exvar.beta[1])==T  & is.na(exvar.nu[1])==F  & is.na(exvar.rho[1])==T ){
      source("imkreg0Z1.R")
      imkreg0Z1(y=y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,tau=tau,graph=graph,print=print,check=check,link=link)
    }
  }else{
    source("mkreg.R")
    mkreg01(y=y,exvar.beta=exvar.beta,exvar.nu=NA,exvar.rho=NA,tau=tau,graph=graph,print=print,check=check,link=link)
  }
}