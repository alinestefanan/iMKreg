#GENERATE iMKreg SAMPLE 
source("sample.R")
n=50
x2<-c(rep(0,2*n%/%3),rep(1,n-2*n%/%3))
x3<-c(rep(1,n-2*n%/%3),rep(0,2*n%/%3))
exvar.beta<-cbind(x2,x3)
exvar.nu<-NA
exvar.rho<-NA
y<-sample.mkreg(n,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho,beta=c(0.6812028,0.1722346,-0.7170177),nu=c(-2.6389866),rho=NA,alpha=11.3205890)

#iMKreg FIT
source("fit.R")
iMK<-mkreg(y,exvar.beta=exvar.beta,exvar.nu=exvar.nu,exvar.rho=exvar.rho)


