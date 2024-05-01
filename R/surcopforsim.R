
# This is an  function named 'Bivariate Copula based Seemingly unrelated Regression'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#library(copula)
#library(VineCopula)
#library(fGarch)

# Support funtions

normallike<-function(param,yyy,Xmat){
#### yyy is the response, Xmat is the covariates
   H=length(param)
   beta=param[1:(H-1)]
   sigma=param[H]

   mu = Xmat%*%beta
   e=yyy-mu
   llk <- dnorm( e, mean=0,sd = abs(sigma),log=TRUE)

   u = pnorm(e, mean=0,sd = abs(sigma))

   u<-ifelse (u > 0.999999999999999, 0.999999999999999, u)
   Results <- list(llk, u)
   return(Results)
   }

tlike<-function(param,yyy,Xmat){
   H=length(param)
   v1=param[H]
   beta=param[1:(H-2)]
   sigma=param[H-1]
   mu = Xmat%*%beta
   e=yyy-mu
llk <- dstd(e, mean=0,nu=v1,sd = abs(sigma),log=TRUE)  # student t       #2) ### Adapt likelihood here####
# student t marginal
   u = pstd(e,mean=0,sd= abs(sigma),nu=v1)

   u<-ifelse (u > 0.999999999999999, 0.999999999999999, u)
   Results <- list(llk, u)
   return(Results)
   }


######################################################################
# 3) DEFINE THE COPULA LIKELIHOOD

loglikCopula <-function(theta, y1, xmat1, y2, xmat2,Copula,type){

# 1 = G-G
# 2 = G-T
# 3 = T-G
# 4 = T-T
  Xmat1=cbind(1,xmat1)
  Xmat2=cbind(1,xmat2)
if (type==1) {
    F=length(theta)
    K1=ncol(Xmat1)+1
    K2=ncol(Xmat2)+1
    param1 = theta[1:K1]
    param2 = theta[(K1+1):(F-1)]
    rho=theta[F]}
if (type==2) {
    F=length(theta)
    K1=ncol(Xmat1)+1
    K2=ncol(Xmat2)+2
    param1 = theta[1:K1]
    param2 = theta[(K1+1):(F-1)]
    rho=theta[F]}
if (type==3) {
    F=length(theta)
    K1=ncol(Xmat1)+2
    K2=ncol(Xmat2)+1
    param1 = theta[1:K1]
    param2 = theta[(K1+1):(F-1)]
    rho=theta[F]}
if (type==4) {
    F=length(theta)
    K1=ncol(Xmat1)+2
    K2=ncol(Xmat2)+2
    param1 = theta[1:K1]
    param2 = theta[(K1+1):(F-1)]
    rho=theta[F]}


if (type==1) {
    temp1 <-  normallike( param1,yyy=y1,Xmat=Xmat1)
    temp2 <-  normallike( param2,yyy=y2,Xmat=Xmat2) }

if (type==2) {
    temp1 <-  normallike( param1,yyy=y1,Xmat=Xmat1)
    temp2 <-  tlike( param2,yyy=y2,Xmat=Xmat2) }

if (type==3) {
    temp1 <-  tlike( param1,yyy=y1,Xmat=Xmat1)
    temp2 <-  normallike( param2,yyy=y2,Xmat=Xmat2) }

if (type==4) {
    temp1 <-  tlike( param1,yyy=y1,Xmat=Xmat1)
    temp2 <-  tlike( param2,yyy=y2,Xmat=Xmat2) }



      uu = cbind(temp1[[2]],temp2[[2]])
      demandloglik = temp1[[1]]
      supplyloglik = temp2[[1]]
      u=uu[,1]
      v=uu[,2]
if (Copula=="frank") {
rcltcopula=BiCopPDF(u, v, family=5, par=rho, par2=0)+0.000001}

if (Copula=="joe"){
rcltcopula=BiCopPDF(u, v, family=6, par=rho, par2=0)+0.000001}

if (Copula=="clayton"){
rcltcopula=BiCopPDF(u, v, family=3, par=rho, par2=0)+0.000001}

if (Copula=="clayton90"){
rcltcopula=BiCopPDF(u, v, family=23, par=rho, par2=0)+0.000001}

if (Copula=="gumbel"){
rcltcopula=BiCopPDF(u, v, family=4, par=rho, par2=0)+0.00000001}

if (Copula=="gumbel90"){
rcltcopula=BiCopPDF(u, v, family=24, par=rho, par2=0)+0.0000001}

if (Copula=="studentt"){
rcltcopula=BiCopPDF(u, v, family=2, par=rho, par2=4)+0.000000001 }

if (Copula=="normal"){
rcltcopula=BiCopPDF(u, v, family=1, par=rho, par2=0)+0.00000001}


if (Copula=="clayton270"){
rcltcopula=BiCopPDF(u, v, family=33, par=rho, par2=0)+0.000000001}

if (Copula=="gumbel270"){
rcltcopula=BiCopPDF(u, v, family=34, par=rho, par2=0)+0.000000001}



      Loglik <- sum(demandloglik) + sum(supplyloglik) + sum(log(rcltcopula))

   LL=(-Loglik )
   if (is.infinite(LL))  # control for optimization
  LL=-1000000
  if (is.nan(LL))  # control for optimization
   LL=-1000000

   return(LL)
   }
#### 4) estimate system linear copula #################################


BicopSUR=function (Y1,X1,Y2,X2,Copula,type,rho,Lcop,Ucop){

Eq1=lm(Y1~X1)
Eq2=lm(Y2~X2)
sigma=1
v=5
beta1init=c(coef(Eq1),sigma1=sigma)
beta2init=c(coef(Eq2),sigma2=sigma)
K1 = length(coef(Eq1))
K2 = length(coef(Eq2))
if (type == 1) {
  theta = c(beta1init, beta2init, rho = rho)
  lower = c(rep(-Inf, K1), 1e-04, rep(-Inf, K2), 0.001, Lcop)
  upper = c(rep(Inf, K1), Inf, rep(Inf, K2), Inf, Ucop)
}
if (type == 2) {
  theta = c(beta1init, beta2init, df = v, rho = rho)
  lower = c(rep(-Inf, K1), 1e-04, rep(-Inf, K2), 0.001,df=2.1, Lcop)
  upper = c(rep(Inf, K1), Inf, rep(Inf, K2), Inf, 30,Ucop)
}
if (type == 3) {
  theta = c(beta1init, df = v, beta2init, rho = rho)
  lower = c(rep(-Inf, K1), 1e-04,df=2.1, rep(-Inf, K2), 0.001, Lcop)
  upper = c(rep(Inf, K1), Inf,30, rep(Inf, K2), Inf, Ucop)
}
if (type == 4) {
  theta = c(beta1init, df = v, beta2init, df = v, rho = rho)
  lower = c(rep(-Inf, K1), 1e-04,df=2.1, rep(-Inf, K2), 0.001, df=2.1,Lcop)
  upper = c(rep(Inf, K1), Inf,30, rep(Inf, K2), Inf,df=30, Ucop)
}
=======
 v=5
beta1init = c(coef(Eq1), sigma1 = sigma)
beta2init = c(coef(Eq2), sigma2 = sigma)

 K1 = length(coef(Eq1))
    K2 = length(coef(Eq2))
    if (type == 1) {
        theta = c(beta1init, beta2init, rho = rho)
    lower = c(rep(-Inf, K1), 1e-04, rep(-Inf, K2), 0.001, Lcop)
    upper = c(rep(Inf, K1), Inf, rep(Inf, K2), Inf, Ucop)
    }
    if (type == 2) {
        theta = c(beta1init, beta2init, df = v, rho = rho)
    lower = c(rep(-Inf, K1), 1e-04, rep(-Inf, K2), 0.001,df=2.1, Lcop)
    upper = c(rep(Inf, K1), Inf, rep(Inf, K2), Inf, 30,Ucop)
    }
    if (type == 3) {
        theta = c(beta1init, df = v, beta2init, rho = rho)
    lower = c(rep(-Inf, K1), 1e-04,df=2.1, rep(-Inf, K2), 0.001, Lcop)
    upper = c(rep(Inf, K1), Inf,30, rep(Inf, K2), Inf, Ucop)
    }
    if (type == 4) {
        theta = c(beta1init, df = v, beta2init, df = v, rho = rho)
    lower = c(rep(-Inf, K1), 1e-04,df=2.1, rep(-Inf, K2), 0.001, df=2.1,Lcop)
    upper = c(rep(Inf, K1), Inf,30, rep(Inf, K2), Inf,df=30, Ucop)
    }


## Step 3 ) change Copula, Type
model <- optim(theta,loglikCopula, y1=Y1,xmat1=X1,
              y2=Y2,xmat2=X2,Copula=Copula,type=type,
              method=c("L-BFGS-B"),lower =lower,upper = upper,
              control=list(trace=1,maxit=100000),hessian=TRUE )

# goodman(1996)
MBF2=function(t){
if (length(t)==1){
mbf=exp(-0.5*(t^2))
}

if (length(t)>1){
mbf=c()
for ( i in 1:length(t)){
mm=exp(-0.5*(t[i]^2))
mbf[i]=round(mm,4)}
}
mbf
}

# table of results
n=length(Y1)
coef<- model$par
model$se <- sqrt(diag(solve(-model$hessian)))
(paramsWithTs = cbind (model$par , model$par/model$se ) )

for ( i in 1:length(coef)){
if (is.nan(model$se[i]))
model$se[i]=sqrt(diag(solve(model$hessian)))[i]
}
S.E.=model$se
stat=coef/S.E.
mbf=MBF2(stat)
pvalue <- 2*(1 - pnorm(abs(stat)))
result <- cbind(coef,S.E.,stat,pvalue,mbf)
BIC= -2*model$value+ (log(n)*length(coef))
AIC = -2*model$value + 2*length(coef)

output=list(
  result=result,
  AIC=AIC,
  BIC=BIC
)
output
}
#########################################################
## SUR estimation
# library(systemfit)
#fitols <- systemfit( system, data=Kmenta,method = "SUR" )
#print( fitols )


#normal (-.99, .99)
#studentt (-.99, .99)
#clayton (0.1, Inf)
#gumbel [1.1, Inf)
#frank {Inf}
#joe (1.1, Inf)
#clayton270  (90 and 270 degrees) (-Inf, -0.01)
#gumbel270 (90 and 270 degrees) (-Inf, -1.1]

# Simulate Gaussian Copula
#n=1000
#sim = BiCopSim(n,family=1,0.5)
#e1=qnorm(sim[,1])
#e2=qnorm(sim[,2])
#X1=rnorm(n)
#X2=rnorm(n)
#Y1=1+2*X1 +e1
#Y2=1+2*X2 +e2
#type=1
#Copula="normal"
#rho=0.5
#Lcop=-0.99
#Ucop=0.99

#fit=BicopSUR(Y1,X1,Y2,X2,Copula,type,rho,Lcop,Ucop)

