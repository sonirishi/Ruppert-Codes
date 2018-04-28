rm(list=ls(all=T))
setwd("E:/Documents/Practice")

library(ggplot2)
library(dplyr)

gamma = rgamma(10000,scale = 2,shape = 3)

density_gamma = density(gamma)

plot(density_gamma$y)

scaledgamma = 1/2*(density(gamma/2)$y)  ## same as the original

print(sd(gamma))

scaledgamma_1 = (density(gamma/2)$y)

logdata = rlnorm(1000, meanlog = 0, sdlog = 1)

plot(density(logdata)$y)

plot(density(log(logdata))$y)

norm = rnorm(999,0,1) + jitter(999)

plot(density(norm)$y)

library(tseries)

jarque.bera.test(norm)  ## Tests for kurtosis and Skewness.

library(rmutil)

laplace = as.data.frame(rlaplace(1000,0,1)); colnames(laplace) = "laplace"

library(ggplot2)

norm = as.data.frame(rnorm(1000)); colnames(norm) = "norm"

library(EnvStats)

pareto = as.data.frame(rpareto(10000,location=1)); colnames(pareto) = "pareto"

ggplot() + geom_density(data=laplace,aes(x=laplace, colour = "laplace")) + 
  geom_density(data=norm,aes(x=norm, colour = "norm")) + 
  geom_density(data=pareto,aes(x=pareto, colour = "pareto"))

## Laplace has thicker tails

library(WRS)

norm = rnorm(1000); unif = runif(1000)

mixture = function(uni,norm){
  a = ifelse(uni < 0.9,norm,5*norm)
}

ndata = mapply(mixture, unif, norm)

plot(density(ndata)$y,col="red")

qqnorm(ndata,datax=T)

shapiro.test(ndata)

library(moments)

kurtosis(ndata)

library(sn)

tdata = rt(1000,10)

eps = 2

density_tdata = density(tdata)

density_eps = density(tdata*eps)
density_by_eps = density(tdata/eps)

a = approx(density_eps$x, density_eps$y, density_tdata$x, method = "linear")

b = approx(density_by_eps$x, density_by_eps$y, density_tdata$x, method = "linear")

skewed = function(y,density_eps,density_by_eps){
  skew_t = ifelse(y < 0, density_eps, density_by_eps)
}

skew_dense = mapply(skewed,density_tdata$x,a$y,b$y)

plot(skew_dense)

library(GoFKernel)

#inv_norm = inverse(norm)

(qnorm(0.7) - qnorm(0.3))/(pnorm(0.7) - pnorm(0.3))

###############################################################

ford = read.csv("ford.csv")

rets = diff(ford$FORD) 

library(fGarch)
library(MASS)
library(Ecdat)

data("Capm")

diffrf = diff(Capm$rf)

fit1 = stdFit(diffrf)  ## default fit is t distribution

fit2 = fitdistr(diffrf,"t") # supply the distribution you want to fit

chis = rchisq(1000,5)

gamma = rgamma(1000,shape=2.5,scale=2)

fit3 = sstdFit(diffrf)

fit4 = gedFit(diffrf)

normalfit = function(series){

  loglik = function(par) {
    f = -sum(log(dnorm(series, par[1], par[2])))
    f
  }
  results = optim(par=c(mean(series),sqrt(var(series))), loglik, method = "Nelder-Mead")
  return(results)
}

tfit = function(series){
  
  start = c(mean(series), sd(series), 2.1)
  
  loglik = function(par) {
    f = -sum(log(dstd(series, par[1], par[2], par[3])))
    f
  }
  
  lower = c(-1000,0.001,2.1)
  upper = c(1000000,1000000,60)
  results = optim(start, loglik, method = "L-BFGS-B",lower=lower,upper=upper)
  return(results)
}

tfit(diffrf)

# results = optim(par=c(mean(diffrf),sqrt(var(diffrf))), loglik, lower=c(-100,0.02), upper = c(100,2000))

tmix = function(series){
  
  start = c(mean(series), sd(series), 0.01, 2.1, 0.01)
  
  loglik = function(par) {
    f = -sum(log(par[5]*dstd(diffrf, par[1], par[2], par[4]) + 
                   (1-par[5])*dstd(diffrf, par[1], par[2] + par[3], par[4])))
    f
  }
  lower = c(-1000,0.001,0.001,2.1,0.01)
  upper = c(1000,1000,1000,60,0.99)
  results = optim(start, loglik, method = "L-BFGS-B", lower=lower, upper = upper)
  
  return(results)
}

tmix_res = tmix(diffrf)

nmix = function(series){
  
  start = c(mean(series), sd(series), 0.001, 0.01)
  
  loglik = function(par) {
    f = -sum(log(par[4]*dnorm(diffrf, par[1], par[2]) + 
                   (1-par[4])*dnorm(diffrf, par[1], (par[2] + par[3]))))
    f
  }
  lower = c(0,0.01,0.01,0.01)
  upper = c(1,1,1,0.99)
  results = optim(start, loglik, method = "L-BFGS-B", lower=lower, upper = upper)
  
  return(results)
}

nmis_res = nmix(diffrf)

setwd("E:/Documents/Practice")

gas = read.csv("GasFlowData.csv")

gasfit = sstdFit(gas$Flow1)

library(sn)

library(fGarch)

data = rbinom(10000,1,0.02)

dbetabinom <- function(k, n, a, b) {
  n2 <- ifelse(n > 100, 100, n)
  k2 <- round(k * n2 / n)
  beta(k2 + a, n2 - k2 + b) / beta(a, b)
}

betabinom_ll <- function(k, n, par) {
  sum(-log(dbetabinom(k, n, par[1], par[2])))
}

beta_mle <- function(...){
  par <- optim(par = c(1,1), fn=betabinom_ll, method="L-BFGS-B", lower=c(0.5,0.5), upper=c(500,500), ...)$par
  return(data.frame(a = par[1], b = par[2]))
}

#beta_mle(227,10000)

library(rmutil)

betabinom = function(wins, total){
  
  #start = c(mean(series), sd(series), 0.001, 0.01)
  
  loglik = function(par) {
    f = -sum(log(dbetabinom(wins,total,par[1],par[2])))
    f
  }
  lower = c(1,1)
  upper = c(200,200)
  results = optim(c(1,1), loglik, method = "L-BFGS-B", lower=lower, upper = upper)
  
  return(results)
}

betabinom(c(100,98,99,102,105),rep(200,5))  # this basically works on trials of multiple batches

betabinom(c(120,39,150,240,10),c(200,100,400,400,50))

### Bad attempt at box cox :(

box_n = function(series){
  
  start = c(1)
  
  loglik = function(par) {
    
    new_series = (series^par-1)/par
    
    mean = mean(new_series); std = sd(new_series)
    
    temp = density(new_series)
    
    density_norm = dnorm(temp$x,mean(temp$x),sd(temp$x))
    
    return(abs(sum(temp$y-density_norm)))
  }
  
  lower = c(1)
  upper = c(100)
  results = optim(start, loglik, method = "L-BFGS-B",lower=lower, upper=upper)
  
  return(results)
}

gasfit = box_n(gas$Flow1)

prof_log_lik=function(a){
     b=(optim(1,function(z) -sum(log(dgamma(x,a,z)))))$par
     return(-sum(log(dgamma(x,a,b))))
  }

library(caret)

BoxCoxTrans(gas$Flow1)

########################################
rm(list=ls(all=T))
setwd("E:/Documents/Practice")

library(ggplot2)
library(dplyr)
library(Ecdat)
library(data.table)

data("CPSch3")

dimnames(CPSch3)[[2]]

male.earnings = CPSch3 %>% filter(sex =="male") %>% select(ahe)

plot(density(male.earnings$ahe))

plot(density(log(male.earnings$ahe)))

plot(density(sqrt(male.earnings$ahe)))

qqnorm(male.earnings$ahe,datax = T)
qqline(male.earnings$ahe,datax = T)

qqnorm(log(male.earnings$ahe),datax = T)
qqline(log(male.earnings$ahe),datax = T)

qqnorm(sqrt(male.earnings$ahe),datax = T)
qqline(sqrt(male.earnings$ahe),datax = T)  # Best normal mimic with heavy tail

boxplot(male.earnings$ahe)
boxplot(log(male.earnings$ahe))
boxplot(sqrt(male.earnings$ahe))

plot(male.earnings$ahe,lwd=0.01)
abline(lsfit(1:nrow(male.earnings), male.earnings$ahe),col="red")  #least square fit with index

library(caret)

fitbx = BoxCoxTrans(male.earnings$ahe)

box.male.earn = predict(fitbx,male.earnings$ahe)

qqnorm(sqrt(box.male.earn),datax = T)
qqline(sqrt(box.male.earn),datax = T)  # not so good

library(tseries)

jarque.bera.test(box.male.earn)

jarque.bera.test(sqrt(male.earnings$ahe))

plot(density(box.male.earn))

library(MASS)

boxcox(male.earnings$ahe~1,lambda=seq(.3,.45,1/100))

bc = boxcox(male.earnings$ahe~1,lambda=seq(.3,.45,1/100),interp=F)

ind = (bc$y==max(bc$y))

ind2 = (bc$y > (max(bc$y) - qchisq(.95,df=1)/2))  # boxcox regression chi square distribution

bc2 = boxcox(male.earnings$ahe~1,lambda=seq(.3,.45,1/100),interp=T)

library(fGarch)

fit = sstdFit(male.earnings$ahe)

tmix = function(series){
  
  start = c(mean(series), sd(series), 0.01, 2.1, 0.01)
  
  loglik = function(par) {
    f = -sum(log(par[5]*dstd(series, par[1], par[2], par[4]) + 
                   (1-par[5])*dstd(series, par[1], par[2] + par[3], par[4])))
    f
  }
  lower = c(-1000,0.001,0.001,2.1,0.01)
  upper = c(1000,1000,1000,60,0.99)
  results = optim(start, loglik, method = "L-BFGS-B", lower=lower, upper = upper)
  
  return(results)
}

tmix_res = tmix(male.earnings$ahe)

plot(density(male.earnings$ahe))
par(new=T)
plot(male.earnings$ahe,
     dsstd(male.earnings$ahe,fit$estimate[1],fit$estimate[2],fit$estimate[3],
           fit$estimate[4]),col="red")
par(new=T)

skged <- function(series){
  start <- c(10,7,2.01)
  
  optim_func <- function(par){
    loglik <- -sum(log(dged(series,par[1],par[2],par[3])))
    return(loglik)
  }
  lower1 <- c(-10,0,1)
  upper1 <- c(1000,1000,100)
  
  vals <- optim(start, fn = optim_func,lower=lower1,upper = upper1,method="L-BFGS-B")
  return(vals)
}

gedfit <- skged(male.earnings$ahe)

r_gedfit <- fGarch::gedFit(male.earnings$ahe)

skged.s <- function(series){
  start <- c(10,7,2.01,0.1)
  
  optim_func <- function(par){
    loglik <- -sum(log(dsged(series,par[1],par[2],par[3],par[4])))
    return(loglik)
  }
  lower1 <- c(-10,0,1,1)
  upper1 <- c(1000,1000,100,100)
  
  vals <- optim(start, fn = optim_func,lower=lower1,upper = upper1,method="L-BFGS-B")
  return(vals)
}

gedfit.s <- skged.s(male.earnings$ahe)

r_gedfit.s <- fGarch::sgedFit(male.earnings$ahe)

plot(density(male.earnings$ahe))
par(new=T)
plot(male.earnings$ahe,
     dsged(male.earnings$ahe,gedfit.s$par[1],gedfit.s$par[2],gedfit.s$par[3],
           gedfit.s$par[4]),col="red")
par(new=T)
plot(male.earnings$ahe,
     dsstd(male.earnings$ahe,fit$estimate[1],fit$estimate[2],fit$estimate[3],
           fit$estimate[4]),col="blue")

data("Garch")

data("EuStockMarkets")

y <- diff(log(EuStockMarkets[,1]))

library(zoo)

rollmean <- rollapply(EuStockMarkets[,1],25,mean)
rollsd <- rollapply(EuStockMarkets[,1],25,sd)

plot(rollmean,col="red")
par(new=T)
plot(rollsd,col="blue")

tfit <- function(series){
  
  start <- c(mean(series), sd(series), 2.1)
  
  loglik <- function(par) {
    f <- -sum(log(dstd(series, par[1], par[2], par[3])))
    f
  }
  
  lower <- c(-1000,0.001,2.1)
  upper <- c(1000000,1000000,60)
  results <- optim(start, loglik, method = "L-BFGS-B",lower=lower,upper=upper,hessian=T)
  return(results)
}

tfit_euro <- tfit(y)

### TKDE -- Bad attempt

t.euro.param <- tfit_euro$par

y.n = qnorm(pstd(y,t.euro.param[1],t.euro.param[2],t.euro.param[3]))

y.n.dense = density(y.n)

denom = dnorm(qnorm(pstd(y,t.euro.param[1],t.euro.param[2],t.euro.param[3])))

dense.y = dstd(y,t.euro.param[1],t.euro.param[2],t.euro.param[3])

deriv = dense.y/denom

orig.dense = dnorm(y.n)


library(Ecdat)

data("CRSPday")

r <- CRSPday[,5]

plot(r)
mode(r); class(r)

cov(CRSPday[,4:6])
cor(CRSPday[,4:6])
apply(CRSPday[,4:6],2,mean)

gas <- read.csv("GasFlowData.csv")

ssdFit(gas$Flow1)

library(evir)

data("bmw")

bmw.tfit <- tfit(bmw)

fisher = solve(bmw.tfit$hessian)  # this calculates the inverse which is fisher information

eigen(bmw.tfit$hessian)  # Positive definite matrix

data("siemens")

plot(siemens)

siemens.fit <- tfit(siemens)
