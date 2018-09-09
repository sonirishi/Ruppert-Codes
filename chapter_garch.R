library(fGarch)
library(dplyr)

w <- 1; alpha1 <- 0.08; beta1 <- 0.9; phi <- 0.8
n <- 200; mu <- 0.1

a <- rep(0,n)
sigma <- rep(0,n)
u <- rep(0,n)

a[1] <- 0.5
sigma[1] <- 0.3
u[1] <- 0.5

whitenoise <- rnorm(n)

for(i in 2:n){
  sigma[i] <- sqrt(w + alpha1*a[i-1]^2 + beta1*sigma[i-1]^2)
  a[i] <- sigma[i]*whitenoise[i]
  u[i] <- mu + phi*(u[i-1]-mu) + a[i]
}

plot(u,type='l')
lines(whitenoise,type='l',col='green')

acf(a)

acf(a^2)

acf(u)

acf(diff(u))

##########

setwd("E:/Documents/Practice")

bmwretrun <- read.csv("bmwRet.csv")

garchmodel <- garchFit(formula = ~arma(1,0) + garch(1,1), data = bmwretrun$BMW.RET)

garchmodel@fit$coef

garchmodel@fit$matcoef   ### This is the standard error of the coefficients

std_residuals <- garchmodel@residuals/garchmodel@sigma.t

n <- nrow(bmwretrun)

library(MASS)

fit_t_residual <- fitdistr(std_residuals, densfun = "t")

qqnorm(std_residuals,datax = T)
qqline(std_residuals,datax = T)  ## This proves residuals are not normal

qqplot(std_residuals,qt((1:n)/(n+1),df=4))  ## Seems like a straight line so t is a better fit

summary(garchmodel)

#############

garchmodel_t <- garchFit(formula = ~arma(1,0) + garch(1,1), 
                       data = bmwretrun$BMW.RET, cond.dist = "std")

summary(garchmodel_t)

summary(garchmodel_t)  ## Best model as lower AIC the better

##############3

library(forecast)

setwd("E:/Documents/Practice")

bmwretrun <- read.csv("bmwRet.csv")

ar_model <- forecast::Arima(bmwretrun$BMW.RET,order=c(1,0,0))

rsidual <- ar_model$residuals

acf(rsidual^2)

###########

mean(rnorm(100000))

mean(rnorm(10000)^2)  ## expectation of chi square given normal

#### delta parameter be included as well for optimization; delta is the aparch parameter

aparchmodel_t <- garchFit(formula = ~arma(1,0) + aparch(1,1), 
                         data = bmwretrun$BMW.RET, cond.dist = "std",include.delta = T)

aparchmodel_t@fit$coef

aparchmodel_t@fit$se.coef  ### Std Error of coefficients

aparchmodel_t@fit$matcoef

