rm(list=ls(all=T))

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

##############

nelsonplosser <- read.csv("nelsonplosser.csv")

final_data <- nelsonplosser %>% filter(!is.na(bnd)) %>% select(sp,ip,bnd)
  
final_data <- final_data %>% mutate(logsp = log(sp),logip = log(ip)) %>%
  mutate(laglogsp = lag(logsp),laglogip = lag(logip), lagbnd = lag(bnd)) %>% 
  mutate(target = logsp -laglogsp, x1 = logip - laglogip, x2 = bnd - lagbnd)

final_data <- final_data[-1,]

lm_model <- lm(target~x1+x2,final_data)

summary(lm_model)

acf(lm_model$residuals)

acf(lm_model$residuals^2)

ts_residual <- auto.arima(lm_model$residuals)  ## ARIMA(1,0,1)

garch_model <- garchFit(formula = ~arma(0,1) + garch(1,1),data = lm_model$residuals)

garch_model@fit$coef

acf(garch_model@residuals/garch_model@sigma.t)

std_res <- garch_model@residuals/garch_model@sigma.t

qqnorm(std_res,datax = T)
qqline(std_res,datax = T)

new_lm_model <- lm(target~x1+x2,final_data,weights = 1/garch_model@sigma.t^2)  ## reciiprocal of conditional variance as wts

summary(new_lm_model)

acf(new_lm_model$residuals)
