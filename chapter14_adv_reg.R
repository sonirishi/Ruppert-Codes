source("E:/Documents/Practice/standard_code.R")

library(Ecdat)

data("Icecream")

model <- lm(cons~income+price+temp,data=Icecream)

summary(model)

library(car)

durbinWatsonTest(model)

acf(resid(model))

acf(rstudent(model))

model_1 <- arima(x=Icecream$cons,order = c(1,0,0),
                 xreg = cbind(Icecream$income,Icecream$price,Icecream$temp))

print(model_1)

model_2 <- arima(x=Icecream$cons,order = c(0,0,1),
                 xreg = cbind(Icecream$income,Icecream$price,Icecream$temp))

print(model_2)

########

data_1 <- Icecream %>% mutate(conslag = lag(cons))

model_new <- lm(cons~conslag+income+price+temp,data=data_1)

model_1 <- arima(x=Icecream$cons,order = c(1,0,0))

model_2 <- lm(cons~conslag,data=data_1)  

### xreg basically is linear regression but with ARMA errors. 
# Lagged dependent variable seems equivalent to ARIMAX

phi <- 0.5

X <- seq(-10,10,1)

sigma_mat <- toeplitz(phi^(0:20))

XprimeXinv <- solve(t(X)%*%X)

CovB <- XprimeXinv%*%t(X)%*%sigma_mat%*%X%*%XprimeXinv

ratio_SE <- sqrt(CovB)/sqrt(XprimeXinv)
