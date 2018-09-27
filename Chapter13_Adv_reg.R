source("E:/Documents/Practice/standard_code.R")

dat = read.table(file="WeekInt.txt",header=T)

dat <- dat %>% mutate(cm10_diff = cm10 - lag(cm10), aaa_diff = aaa - lag(aaa),
                      ff_diff = ff - lag(ff), cm30_diff = cm30 - lag(cm30))

model1 <- lm(aaa_diff~cm10_diff,data=dat)

model2 <- lm(aaa_diff~cm30_diff,data=dat)

residual_1 <- model1$residuals

residual_2 <- model2$residuals

acf(residual_1)

acf(residual_2)

model3 <- lm(aaa_diff~cm10_diff+cm30_diff,data=dat)

residual_3 <- model3$residuals

acf(residual_3)  ## Autocorrelation of level 1

qqnorm(residual_3,datax = T)
qqline(residual_3,datax = T)

library(car)

dw_result <- car::durbinWatsonTest(model3)

plot(y=residual_3,x=model3$fitted.values,col="green")

library(forecast)

arma_fit <- auto.arima(residual_3)

#############

model_price <- lm(aaa~cm10+cm30,data=dat)

acf(model_price$residuals)

library(MASS)

stud_residual_price <- studres(model_price)

plot(y=stud_residual_price,x=model_price$fitted.values)

library(tseries)

tseries::adf.test(stud_residual_price)

tseries::kpss.test(stud_residual_price)

Box.test(stud_residual_price)  ## Default test  independence hypo rejected

########### AR(1)

error1 <- rnorm(1000,0,1)

error2 <- rnorm(1000,0,1)

y <- rep(0,1000)
x <- rep(0,1000)

y[1] <- 0.5; x[1] <- 0.5

for(i in 2:1000){
  y[i] <- 0.99*y[i-1] + error1[i]
  x[i] <- 0.99*x[i-1] + error2[i]
}

model1 <- lm(y~x)

summary(model1)

# Coefficients:
# (Intercept)  2.65998    0.25607  10.388  < 2e-16 ***
#  x            0.14407    0.03577   4.028 6.06e-05 ***

error1 <- rnorm(1000,0,1)

error2 <- rnorm(1000,0,1)

y <- rep(0,1000)
x <- rep(0,1000)

y[1] <- 0.5; x[1] <- 0.5

for(i in 2:1000){
  y[i] <- 0.99*y[i-1] + error1[i]
  x[i] <- 0.99*x[i-1] + error2[i]
}

model1 <- lm(y~x)

summary(model1)

# Coefficients:
# (Intercept)  1.26782    0.19599   6.469 1.54e-10 ***
#  x           -0.18129    0.02481  -7.308 5.55e-13 ***

model2 <- lm(diff(y)~diff(x))

summary(model2)  ### not significant  This is again due to trending kind of concept in time series

########

library(AER)

data("CPS1988")

fitlm1 <- lm(wage~education+experience+ethnicity,data=CPS1988)

resid1 <- rstudent(fitlm1)

plot(fitlm1$fitted.values,resid1)

plot(density(resid1))

plot(fitlm1$fitted.values,abs(resid1))
lines(lowess(fitlm1$fitted.values,abs(resid1)),col="red")  ## not much bad variance

qqnorm(resid1,datax = T)
qqline(resid1,datax = T)   # not normal

fitlm1 <- lm(log(wage)~education+experience+ethnicity,data=CPS1988)

resid1 <- rstudent(fitlm1)

plot(fitlm1$fitted.values,resid1)

plot(density(resid1))

plot(fitlm1$fitted.values,abs(resid1))
lines(lowess(fitlm1$fitted.values,abs(resid1)),col="red")  ## not much bad variance

qqnorm(resid1,datax = T)
qqline(resid1,datax = T)   # not normal

########## Test residuals random forest

library(randomForest)

model2 <- randomForest(aaa_diff~cm30_diff+cm10_diff+ff_diff,data=dat[-1,],ntree=30)

lm_model <- lm(aaa_diff~cm30_diff+cm10_diff+ff_diff,data=dat[-1,])

X <- as.matrix(dat[-1,c("cm30_diff","cm10_diff","ff_diff")])

Y <- matrix(dat[-1,"aaa_diff"])

hat_actual <- lm.influence(lm_model)

hat_new <- hat_actual$hat

ypredict <- model2$predicted

rf_residual <- Y-ypredict

#########

fitlm1 <- lm(sqrt(wage)~education+experience+ethnicity,data=CPS1988)

resid1 <- rstudent(fitlm1)

plot(fitlm1$fitted.values,resid1)

plot(density(resid1))

plot(fitlm1$fitted.values,abs(resid1))
lines(lowess(fitlm1$fitted.values,abs(resid1)),col="red")  ## not much bad variance

qqnorm(resid1,datax = T)
qqline(resid1,datax = T)   # not normal


