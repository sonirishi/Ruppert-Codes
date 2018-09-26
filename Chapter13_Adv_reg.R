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
