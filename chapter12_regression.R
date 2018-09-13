rm(list=ls(all=T))

library(dplyr)

setwd("E:/Documents/Practice")

dat = read.table(file="WeekInt.txt",header=T)

dat <- dat %>% mutate(cm10_diff = cm10 - lag(cm10), aaa_diff = aaa - lag(aaa),
                      ff_diff = ff - lag(ff), cm30_diff = cm30 - lag(cm30))

lm_model1 <- lm(aaa_diff ~ cm10_diff, data = dat)

lm_model2 <- lm(aaa_diff ~ cm10_diff + cm30_diff + ff_diff, data = dat)

lm_model3 <- lm(aaa_diff ~ cm10_diff + cm30_diff, data = dat)

anova(lm_model1)

anova(lm_model2)

## Test that model2 and model 1 are same as in other coeff as 0

residual_diff <- sum(lm_model1$residuals^2) - sum(lm_model2$residual^2)

sigmasq2 <- sum(lm_model2$residuals^2)/(nrow(dat)-3-1)

Fstat <- residual_diff/sigmasq2  ## WHy is my fstat twice the anova

df(Fstat,2,876)

print(anova(lm_model1,lm_model2))

(summary(lm_model2)$sigma)**2 ## This is the error variance

residual_diff <- sum(lm_model3$residuals^2) - sum(lm_model2$residual^2)

Fstat <- residual_diff/sigmasq2

df(Fstat,1,876)

###################

library(car)

vif(lm_model2)

####

nelsonplosser <- read.csv("nelsonplosser.csv")

nelsonplosser <- nelsonplosser %>% select(sp,gnp.r,gnp.pc,ip,cpi,emp,bnd)

nelsonplosser1 <- nelsonplosser[complete.cases(nelsonplosser),]

nelsonplosser1 <- nelsonplosser1 %>% 
  mutate(target = log(sp) - lag(log(sp)), x1 = gnp.r - lag(gnp.r),
         x2 = gnp.pc - lag(gnp.pc), x3 = log(ip) - lag(log(ip)), 
         x4 = log(cpi) - lag(log(cpi)), x5 = emp - lag(emp), x6 = bnd - lag(bnd))

nelsonplosser2 <- nelsonplosser1[complete.cases(nelsonplosser1),]

nelson_model1 <- lm(target ~ x1+x2+x3+x4+x5+x6, data = nelsonplosser2)

vif(nelson_model1)

cor(nelsonplosser2$x1,nelsonplosser2$x2)

library(MASS)

stepAIC(nelson_model1)  ## It stops prematurely

nelson_model2 <- lm(target ~ x3+x6, data = nelsonplosser2)

summary(nelson_model2)

extractAIC(nelson_model2)  ## Best way is to loop through all models

