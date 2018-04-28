rm(list=ls(all=TRUE))

setwd("E:/Documents/Practice")

dat = read.csv("Stock_FX_Bond.csv")

library(ggplot2)

ggplot(data = dat, aes(x = dat$Date, y = dat$GM_AC)) + geom_point() + xlab("Date") + ylab("GM_AC")

ggplot(data = dat, aes(x = dat$Date, y = dat$F_AC)) + geom_point() + xlab("Date") + ylab("F_AC")

GMReturn = dat$GM_AC[2:nrow(dat)]/dat$GM_AC[1:(nrow(dat)-1)]-1

FReturn = dat$F_AC[2:nrow(dat)]/dat$F_AC[1:(nrow(dat)-1)]-1

plot(GMReturn,FReturn)

cor(GMReturn,FReturn)

logGMReturns = log(1+GMReturn)

logFReturns = log(1+FReturn)

plot(logGMReturns,logFReturns)

cor(logGMReturns,logFReturns)

mean_return = 0.05/253
sd_return = 0.23/sqrt(253)

ind_below = list()

for (i in 1:10000){
  log_return = rnorm(45,mean_return,sd_return)
  log_price = log(1e6) + cumsum(log_return)
  min_logprice = min(log_price)
  ind_below[i] = ifelse(min_logprice < log(950000),1,0)
}

print(mean(unlist(ind_below)))

set.seed(2018)

iter = 1e6
days = 100
return = list()
above_hit = 0
below_hit = 0
middle_hit = 0
days_in_trade = list()
stra_return = list()

for(i in 1:iter){
  log_return = rnorm(days,mean_return,sd_return)
  log_price = log(1e6) + cumsum(log_return)  ## geometric random walk
  price = exp(log_price)
  for(j in 1:100){
    if(price[j] < 950000){
      below_hit = below_hit + 1
      return[i] = price[j] - 1000000
      days_in_trade[i] = j
      break
    } else if(price[j] >= 1100000){
      above_hit = above_hit + 1
      return[i] = price[j] - 1000000
      days_in_trade[i] = j
      break
    } else{
      next
    }
  }
  if(j==days){
    middle_hit = middle_hit + 1
    return[i] = price[days] - 1000000
    days_in_trade[i] = days
  }
  stra_return[i] = return[i][[1]]/(50000*days_in_trade[i][[1]])
}

return_new = unlist(return)

mean(above_hit)
mean(below_hit)

length(which(return_new < 0))/iter

mean(return_new)
mean(unlist(stra_return))

#####################################################

## Ques 1

mean = 0.001
std = 0.015
iter = 1e5
stock_to = 1000

score=0
for (i in 1:iter){
  lr = rnorm(1,mean,std)
  stock_t1 = stock_to*exp(lr)
  score = score + ifelse(stock_t1 < 990,1,0)
}

print(score/iter)

score=0
for (i in 1:iter){
  lr = rnorm(5,mean,std)
  stock_t1 = stock_to*exp(sum(lr))
  score = score + ifelse(stock_t1 < 990,1,0)
}

print(score/iter)

mean = 0.1
std = 0.2
stock_to = 100

score=0
for (i in 1:iter){
  lr = rnorm(1,mean,std)
  stock_t1 = stock_to*exp(lr)
  score = score + ifelse(stock_t1 > 110,1,0)
}

print(score/iter)

mean = 0.0002
std = 0.03
stock_to = 97

score=0
for (i in 1:iter){
  lr = rnorm(20,mean,std)
  stock_t1 = stock_to*exp(sum(lr))
  score = score + ifelse(stock_t1 > 100,1,0)
}

print(score/iter)



