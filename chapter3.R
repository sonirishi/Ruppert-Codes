bond_val = function(coup,tm,ym,par_val){
  bv = coup/ym + (par_val - coup/ym)*(1-ym)^(-2*tm)
  return(bv)
}

price = 1200
C = 40
tim = 30
par = 1000

r = seq(0.02,0.05,length=300)

value = bond_val(C,tim,r,par)

y2m = spline(value,r,xout=price)

plot(r,value,type="l")
abline(b=1200)
abline(v=y2m)

y2m = spline(value,r)

##### Fitting linear on it for test and checking the output of y2m ############

lm_model = lm(formula = r~value)

summary(lm_model)

print(4.531e-06*1200 + 3.245e-02)

###################################

uniroot(function(r) r^2-0.5,c(0.7,0.8))

uniroot(function(x) -1200+bond_val(C,tim,x,par),c(0.02,0.05))

uniroot(function(x) -9800+bond_val(280,8,x,10000),c(0.01,0.5))

## Ques 5

library(fEcofin)

if( save_plots ){ postscript("../../WriteUp/Graphics/Chapter3/chap_3_prob_5_short_end.eps", onefile=FALSE, horizontal=FALSE) }
plot(mk.maturity[,1], mk.zero2[5,2:56], type="l", xlab="maturity", ylab="yield", xlim=c(0,3) )
lines(mk.maturity[,1], mk.zero2[6,2:56], lty=2, type="l")
lines(mk.maturity[,1], mk.zero2[7,2:56], lty=3, type="l")
lines(mk.maturity[,1], mk.zero2[8,2:56], lty=4, type="l")

#####################################

integrate(function(x){0.028 + 0.00042*x},lower=0,upper=20)


fw_rate = function(t){
  return(0.022 + 0.005*t - 0.004*t^2 + 0.0003*t^3)
}

DT = sapply(seq(1/2,4,1/2),function(x){exp(-1*integrate(fw_rate,lower=0,upper=x)$value)})

price = c(rep(21,7),1021)*DT

print(sum(price))

sum(price*seq(1/2,4,1/2))/sum(price)

## Ques 17

spot_rate = c(0.025,0.029,0.031,0.035)
tim = seq(1/2,2,1/2)

par = 1000

coupon = 35

par_vector = c(rep(coupon,3),par+coupon)

bond = list()

for (i in 1:4){
  bond[i] = par_vector[i]/(1+spot_rate[i])^tim[i]
}

print(sum(unlist(bond)))

### Ques 3

yield = integrate(function(x){0.032 + 0.001*x + 0.0002*x^2},lower=0,upper=5)$value/5

print(paste0("Price of Bond ", 1000*exp(-5*yield)))

## 18

sp_1 = 1000/980.39 - 1
sp_2 = (1000/957.41)^(1/2) - 1
sp_3 = (1000/923.18)^(1/3) - 1
sp_4 = (1000/888.489)^(1/4) - 1

spot = c(sp_1,sp_2,sp_3,sp_4)

coupon = 21

price =0
for (i in 1:4){
  price = price + ifelse(i<=3,coupon/(1+spot[i])^i,(coupon + 1000)/(1+spot[i])^4)
}

print(price)  # 967

#19
uniroot(function(x){1020 - (21/x + (1000 - 21/x)*(1+x)^(-2*2))},lower = 0.001, upper = 1)

## 7

par = 1000
price = 818

r = log(1000/818)/5

price_new = 1000/exp(4*0.042)

(price_new - 828)/828

## 6

par = 1000
time = 20
coupon = 1000*8.5/200

rate_new = 7.6/100

coupon/rate_new + (par - coupon/rate_new)*(1+rate_new)^(-2*(time-1))

coupon/rate_new + (par - coupon/rate_new)*(1+rate_new)^(-2*(time-1/2))

## Ques 11

r_t = function(t){
  return(ifelse(t-10 > 0,0.03 + 0.001*t - 0.00021*(t-10),0.03 + 0.001*t))
}

y_t = integrate(r_t,lower=0,upper=20)$value/20

d_t = exp(-20*y_t)

price = 100*d_t

print(price)

## Ques 10

r_t = function(t){
  0.035 + 0.0013*t
}

y_t = integrate(r_t,lower=0,upper=15)$value/15

d_t = exp(-15*y_t)

price = 100*d_t

print(price)

## Ques 8

par = 1000
time = 10
coupon = 22

interest_rate = 4/100

22/interest_rate + (1000 - 22/interest_rate)*(1+interest_rate)^(-2*10)

# coupon rate = 4.4 while y2m is 8% hence below par

### Ques 9

par = 1000
time = 7
price = 1050
coupon = 24

uniroot(function(x){price - (coupon/x + (par - coupon/x)*(1+x)^(-2*time))},lower = 0.01, upper = 1)


