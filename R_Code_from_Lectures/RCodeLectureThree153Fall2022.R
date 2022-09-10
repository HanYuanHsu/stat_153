#In the last lecture, we fit a linear trend model to the US population data via linear regression:
uspop.raw = read.csv("POPTHM.csv")
#Data downloaded from FRED. Monthly Data. Data given for each month equals the average of the estimated population on the first day of the month and the first day of the next month. The units are thousands of dollars so 200,000 actually refers to 200 million. 
plot(uspop.raw$POPTHM, type = "l")
uspop.ts = ts(uspop.raw$POPTHM, start = c(1959, 1), end = c(2022, 6), frequency = 12)
plot(uspop.ts, ylab = "Population (in thousands)", xlab = "Time (in months)", main = "US Population")

#The Linear Trend model can be fit to this data in the following way:
t = 1: length(uspop.ts)
lin.model = lm(uspop.ts ~ 1 + t)
plot(t, uspop.ts, type = "l", xlab = "Time (months)", ylab = "US Population", main = "Population of the United States")
points(t, lin.model$fitted, type = "l", col = "red")
summary(lin.model)

n = length(uspop.ts)
#The least squares estimates are: 
bhat = lin.model$coefficients

#Estimate of sigma:
sighat = sqrt((sum((lin.model$residuals)^2))/(n-2)) #this is also denoted by the Residual Standard Error

summary(lin.model)$coefficients #this gives the standard errors also
stderrs = c(193.210975, 0.438743)
#These standard errors are calculated as:
X = model.matrix(lin.model)
X
Sigma.mat = (sighat^2)*solve(t(X) %*% X)
Sigma.mat #the off-diagonal term is negative. Does this make sense? When the slope increases, the intercept decreases which makes sense in this dataset. 
sqrt(diag(Sigma.mat)) #these are the standard errors.

#Posterior densities:
#Let us first plot the univariate t-density
vals = seq(-5, 5, 0.01)
plot(vals, dt(vals, df = n-2), type = "l", xlab = "x", ylab = "density values", main = "Standard t-density")
#posterior density for intercept parameter:
j = 1
vals = c(bhat[j] - stderrs[j]*seq(5, 0, length = 100), bhat[j], bhat[j] + stderrs[j]*seq(0, 5, length = 100))
plot(vals, (1/stderrs[j])*(dt((vals - bhat[j])/(stderrs[j]), df = n-2)), type = "l", xlab = "x", ylab = "density values", main = "Standard t-density")
abline(v = bhat[j])
points(vals, dnorm(vals, mean = bhat[j], sd = stderrs[j]), col = "red", type = "l")
#posterior density for slope parameter:
j = 2
vals = c(bhat[j] - stderrs[j]*seq(5, 0, length = 100), bhat[j], bhat[j] + stderrs[j]*seq(0, 5, length = 100))
plot(vals, (1/stderrs[j])*(dt((vals - bhat[j])/(stderrs[j]), df = n-2)), type = "l", xlab = "x", ylab = "density values", main = "Standard t-density")
abline(v = bhat[j])
points(vals, dnorm(vals, mean = bhat[j], sd = stderrs[j]), col = "red", type = "l")

#Another way to quantify uncertainty
plot(t, uspop.ts, type = "n", xlab = "Time (months)", ylab = "US Population", main = "Population of the United States")
points(t, lin.model$fitted, type = "l", col = "red")
#Draw b0 and b1 from posterior distribution: 
N = 15
library(mvtnorm)
post.samples = (rbind(bhat)[rep(1, N), ]) + (rmvnorm(N, sigma = Sigma.mat))/(sqrt((rchisq(N, df = n-2))/(n-2)))
for(k in 1:N)
{
   abline(a = post.samples[k, 1], b = post.samples[k, 2], col = "blue")
}

#Another dataset: Google Trends (monthly data) for the query "Amazon"
amazon.raw = read.delim("amazontrends25Aug2022.csv")
#the useful data is between row 2 and row 225
amazon.use = amazon.raw[2:225,1]
len = length(amazon.use)
amazon = rep(0, len)
for(i in 1:len)
{
  amazon[i] = as.numeric(unlist(strsplit(as.character(amazon.use[i]), ","))[2])
}
amazon.ts = ts(amazon, start = c(2004, 1), end = c(2022, 8), frequency = 12)
plot(amazon, type = "l", ylab = "Amazon Search Popularity", xlab = "Month", main = "Google Trends Data for Amazon") #this dataset shows both trend and seasonality.

t = 1: length(amazon.ts)
lin.model = lm(amazon.ts ~ 1 + t)
plot(c(t, t), c(amazon.ts, lin.model$fitted), type = "n", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, amazon.ts, type = "l", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, lin.model$fitted, type = "l", col = "red")
summary(lin.model)

#Estimates
bhat = lin.model$coefficients
n = length(amazon.ts)
#Estimate of sigma:
sighat = sqrt((sum((lin.model$residuals)^2))/(n-2)) #this is also denoted by the Residual Standard Error

summary(lin.model)$coefficients #this gives the standard errors also

#These standard errors are calculated as:
X = model.matrix(lin.model)
X
Sigma.mat = (sighat^2)*solve(t(X) %*% X)
Sigma.mat #the off-diagonal term is negative. Does this make sense? When the slope increases, the intercept decreases which makes sense in this dataset. 
stderrs = sqrt(diag(Sigma.mat)) #these are the standard errors.

plot(c(t, t), c(amazon.ts, lin.model$fitted), type = "n", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, amazon.ts, type = "l", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, lin.model$fitted, type = "l", col = "red")
#Draw b0 and b1 from posterior distribution: 
N = 30
library(mvtnorm)
post.samples = (rbind(bhat)[rep(1, N), ]) + (rmvnorm(N, sigma = Sigma.mat))/(sqrt((rchisq(N, df = n-2))/(n-2)))
for(k in 1:N)
{
   abline(a = post.samples[k, 1], b = post.samples[k, 2], col = "blue")
}


#A Simulated dataset: 
n = 400
t = 1:n
sig = 1000
dt = 5 + 0.8*((t-(n/2))^2) + sig*rnorm(n)
plot(t, dt, type = "l", xlab = "Time", ylab = "Data", main = "A Simulated Time Series")

#Let us fit a linear (NOT quadratic) trend model to this dataset. 
t = 1: length(dt)
lin.model = lm(dt ~ 1 + t)
plot(c(t, t), c(dt, lin.model$fitted), type = "n", xlab = "Time", ylab = "Data", main = "A Simulated Dataset")
points(t, dt, type = "l")
points(t, lin.model$fitted, type = "l", col = "red")
summary(lin.model)

#Estimates
bhat = lin.model$coefficients
n = length(dt)
#Estimate of sigma:
sighat = sqrt((sum((lin.model$residuals)^2))/(n-2)) #this is also denoted by the Residual Standard Error

summary(lin.model)$coefficients #this gives the standard errors also

#These standard errors are calculated as:
X = model.matrix(lin.model)
X
Sigma.mat = (sighat^2)*solve(t(X) %*% X)
Sigma.mat #the off-diagonal term is negative. Does this make sense? When the slope increases, the intercept decreases which makes sense in this dataset. 
stderrs = sqrt(diag(Sigma.mat)) #these are the standard errors.

plot(c(t, t), c(dt, lin.model$fitted), type = "n", xlab = "Time", ylab = "Data", main = "Simulated Data") 
points(t, dt, type = "l")
points(t, lin.model$fitted, type = "l", col = "red")
#Draw b0 and b1 from posterior distribution: 
N = 30
library(mvtnorm)
post.samples = (rbind(bhat)[rep(1, N), ]) + (rmvnorm(N, sigma = Sigma.mat))/(sqrt((rchisq(N, df = n-2))/(n-2)))
for(k in 1:N)
{
   abline(a = post.samples[k, 1], b = post.samples[k, 2], col = "blue")
}

#More trend models:  
#Quadratic model:
n = 400
t = 1:n
sig = 1000
dt = 5 + 0.8*((t-(n/2))^2) + sig*rnorm(n)
plot(t, dt, type = "l", xlab = "Time", ylab = "Data", main = "A Simulated Time Series")

#Let us fit a linear (NOT quadratic) trend model to this dataset. 
t = 1: length(dt)
quad.model = lm(dt ~ 1 + t + I(t^2))
plot(c(t, t), c(dt, quad.model$fitted), type = "n", xlab = "Time", ylab = "Data", main = "A Simulated Dataset")
points(t, dt, type = "l")
points(t, quad.model$fitted, type = "l", col = "red")
summary(quad.model)

#Estimates
bhat = quad.model$coefficients
bhat
n = length(dt)
#Estimate of sigma:
sighat = sqrt((sum((quad.model$residuals)^2))/(n-3)) #this is also denoted by the Residual Standard Error (note that the divisor is n-3 and not n-2)
sighat
summary(quad.model)$coefficients #this gives the standard errors also

#These standard errors are calculated as:
X = model.matrix(quad.model)
X
Sigma.mat = (sighat^2)*solve(t(X) %*% X)
Sigma.mat #the off-diagonal term is negative. Does this make sense? When the slope increases, the intercept decreases which makes sense in this dataset. 
stderrs = sqrt(diag(Sigma.mat)) #these are the standard errors.
stderrs
summary(quad.model)$coefficients

plot(c(t, t), c(dt, quad.model$fitted), type = "n", xlab = "Time", ylab = "Data", main = "Simulated Data") 
points(t, dt, type = "l")
points(t, quad.model$fitted, type = "l", col = "red")
#Draw b0 and b1 from posterior distribution: 
N = 30
library(mvtnorm)
post.samples = (rbind(bhat)[rep(1, N), ]) + (rmvnorm(N, sigma = Sigma.mat))/(sqrt((rchisq(N, df = n-3))/(n-3)))
for(k in 1:N)
{
   points(t, post.samples[k,1] + t*post.samples[k,2]+(t^2)*post.samples[k,3], type = "l", col = "blue")
}

#Fitting a seasonal trend with linear regression
#A dataset with seasonality
USAccDeaths
plot(USAccDeaths, type = "o", ylab = "Deaths", main = "Monthly Totals of Accidental Deaths in the US 1973-1978")

#Parametric Functions for Seasonality
#Simplest Periodic Functions of period $d$ are cos(2 pi f t/d) and sin (2 pi f t/d)
d = 12 
f = 2 #Plot with values of f from 1 to 12.
fcos = function(t){cos(2*pi*f*t/d)}
fsin = function(t){sin(2*pi*f*t/d)}
max.t = 72
tme = 0:max.t
par(mfrow=c(2,1))
plot(fcos, 0, max.t)
points(tme, fcos(tme))
plot(fsin, 0, max.t)
points(tme, fsin(tme))

#Linear Combinations of these:

f1 = 1; a1 = 1
f2 = 2; a2 = 5
fval = function(t){a1*cos(2*pi*f1*t/d) + a2*sin(2*pi*f2*t/d)}
plot(fval, 0, max.t)
points(tme, fval(tme))

#Fitting a parametric seasonality function to the US Accidents Data Set.
plot(USAccDeaths, type = "o", ylab = "Deaths", main = "Monthly Totals of Accidental Deaths in the US 1973-1978")

t = 1: length(USAccDeaths)
#Start with the slowest frequency
f = 1
d = 12
v1 = cos(2*pi*f*t/d)
v2 = sin(2*pi*f*t/d)
lin.mod = lm(USAccDeaths ~ 1 + v1 + v2)
plot(t, USAccDeaths, type = "o", xlab = "Time", ylab = "Deaths", main = "Monthly Totals of Accidental Deaths in the US 1973-1978")
points(t, lin.mod$fitted, type = "l", col = "red")
par(mfrow = c(1, 1))

#Residuals:
plot(t, lin.mod$residuals, type = "o", xlab = "Year", ylab = "Residuals after fitting the seasonality function")
sd(lin.mod$residuals)
#There is some trend (quadratic?) still remaining in the residuals. 


#Fit two frequencies
t = 1: length(USAccDeaths)
f1 = 1
f2 = 2
d = 12
v1 = cos(2*pi*f1*t/d)
v2 = sin(2*pi*f1*t/d)
v3 = cos(2*pi*f2*t/d)
v4 = sin(2*pi*f2*t/d)
lin.mod = lm(USAccDeaths ~ 1 + v1 + v2 + v3 + v4)
plot(t, USAccDeaths, type = "o", xlab = "Time", ylab = "Deaths", main = "Monthly Totals of Accidental Deaths in the US 1973-1978")
points(t, lin.mod$fitted, type = "l", col = "red")

#Residuals:
plot(t, lin.mod$residuals, type = "o", xlab = "Year", ylab = "Residuals after fitting the seasonality function")
sd(lin.mod$residuals)


#Fit three frequencies: 
t = 1: length(USAccDeaths)
f1 = 1
f2 = 2
f3 = 3
d = 12
v1 = cos(2*pi*f1*t/d)
v2 = sin(2*pi*f1*t/d)
v3 = cos(2*pi*f2*t/d)
v4 = sin(2*pi*f2*t/d)
v5 = cos(2*pi*f3*t/d)
v6 = sin(2*pi*f3*t/d)
lin.mod = lm(USAccDeaths ~ 1 + v1 + v2 + v3 + v4 + v5 + v6)
plot(t, USAccDeaths, type = "o", xlab = "Time", ylab = "Deaths", main = "Monthly Totals of Accidental Deaths in the US 1973-1978")
points(t, lin.mod$fitted, type = "l", col = "red")

#Residuals:
plot(t, lin.mod$residuals, type = "o", xlab = "Year", ylab = "Residuals after fitting the seasonality function")
sd(lin.mod$residuals)

#Fitting a quadratic plus seasonal trend to this data: 
t = 1: length(USAccDeaths)
f1 = 1
f2 = 2
f3 = 3
d = 12
v1 = cos(2*pi*f1*t/d)
v2 = sin(2*pi*f1*t/d)
v3 = cos(2*pi*f2*t/d)
v4 = sin(2*pi*f2*t/d)
v5 = cos(2*pi*f3*t/d)
v6 = sin(2*pi*f3*t/d)
v7 = t
v8 = t^2
lin.mod = lm(USAccDeaths ~ 1 + v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8)
plot(t, USAccDeaths, type = "o", xlab = "Time", ylab = "Deaths", main = "Monthly Totals of Accidental Deaths in the US 1973-1978")
points(t, lin.mod$fitted, type = "l", col = "red")

plot(t, lin.mod$residuals, type = "o", xlab = "Year", ylab = "Residuals from the final Model")


#Linear/quadratic plus seasonal trend : Google Trends (monthly data) for the query "Amazon"
amazon.raw = read.delim("amazontrends25Aug2022.csv")
#the useful data is between row 2 and row 225
amazon.use = amazon.raw[2:225,1]
len = length(amazon.use)
amazon = rep(0, len)
for(i in 1:len)
{
  amazon[i] = as.numeric(unlist(strsplit(as.character(amazon.use[i]), ","))[2])
}
amazon.ts = ts(amazon, start = c(2004, 1), end = c(2022, 8), frequency = 12)
plot(amazon, type = "l", ylab = "Amazon Search Popularity", xlab = "Month", main = "Google Trends Data for Amazon") #this dataset shows both trend and seasonality.

t = 1: length(amazon.ts)
f1 = 1
d = 12
v1 = cos(2*pi*f1*t/d)
v2 = sin(2*pi*f1*t/d)
lin.model = lm(amazon.ts ~ 1 + t + v1 + v2)
plot(c(t, t), c(amazon.ts, lin.model$fitted), type = "n", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, amazon.ts, type = "l", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, lin.model$fitted, type = "l", col = "red")
summary(lin.model)

t = 1: length(amazon.ts)
f1 = 1
f2 = 2
d = 12
v1 = cos(2*pi*f1*t/d)
v2 = sin(2*pi*f1*t/d)
v3 = cos(2*pi*f2*t/d)
v4 = sin(2*pi*f2*t/d)
lin.model = lm(amazon.ts ~ 1 + t + v1 + v2 + v3 + v4)
plot(c(t, t), c(amazon.ts, lin.model$fitted), type = "n", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, amazon.ts, type = "l", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, lin.model$fitted, type = "l", col = "red")
summary(lin.model)

#quadratic plus seasonal
t = 1: length(amazon.ts)
f1 = 1
f2 = 2
d = 12
v1 = cos(2*pi*f1*t/d)
v2 = sin(2*pi*f1*t/d)
v3 = cos(2*pi*f2*t/d)
v4 = sin(2*pi*f2*t/d)
lin.model = lm(amazon.ts ~ 1 + t + I(t^2) + v1 + v2 + v3 + v4)
plot(c(t, t), c(amazon.ts, lin.model$fitted), type = "n", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, amazon.ts, type = "l", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, lin.model$fitted, type = "l", col = "red")
summary(lin.model)

#quadratic plus seasonal
t = 1: length(amazon.ts)
f1 = 1
f2 = 2
f3 = 3
d = 12
v1 = cos(2*pi*f1*t/d)
v2 = sin(2*pi*f1*t/d)
v3 = cos(2*pi*f2*t/d)
v4 = sin(2*pi*f2*t/d)
v5 = cos(2*pi*f3*t/d)
v6 = sin(2*pi*f3*t/d)
lin.model = lm(amazon.ts ~ 1 + t + I(t^2) + v1 + v2 + v3 + v4 + v5 + v6)
plot(c(t, t), c(amazon.ts, lin.model$fitted), type = "n", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, amazon.ts, type = "l", xlab = "Time (months)", ylab = "Trend Popularity", main = "Google trends for query Amazon")
points(t, lin.model$fitted, type = "l", col = "red")
summary(lin.model)

