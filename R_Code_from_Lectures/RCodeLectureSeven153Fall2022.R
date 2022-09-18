#Change-point estimation:
n = 10000
mu1 = 0
mu2 = 0.2
dt = c(rep(mu1, n/2), rep(mu2, n/2)) + rnorm(n)
plot(dt, type = "l")

tme = 1:n
n = length(tme)
grid.res = 1
chnge.val = seq(3, n-3, grid.res)
X = matrix(1, nrow = n, ncol = 2)
expos = rep(-1, length(chnge.val)) #exact marginal posterior for omega
log.values = rep(-1, length(chnge.val))
for(i in 1:length(chnge.val))
{
    X[,2] = (tme > chnge.val[i])
    mod = lm(dt ~ X[,2])
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    log.values[i] = log.value
}
log.values = log.values - max(log.values) #scaling to remove large values 
expos = exp(log.values)
expos = (expos/sum(expos))/grid.res
plot(chnge.val, expos, type = "l")
ind.max = which.max(expos)
est.cp = chnge.val[ind.max] #posterior maximizer
est.cp


library(mvtnorm)
#Inference for the mean values before and after the change-point
N = 2000 #number of posterior samples
post.samples = matrix(-1, N, 4)
post.samples[,1] = sample(chnge.val, N, replace = T, prob = expos)
for(i in 1:N)
{
    cp = post.samples[i,1]
    X = matrix(1, nrow = n, ncol = 2)
    X[,2] = (tme > cp)
    lin.model = lm(dt ~ X[,2])
    bhat = lin.model$coefficients
    sighat = sqrt((sum((lin.model$residuals)^2))/(n-2)) #this is also denoted by the Residual Standard Error
    Sigma.mat = (sighat^2)*solve(t(X) %*% X)
    chiran = (rchisq(1, df = n-2))
    mu.samples = bhat + (rmvnorm(1, sigma = Sigma.mat))/(sqrt(chiran/(n-2)))
    mu.samples[2] = mu.samples[1] + mu.samples[2]
    sig.sample = sqrt((sum((lin.model$residuals)^2))/chiran)
    post.samples[i,2:3] = mu.samples
    post.samples[i, 4] = sig.sample
}
summary(post.samples[,2])
summary(post.samples[,3])
par(mfrow = c(3, 1))
hist(post.samples[,1], breaks = 100)
hist(post.samples[,2], breaks = 100)
hist(post.samples[,3], breaks = 100)
par(mfrow = c(1, 1))
summary(post.samples[,4]) #true value of sigma equals 1

#Plotting some posterior samples
plot(dt, type = "l")
abline(v = 5000, col = "blue")
for(i in 1:25)
{
    points(1:post.samples[i,1], rep(post.samples[i,2], post.samples[i,1]), col = "red", type = "l")
    points((1 + post.samples[i,1]):n, rep(post.samples[i,3], n-post.samples[i,1]), col = "red", type = "l")
}

#Multiple Change-Points
n = 400
th = c(rep(0, 100), rep(1.5, 100), rep(-1, 100), rep(0, 100))
th
plot(th, type = "l")

#Generate data: 
sig = 1
tme = 1:n
yy = th + sig*rnorm(n)
plot(tme, yy, type = "l")
#actual changepoints: 
cp.true = c(100, 200, 300)

#Let us first fit the single changepoint model to this data: 
tme = 1:n
n = length(tme)
grid.res = 1
chnge.val = seq(3, n-3, grid.res)
X = matrix(1, nrow = n, ncol = 2)
expos = rep(-1, length(chnge.val)) #exact marginal posterior for omega
log.values = rep(-1, length(chnge.val))
for(i in 1:length(chnge.val))
{
    X[,2] = (tme > chnge.val[i])
    mod = lm(yy ~ X[,2])
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    log.values[i] = log.value
}
log.values = log.values - max(log.values) #scaling to remove large values 
expos = exp(log.values)
expos = (expos/sum(expos))/grid.res
plot(chnge.val, expos, type = "l")
ind.max = which.max(expos)
est.cp = chnge.val[ind.max] #posterior maximizer
est.cp

#Inference for the mean values before and after the change-point
N = 2000 #number of posterior samples
post.samples = matrix(-1, N, 4)
post.samples[,1] = sample(chnge.val, N, replace = T, prob = expos)
for(i in 1:N)
{
    cp = post.samples[i,1]
    X = matrix(1, nrow = n, ncol = 2)
    X[,2] = (tme > cp)
    lin.model = lm(yy ~ X[,2])
    bhat = lin.model$coefficients
    sighat = sqrt((sum((lin.model$residuals)^2))/(n-2)) #this is also denoted by the Residual Standard Error
    Sigma.mat = (sighat^2)*solve(t(X) %*% X)
    chiran = (rchisq(1, df = n-2))
    mu.samples = bhat + (rmvnorm(1, sigma = Sigma.mat))/(sqrt(chiran/(n-2)))
    mu.samples[2] = mu.samples[1] + mu.samples[2]
    sig.sample = sqrt((sum((lin.model$residuals)^2))/chiran)
    post.samples[i,2:3] = mu.samples
    post.samples[i, 4] = sig.sample
}
summary(post.samples[,2])
summary(post.samples[,3])
par(mfrow = c(3, 1))
hist(post.samples[,1], breaks = 100)
hist(post.samples[,2], breaks = 100)
hist(post.samples[,3], breaks = 100)
par(mfrow = c(1, 1))
summary(post.samples[,4]) #true value of sigma equals 1

#Plotting some posterior samples
plot(yy, type = "l")
abline(v = 200, col = "blue")
for(i in 1:25)
{
    points(1:post.samples[i,1], rep(post.samples[i,2], post.samples[i,1]), col = "red", type = "l")
    points((1 + post.samples[i,1]):n, rep(post.samples[i,3], n-post.samples[i,1]), col = "red", type = "l")
}

#Posterior for Multiple Change Points:
logpost = function(cps)
{
    cps = sort(cps)
    k = length(cps)
    X = matrix(1, nrow = n, ncol = (k+1))
    for(j in 1:k)
    {
        X[,(j+1)] = as.numeric(tme > cps[j])
    }
    mod = lm(yy ~ X)
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    return(log.value)
}
logpost(cp.true)
logpost(200)

cp1.gr = seq(75, 125, by = 1)
cp2.gr = seq(175, 225, by = 1)
cp3.gr = seq(275, 325, by = 1)
g = expand.grid(x=cp1.gr, y=cp2.gr, z = cp3.gr)
head(g)
tail(g)
for(i in 1:nrow(g)) g$lp[i] = logpost(c(g$x[i], g$y[i], g$z[i]))
head(g)
tail(g)
#Convert LogPosterior to Posterior
g$p = exp(g$lp - max(g$lp))
g$p = g$p/(sum(g$p))
head(g)
tail(g)
#Posterior mode:
ind.max = which.max(g$p)
est.cp = g[ind.max,] #posterior maximizer
est.cp

#marginal distributions of each changepoint:
mar.cp1 = rep(-1, length(cp1.gr))
for(i in 1:length(cp1.gr))
{
  ind = (g$x == cp1.gr[i])
  mar.cp1[i] = sum(g$p[ind])
}
par(mfrow = c(3, 1))
plot(cp1.gr, mar.cp1, type = "h")

mar.cp2 = rep(-1, length(cp2.gr))
for(i in 1:length(cp2.gr))
{
  ind = (g$y == cp2.gr[i])
  mar.cp2[i] = sum(g$p[ind])
}
plot(cp2.gr, mar.cp2, type = "h")

mar.cp3 = rep(-1, length(cp3.gr))
for(i in 1:length(cp3.gr))
{
  ind = (g$z == cp3.gr[i])
  mar.cp3[i] = sum(g$p[ind])
}
plot(cp3.gr, mar.cp3, type = "h")
par(mfrow = c(1, 1))

#Generating posterior samples:
N = 2000 #number of posterior samples
post.samples = matrix(-1, N, 8) #there are 8 parameters in total (3 changepoints + 4 levels + sigma)
#First sample the change points:
smp.cp = g[sample(1:nrow(g), N, replace = T, prob = g$p),1:3]
rownames(smp.cp) = NULL
post.samples = cbind(smp.cp, matrix(-1, N, 5))
for(i in 1:N)
{
    cp = post.samples[i,1:3]
    k = length(cp)
    X = matrix(1, nrow = n, ncol = (k+1))
    for(j in 1:k)
    {
        X[,(j+1)] = as.numeric(tme > as.numeric(cp[j]))
    }
    lin.model = lm(yy ~ -1 + X)
    bhat = lin.model$coefficients
    sighat = sqrt((sum((lin.model$residuals)^2))/(n-4)) #this is also denoted by the Residual Standard Error
    Sigma.mat = (sighat^2)*solve(t(X) %*% X)
    chiran = (rchisq(1, df = n-4))
    levels = cumsum(bhat + (rmvnorm(1, sigma = Sigma.mat))/(sqrt(chiran/(n-4))))
    sig.sample = sqrt((sum((lin.model$residuals)^2))/chiran)
    post.samples[i,4:7] = levels
    post.samples[i, 8] = sig.sample
}

#Summary of the posterior values of sigma:
summary(post.samples[,8]) #true value of sigma is 1

#Posterior values of the 4 levels:
summary(post.samples[,4])
summary(post.samples[,5])
summary(post.samples[,6])
summary(post.samples[,7])
        
#Plotting some posterior fitted functions:
plot(tme, yy, type = "l")
#abline(v = cp.true[1], col = "blue")
#abline(v = cp.true[2], col = "blue")
#abline(v = cp.true[3], col = "blue")

for(i in 1:55)
{
    cp1 = post.samples[i,1]
    cp2 = post.samples[i,2]
    cp3 = post.samples[i,3]
    abline(v = cp1, col = "blue")
    abline(v = cp2, col = "blue")
    abline(v = cp3, col = "blue")
    lev1 = post.samples[i,4]
    lev2 = post.samples[i,5]
    lev3 = post.samples[i,6]
    lev4 = post.samples[i,7]
    points(1:cp1, rep(lev1, cp1), col = "red", type = "l")
    points((cp1+1):cp2, rep(lev2, cp2 - cp1), col = "red", type = "l")
    points((cp2+1):cp3, rep(lev3, cp3 - cp2), col = "red", type = "l")
    points((cp3+1):n, rep(lev4, n-cp3), col = "red", type = "l")
}

#Change of slope:
n = 1000
tme = 1:n
th = (c(0.5*tme[1:(n/2)], 0.9*tme[((n/2)+1):n] - (0.9 - 0.5)*(n/2)))/10
sig = 10
yy = th + sig*rnorm(n)
plot(tme, yy, type = "l")

#cp here will denote change of slope point
logpost = function(cps) 
{
    cps = sort(cps)
    k = length(cps)
    X = matrix(1, nrow = n, ncol = (k+2))
    X[,2] = tme
    for(j in 1:k)
    {
        X[,(j+2)] = pmax(tme - cps[j], 0)
    }
    mod = lm(yy ~ X)
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    return(log.value)
}
logpost(500)

#Posterior for the change of slope point:
csp.vals = 3:(n-3)
lp.vals = rep(100, length(csp.vals))
for(i in 1:length(csp.vals))
{
      lp.vals[i] = logpost(csp.vals[i])
}
expos = exp(lp.vals - max(lp.vals))
expos = expos/(sum(expos))
plot(csp.vals, expos, type = "h")

#Posterior samples for all the parameters:
N = 2000 #number of posterior samples
post.samples = matrix(-1, N, 5)
post.samples[,1] = sample(csp.vals, N, replace = T, prob = expos)
for(i in 1:N)
{
    cps = post.samples[i,1]
    k = length(cps)
    X = matrix(1, nrow = n, ncol = (k+2))
    X[,2] = tme
    for(j in 1:k)
    {
        X[,(j+2)] = pmax(tme - cps[j], 0)
    }
    lin.model = lm(yy ~ -1 + X)
    bhat = lin.model$coefficients
    sighat = sqrt((sum((lin.model$residuals)^2))/(n-3)) #this is also denoted by the Residual Standard Error
    Sigma.mat = (sighat^2)*solve(t(X) %*% X)
    chiran = (rchisq(1, df = n-3))
    beta.samples = bhat + (rmvnorm(1, sigma = Sigma.mat))/(sqrt(chiran/(n-3)))
    sig.sample = sqrt((sum((lin.model$residuals)^2))/chiran)
    post.samples[i,c(2, 3, 4)] = beta.samples
    post.samples[i, 5] = sig.sample
}
summary(post.samples[,5]) #true sigma is 10

#Plotting the fitted curves:
plot(tme, yy, type = "l")
for(i in 1:100)
{
    csp = post.samples[i,1]
    b0 = post.samples[i,2]
    b1 = post.samples[i,3]
    b2 = post.samples[i,4]
    points(tme, b0 + b1*tme + b2*(pmax(tme - csp, 0)), type = "l", col = "red")
    abline(v = csp, col = "blue")
}

#Real Dataset: US Population
uspop.raw = read.csv("POPTHM.csv")
#Data downloaded from FRED. Monthly Data. Data given for each month equals the average of the estimated population on the first day of the month and the first day of the next month. The units are thousands of dollars so 200,000 actually refers to 200 million. 
plot(uspop.raw$POPTHM, type = "l")
uspop.ts = ts(uspop.raw$POPTHM, start = c(1959, 1), end = c(2022, 6), frequency = 12)
plot(uspop.ts, ylab = "Population (in thousands)", xlab = "Time (in months)", main = "US Population")

yy = uspop.ts
n = length(yy)
tme = 1:n

#cp here will denote change of slope point
logpost = function(cps) 
{
    cps = sort(cps)
    k = length(cps)
    X = matrix(1, nrow = n, ncol = (k+2))
    X[,2] = tme
    for(j in 1:k)
    {
        X[,(j+2)] = pmax(tme - cps[j], 0)
    }
    mod = lm(yy ~ X)
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    return(log.value)
}
logpost(500) #example

#Posterior for the change of slope point:
csp.vals = 3:(n-3)
lp.vals = rep(100, length(csp.vals))
for(i in 1:length(csp.vals))
{
      lp.vals[i] = logpost(csp.vals[i])
}
expos = exp(lp.vals - max(lp.vals))
expos = expos/(sum(expos))
plot(csp.vals, expos, type = "h")
ind = which.max(expos)
csp.est = csp.vals[ind]
csp.est
abline(v = csp.est, col = "blue")

#Posterior samples for all the parameters:
N = 2000 #number of posterior samples
post.samples = matrix(-1, N, 5)
post.samples[,1] = sample(csp.vals, N, replace = T, prob = expos)
for(i in 1:N)
{
    cps = post.samples[i,1]
    k = length(cps)
    X = matrix(1, nrow = n, ncol = (k+2))
    X[,2] = tme
    for(j in 1:k)
    {
        X[,(j+2)] = pmax(tme - cps[j], 0)
    }
    lin.model = lm(yy ~ -1 + X)
    bhat = lin.model$coefficients
    sighat = sqrt((sum((lin.model$residuals)^2))/(n-3)) #this is also denoted by the Residual Standard Error
    Sigma.mat = (sighat^2)*solve(t(X) %*% X)
    chiran = (rchisq(1, df = n-3))
    beta.samples = bhat + (rmvnorm(1, sigma = Sigma.mat))/(sqrt(chiran/(n-3)))
    sig.sample = sqrt((sum((lin.model$residuals)^2))/chiran)
    post.samples[i,c(2, 3, 4)] = beta.samples
    post.samples[i, 5] = sig.sample
}
summary(post.samples[,5]) #true sigma is 10

#Plotting the fitted curves:
plot(tme, yy, type = "l")
for(i in 1:100)
{
    csp = post.samples[i,1]
    b0 = post.samples[i,2]
    b1 = post.samples[i,3]
    b2 = post.samples[i,4]
    points(tme, b0 + b1*tme + b2*(pmax(tme - csp, 0)), type = "l", col = "red")
    abline(v = csp, col = "blue")
}
summary(post.samples[,3]) #mean slope before breakpoint is 190.5
summary(post.samples[,4]) #mean slope after breakpoint is 190.5 + 38.7 = 229.2

#Sunspots Data:
sunspots.data = read.delim("SN_y_tot_V2.0_25Aug2022.txt", header = F, sep = "")
head(sunspots.data)
sunspots = sunspots.data[,1:2]
plot(sunspots[,1], sunspots[,2], xlab = "Year (1700 to 2021)", ylab = "Yearly Sunspot Numbers", type = "l", main = "Sunspot Data")

tme = 1700:2021
n = length(tme)

#Periodogram:
plot(1:(n/2), abs(fft(sunspots[,2])[2:((n/2)+1)])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)

#Bayesian posterior
grid.res = 0.001
f.val = seq(0.01, 0.49, grid.res)
X = matrix(1, nrow = n, ncol = 3)
expos = rep(-1, length(f.val)) #exact marginal posterior for omega
log.values = rep(-1, length(f.val))
log.det.term = rep(-1, length(f.val))
for(i in 1:length(f.val))
{
    X[,2] = cos(2*pi*f.val[i]*tme)
    X[,3] = sin(2*pi*f.val[i]*tme)
    mod = lm(sunspots[,2] ~ X)
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    log.det.term[i] = (0.5*(log(det(t(X) %*% X))))
    log.values[i] = log.value
}
log.values = log.values - max(log.values) #scaling to remove large values 
expos = exp(log.values)
expos = (expos/sum(expos))/grid.res
plot(f.val, expos, type = "l")

#Plotting posterior samples of the fitted function
N = 2000 #number of posterior samples
post.samples = matrix(-1, N, 5)
post.samples[,1] = sample(f.val, N, replace = T, prob = expos)
for(i in 1:N)
{
    fr = post.samples[i,1]
    X = matrix(1, nrow = n, ncol = 3)
    X[,2] = cos(2*pi*fr*tme)
    X[,3] = sin(2*pi*fr*tme)
    lin.model = lm(sunspots[,2] ~ -1 + X)
    bhat = lin.model$coefficients
    sighat = sqrt((sum((lin.model$residuals)^2))/(n-3)) #this is also denoted by the Residual Standard Error
    Sigma.mat = (sighat^2)*solve(t(X) %*% X)
    chiran = (rchisq(1, df = n-3))
    beta.samples = bhat + (rmvnorm(1, sigma = Sigma.mat))/(sqrt(chiran/(n-3)))
    sig.sample = sqrt((sum((lin.model$residuals)^2))/chiran)
    post.samples[i,2:4] = beta.samples
    post.samples[i, 5] = sig.sample
}
summary(post.samples[,5]) 

#Plotting some posterior samples
plot(tme, sunspots[,2], type = "l")
for(i in 1:50)
{
    fr = post.samples[i,1]
    b0 = post.samples[i,2]
    b1 = post.samples[i, 3]
    b2 = post.samples[i, 4]
    points(tme, b0 + b1*cos(2*pi*fr*tme) + b2*sin(2*pi*fr*tme), type = "l", col = "red")                        
}

#While this analysis gives the 11 year period for the solar cycle, it is not fully satisfactory simply because the single sinusoid + noise model is not a good fit for the data. 
#Indeed, data generated from the model will not look like the observed sunspots data. 
sig.est = 52.28 #posterior mean for sigma from the model
b0 = mean(post.samples[,2])
b1 = mean(post.samples[,3])
b2 = mean(post.samples[,4])
simul.data = b0 + b1*cos(2*pi*tme) + b2*sin(2*pi*tme) + sig.est*rnorm(length(tme))
plot(tme, simul.data, type = "l")
points(tme, b0 + b1*cos(2*pi*tme) + b2*sin(2*pi*tme), col = "red", type = "l")
#The sunspots data obviously does not look like this:
par(mfrow = c(4, 4))
for(i in 1:10){plot(tme, b0 + b1*cos(2*pi*tme) + b2*sin(2*pi*tme) + sig.est*rnorm(length(tme)), type = "l", xlab = "Time", ylab = "Data")} 
plot(tme, sunspots[,2], type = "l", xlab = "Time", ylab = "Data")
for(i in 12:16){plot(tme, b0 + b1*cos(2*pi*tme) + b2*sin(2*pi*tme) + sig.est*rnorm(length(tme)), type = "l", xlab = "Time", ylab = "Data")} 
par(mfrow = c(1, 1))

#Two sinusoid recovery:
n = 500
tme = 0:(n-1)
ftrue = c(0.201, 0.401)
amp = c(2, 2, 3, 1.5)
truth = amp[1]*cos(2*pi*ftrue[1]*tme) + amp[2]*sin(2*pi*ftrue[1]*tme) + amp[3]*cos(2*pi*ftrue[2]*tme) + amp[4]*sin(2*pi*ftrue[2]*tme)
plot(truth, type = "l")

sig = 10
yy = truth + sig*rnorm(n)
plot(yy, type = "l")

#Usual periodogram: 
plot(1:(n/2), abs(fft(yy)[2:((n/2)+1)])^2/n, type = "h", ylab = "Periodogram")

#Periodogram on dense scale: 
grid.res = 0.001
f.val = seq(0.001, 0.5, grid.res)
pgram = rep(-1, length(f.val))
for(i in 1:length(f.val))
{
    pgram[i] = (((sum(yy * cos(2*pi*f.val[i]*tme)))^2) + ((sum(yy * sin(2*pi*f.val[i]*tme)))^2))/n
}
plot(f.val, pgram, type = "l")
f.val[which.max(pgram)]

#negative logarithm of unnormalized posterior function: 
logpost = function(freq)
{
    X = matrix(1, nrow = n, ncol = 5)
    X[,2] = cos(2*pi*freq[1]*tme)
    X[,3] = sin(2*pi*freq[1]*tme)
    X[,4] = cos(2*pi*freq[2]*tme)
    X[,5] = sin(2*pi*freq[2]*tme)
    mod = lm(yy ~ X)
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    return(log.value)
}
logpost(ftrue)

#Let us plot this function (this will be a 3D plot)
library(lattice)
f1.gr = seq(0.11, 0.29, len = 200)
f2.gr = seq(0.31, 0.49, len = 200)
g = expand.grid(x=f1.gr, y=f2.gr)
for(i in 1:nrow(g)) g$z[i] = logpost(c(g$x[i], g$y[i]))
g$z = exp(g$z - max(g$z))
g$z = g$z/(sum(g$z))
wireframe(z~x*y, g)
g[which.max(g$z),]

#marginalize:
#Marginal distributions: 
mar.f1 = rep(-1, length(f1.gr))
for(i in 1:length(f1.gr))
{
  ind = (g$x == f1.gr[i])
  mar.f1[i] = sum(g$z[ind])
}
par(mfrow = c(2, 1))
plot(f1.gr, mar.f1, type = "h")
mar.f2 = rep(-1, length(f2.gr))
for(i in 1:length(f2.gr))
{
  ind = (g$y == f2.gr[i])
  mar.f2[i] = sum(g$z[ind])
}
plot(f2.gr, mar.f2, type = "h")
par(mfrow = c(1, 1))







