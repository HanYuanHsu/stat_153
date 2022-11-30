library(Sleuth3)
dt = ex1029
help(ex1029)
names(dt)
n = 500
dt = dt[sample(nrow(dt), n),]
y = log(dt$WeeklyEarnings)
x = dt$Exper
plot(x, y, xlab = "Experience in years", ylab = "log(Earnings)")
n = length(y)
kmax = 63
X = matrix(1, n, (kmax+2))
X[,2] = x
for(j in 1:kmax)
{
  X[,(j+2)] = pmax(x-j,0)
}
md = lm(y ~ -1 + X)
summary(md)

#plot(x, y, type = "n")
plot(x, y, xlab = "Experience in years", ylab = "log(Earnings)")
points(x[order(x)], md$fitted.values[order(x)], col = "red", type = "l")

#Posterior density of beta (given tau and sigma)
C = 10^6
tau = 0.0001
T = diag(c(C, C, rep(tau^2, kmax)))
sig = 0.6
TempMat = solve((sig^2)*(solve(T)) + (t(X) %*% X))
pSigmat = (sig^2)*TempMat
pm = TempMat%*%(t(X)%*%(matrix(y,n,1)))
muhat = X %*% pm
points(x[order(x)], muhat[order(x)], col = "blue", type = "l")

#Marginal or integrated density calculation:
C = 10^6
library(mvtnorm)
logmarg = function(tau, sig)
{
    T = diag(c(C, C, rep(tau^2, kmax)))
    Sigmat = (X%*%T%*%(t(X))) + diag(rep(sig^2, n))
    ans = dmvnorm(y, mean = rep(0, n), sigma = Sigmat, log = TRUE) 
    return(ans)
}
#grid based minimization of logmarg
taugrid = seq(0.001, 0.1, length.out = 50)
siggrid = seq(0.1, 1, length.out = 50)
g = expand.grid(tau=taugrid, sig=siggrid)
for(i in 1:nrow(g)) {g$lm[i] = logmarg(g$tau[i], g$sig[i]); print(i)}
ind.max = which.max(g$lm)
g[ind.max,]                                      

tau = g$tau[ind.max]
sig = g$sig[ind.max]
c(tau, sig)
T = diag(c(C, C, rep(tau^2, kmax)))
TempMat = solve((sig^2)*(solve(T)) + (t(X) %*% X))
pSigmat = (sig^2)*TempMat
pm = TempMat%*%(t(X)%*%(matrix(y,n,1)))
muhat = X %*% pm
#plot(x, y, type = "n", xlab = "Experience in years", ylab = "Earnings")
plot(x, y, xlab = "Experience in years", ylab = "Earnings", col = "gray")
points(x[order(x)], md$fitted.values[order(x)], col = "red", type = "l")
points(x[order(x)], muhat[order(x)], col = "blue", type = "l")






#Google Trends for Yahoo
yahoo.raw = read.delim("YahooTrends19Nov2022.csv")
#the useful data is between row 2 and row 228
yahoo.use = yahoo.raw[2:228,1]
len = length(yahoo.use)
yahoo = rep(0, len)
for(i in 1:len)
{
  yahoo[i] = as.numeric(unlist(strsplit(as.character(yahoo.use[i]), ","))[2])
}
plot(yahoo, type = "l")
n = length(yahoo)
x = 1:n
X = matrix(1, n, n)
X[,2] = x
for(j in 2:(n-1))
{
  X[,(j+1)] = pmax(x-j,0)
}

#Posterior density of beta (given tau and sigma)
C = 10^6
tau = 0.05
T = diag(c(C, C, rep(tau^2, (n-2))))
sig = 1
TempMat = solve((sig^2)*(solve(T)) + (t(X) %*% X))
pSigmat = (sig^2)*TempMat
pm = TempMat%*%(t(X)%*%(matrix(yahoo,n,1)))
muhat = X %*% pm

plot(x, yahoo, type = "l", col = "gray")
points(x, muhat, col = "blue", type = "l")

#Marginal density calculation:
C = 10^6
y = yahoo
library(mvtnorm)
logmarg = function(tau, sig)
{
    T = diag(c(C, C, rep(tau^2, (n-2))))
    Sigmat = (X%*%T%*%(t(X))) + diag(rep(sig^2, n))
    ans = dmvnorm(y, mean = rep(0, n), sigma = Sigmat, log = TRUE) 
    return(ans)
}
#grid based minimization of logmarg
taugrid = seq(0.005, 1, length.out = 100)
siggrid = seq(1, 5, length.out = 100)
g = expand.grid(tau=taugrid, sig=siggrid)
for(i in 1:nrow(g)) {g$lm[i] = logmarg(g$tau[i], g$sig[i]); print(i)}
ind.max = which.max(g$lm)
g[ind.max,]                                      

tau = g$tau[ind.max]
sig = g$sig[ind.max]
c(tau, sig)
T = diag(c(C, C, rep(tau^2, (n-2))))
TempMat = solve((sig^2)*(solve(T)) + (t(X) %*% X))
pSigmat = (sig^2)*TempMat
pm = TempMat%*%(t(X)%*%(matrix(y,n,1)))
muhat = X %*% pm
plot(x, y,col = "gray")
points(x[order(x)], muhat[order(x)], col = "blue", type = "l")

 
#US Population Dataset
unrate.raw = read.csv("UNRATE22Nov2022.csv")
#Data downloaded from FRED. Yearly (seasonally adjusted) Data. Rate represents the number of unemployed as a percentage of the labor force
y = unrate.raw$UNRATE
plot(y, type = "l", xlab = "Year", ylab = "Unemployment Rate")
n = length(y)
x = 1:n
X = matrix(1, n, n)
X[,2] = x
for(j in 2:(n-1))
{
  X[,(j+1)] = pmax(x-j,0)
}

#Posterior density of beta (given tau and sigma)
C = 1e6
tau = 0.5
T = diag(c(C, C, rep(tau^2, (n-2))))
sig = 1
TempMat = solve((sig^2)*(solve(T)) + (t(X) %*% X))
pSigmat = (sig^2)*TempMat
pm = TempMat%*%(t(X)%*%(matrix(y,n,1)))
muhat = X %*% pm

plot(x, y, type = "l", col = "black")
points(x, muhat, col = "blue", type = "l")

#Marginal density calculation:
C = 100
library(mvtnorm)
logmarg = function(tau, sig)
{
    T = diag(c(C, C, rep(tau^2, (n-2))))
    Sigmat = (X%*%T%*%(t(X))) + diag(rep(sig^2, n))
    ans = dmvnorm(y, mean = rep(0, n), sigma = Sigmat, log = TRUE) 
    return(ans)
}
#grid based minimization of logmarg
taugrid = seq(0.02, 0.1, length.out = 10)
siggrid = seq(0.1, 3, length.out = 20)
g = expand.grid(tau=taugrid, sig=siggrid)
for(i in 1:nrow(g)) {g$lm[i] = logmarg(g$tau[i], g$sig[i]); print(i)}
ind.max = which.max(g$lm)
g[ind.max,]                                      

tau = g$tau[ind.max]
sig = g$sig[ind.max]
c(tau, sig)
T = diag(c(C, C, rep(tau^2, (n-2))))
TempMat = solve((sig^2)*(solve(T)) + (t(X) %*% X))
pSigmat = (sig^2)*TempMat
pm = TempMat%*%(t(X)%*%(matrix(y,n,1)))
muhat = X %*% pm
plot(x, y,col = "gray", type = "l")
points(x[order(x)], muhat[order(x)], col = "blue", type = "l")


