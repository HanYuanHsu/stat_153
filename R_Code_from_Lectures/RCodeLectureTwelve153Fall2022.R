#Model selection of the number of sinuosoids:
#Example (two sinusoids)
n = 500
tme = 0:(n-1)
ftrue = c(0.201, 0.401)
amp = c(2, 2, 3, 1.5)
truth = amp[1]*cos(2*pi*ftrue[1]*tme) + amp[2]*sin(2*pi*ftrue[1]*tme) + amp[3]*cos(2*pi*ftrue[2]*tme) + amp[4]*sin(2*pi*ftrue[2]*tme)
plot(truth, type = "l")

sig = 3
yy = truth + sig*rnorm(n)
plot(yy, type = "l")

#Usual periodogram: 
perdgrm = abs(fft(yy)[2:((n/2)+1)])^2/n
plot(1:(n/2), perdgrm, type = "h", ylab = "Periodogram")

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

#Least squares minimization: 
squares = function(freq)
{
    k = length(freq)
    X = matrix(1, nrow = n, ncol = 1+(2*k))
    for(j in 1:k)
    {
        X[,(2*j)] = cos(2*pi*freq[j]*tme)
        X[,(1+2*j)] = sin(2*pi*freq[j]*tme)
    }
    mod = lm(yy ~ X)
    return(sum(mod$residuals^2))
}
#Optim minimizes functions iteratively starting from an initializer:
fhat = optim(ftrue, squares)$par
fhat
optim(fhat, squares) #a value of 0 for "convergence" indicates successful completion of the algorithm.

#Grid based calculation of the least squares minimizer: 
f1.gr = seq(0.11, 0.29, len = 200)
f2.gr = seq(0.31, 0.49, len = 200)
g = expand.grid(x=f1.gr, y=f2.gr)
for(i in 1:nrow(g)) g$sq[i] = squares(c(g$x[i], g$y[i]))
g[which.min(g$sq),]

fhatgridinit = optim(g[which.min(g$sq),1:2], squares)$par
fhatgridinit
#Compare to:
optim(ftrue, squares)$par

#Hessian calculation
library(numDeriv)
H = hessian(squares, fhatgridinit)
log(det(H))

#Evidence calculation:
k = 2 #number of sinusoids
p = 2*k + 1 #rank of the X matrix

Xfhat = matrix(1, nrow = n, ncol = 1+(2*k))
freq = fhatgridinit
for(j in 1:k)
{
   Xfhat[,(2*j)] = cos(2*pi*freq[j]*tme)
   Xfhat[,(1+2*j)] = sin(2*pi*freq[j]*tme)
}
mod = lm(yy ~ Xfhat)

t1 = lgamma(p/2) - (p/2)*(log(sum(mod$fitted.values^2)))
t2 = lgamma((n - p - k - 1)/2) - ((n-p-k)/2)*(log(sum(mod$residuals^2)))
t3 = lgamma(k/2) - (k/2)*(log(sum(fhatgridinit^2)))
t4 = (-0.5)*(log(det(0.5*H)))
log.evid.2 = t1 + t2 + t3 + t4
log.evid.2

#Evidence for model with one sinusoid: 
grid.res = 0.001
f.val = seq(0.01, 0.49, grid.res)
lsobj = lapply(f.val, squares)
plot(f.val, lsobj, type = "l")
fhat = f.val[which.min(lsobj)]
fhat
#Hessian calculation
H = hessian(squares, fhat)
log(det(H))

#Evidence calculation:
k = 1 #number of sinusoids
p = 2*k + 1 #rank of the X matrix

Xfhat = matrix(1, nrow = n, ncol = 1+(2*k))
freq = fhat
for(j in 1:k)
{
   Xfhat[,(2*j)] = cos(2*pi*freq[j]*tme)
   Xfhat[,(1+2*j)] = sin(2*pi*freq[j]*tme)
}
mod = lm(yy ~ Xfhat)

t1 = lgamma(p/2) - (p/2)*(log(sum(mod$fitted.values^2)))
t2 = lgamma((n - p - k - 1)/2) - ((n-p-k)/2)*(log(sum(mod$residuals^2)))
pt3 = lgamma(k/2) - (k/2)*(log(sum(fhat^2)))
t4 = (-0.5)*(log(det(0.5*H)))
log.evid.1 = t1 + t2 + t3 + t4

c(log.evid.1, log.evid.2)


#Sunspots Data:
sunspots.data = read.delim("SN_y_tot_V2.0_25Aug2022.txt", header = F, sep = "")
head(sunspots.data)
sunspots = sunspots.data[,1:2]
plot(sunspots[,1], sunspots[,2], xlab = "Year (1700 to 2021)", ylab = "Yearly Sunspot Numbers", type = "l", main = "Sunspot Data")

tme = 1700:2021
n = length(tme)

#Periodogram:
perdgrm = abs(fft(sunspots[,2])[2:((n/2)+1)])^2/n
plot(1:(n/2), perdgrm, type = "h", ylab = "Periodogram")
abline(h = 0)
#top three Fourier frequencies:
(which(perdgrm > 50000))/n

#Calculate Evidence for model with one sinusoid:
yy = sunspots[,2]
squares = function(freq)
{
    k = length(freq)
    X = matrix(1, nrow = n, ncol = 1+(2*k))
    for(j in 1:k)
    {
        X[,(2*j)] = cos(2*pi*freq[j]*tme)
        X[,(1+2*j)] = sin(2*pi*freq[j]*tme)
    }
    mod = lm(yy ~ X)
    return(sum(mod$residuals^2))
}

grid.res = 0.001
f.val = seq(0.01, 0.49, grid.res)
lsobj = lapply(f.val, squares)
plot(f.val, lsobj, type = "l")
fhat = f.val[which.min(lsobj)]
fhat
fhat = optim(fhat, squares)$par
fhat

#Hessian calculation
H = hessian(squares, fhat)
log(det(H))
k = 1 #number of sinusoids
p = 2*k + 1 #rank of the X matrix

Xfhat = matrix(1, nrow = n, ncol = 1+(2*k))
freq = fhat
for(j in 1:k)
{
   Xfhat[,(2*j)] = cos(2*pi*freq[j]*tme)
   Xfhat[,(1+2*j)] = sin(2*pi*freq[j]*tme)
}
mod = lm(yy ~ Xfhat)

t1 = lgamma(p/2) - (p/2)*(log(sum(mod$fitted.values^2)))
t2 = lgamma((n - p - k - 1)/2) - ((n-p-k)/2)*(log(sum(mod$residuals^2)))
t3 = lgamma(k/2) - (k/2)*(log(sum(fhat^2)))
t4 = (-0.5)*(log(det(0.5*H)))
log.evid.1 = t1 + t2 + t3 + t4
log.evid.1

#Two sinusoids
cbind(which(perdgrm > 50000)/n, perdgrm[which(perdgrm > 50000)])
f1.gr = seq(0.085, 0.095, len = 200)
f2.gr = seq(0.096, 0.105, len = 200)
g = expand.grid(x=f1.gr, y=f2.gr)
for(i in 1:nrow(g)) g$sq[i] = squares(c(g$x[i], g$y[i]))
g[which.min(g$sq),]
fhatgridinit = optim(g[which.min(g$sq),1:2], squares)$par
fhatgridinit
1/fhatgridinit #these correspond to the periods of the two cycles
squares(fhatgridinit)

#Top two periodogram maximizers:
k = 2
f2 = ((order(perdgrm, decreasing = T))[1:k])/n
fhatf2 = optim(f2, squares)$par
c(squares(fhatf2), squares(fhatgridinit))

#Hessian calculation
library(numDeriv)
H = hessian(squares, fhatgridinit)
log(det(H))

#Evidence calculation:
k = 2 #number of sinusoids
p = 2*k + 1 #rank of the X matrix

Xfhat = matrix(1, nrow = n, ncol = 1+(2*k))
freq = fhatgridinit
for(j in 1:k)
{
   Xfhat[,(2*j)] = cos(2*pi*freq[j]*tme)
   Xfhat[,(1+2*j)] = sin(2*pi*freq[j]*tme)
}
mod = lm(yy ~ Xfhat)

t1 = lgamma(p/2) - (p/2)*(log(sum(mod$fitted.values^2)))
t2 = lgamma((n - p - k - 1)/2) - ((n-p-k)/2)*(log(sum(mod$residuals^2)))
t3 = lgamma(k/2) - (k/2)*(log(sum(fhatgridinit^2)))
t4 = (-0.5)*(log(det(0.5*H)))
log.evid.2 = t1 + t2 + t3 + t4
log.evid.2

c(log.evid.1, log.evid.2)

#Three Sinusoids
k = 3
f3 = ((order(perdgrm, decreasing = T))[1:k])/n
fhatf3 = optim(f3, squares)$par
c(squares(fhatf2), squares(fhatf3))

1/fhatf3 #these correspond to the periods of the three cycles
H = hessian(squares, fhatf3)
log(det(H))
k = 3 #number of sinusoids
p = 2*k + 1 #rank of the X matrix

Xfhat = matrix(1, nrow = n, ncol = 1+(2*k))
freq = fhatf3
for(j in 1:k)
{
   Xfhat[,(2*j)] = cos(2*pi*freq[j]*tme)
   Xfhat[,(1+2*j)] = sin(2*pi*freq[j]*tme)
}
mod = lm(yy ~ Xfhat)

t1 = lgamma(p/2) - (p/2)*(log(sum(mod$fitted.values^2)))
t2 = lgamma((n - p - k - 1)/2) - ((n-p-k)/2)*(log(sum(mod$residuals^2)))
t3 = lgamma(k/2) - (k/2)*(log(sum(fhat^2)))
t4 = (-0.5)*(log(det(0.5*H)))
log.evid.3 = t1 + t2 + t3 + t4
log.evid.3

logevid = c(log.evid.1, log.evid.2, log.evid.3)
logevid

#Evidences:
evid = exp(logevid - max(logevid))
evid = evid/sum(evid)
evid







