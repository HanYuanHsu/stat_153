#Orthogonality of Sinusoids
n = 100 #even
tme = 0:(n-1)
f1 = 7/n
cf1 = cos(2*pi*f1*tme)
sf1 = sin(2*pi*f1*tme)
plot(cf1, type = "l")
points(sf1, type = "l", col = "red")
sum(cf1*sf1)
f2 = 8/n
cf2 = cos(2*pi*f2*tme)
points(cf2, type = "l", col = "blue")
sum(cf1*cf2)

#The Discrete Fourier Transform and Periodogram:
n = 11
x = rnorm(n)
fft(x)
pgram = abs(fft(x)[1:6])^2/n
plot(pgram, type = "h")

n = 12
x = rnorm(n)
fft(x)
pgram = abs(fft(x)[1:7])^2/n
plot(pgram, type = "h")

#Main Utility of the Periodogram: picks out sinusoidal components at Fourier Coefficients magically
n = 11
f = 1/n
tme = 0:(n-1)
x = sin(2*pi*f*tme)
plot(tme, x, type = "o")
plot(1:5, abs(fft(x)[2:6])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)

f = 1/n
x = 10*sin(2*pi*f*tme) + 3*cos(2*pi*f*tme)
plot(tme, x, type = "o")
plot(1:5, abs(fft(x)[2:6])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)

#Larger n
n = 100
tme = 0:(n-1)
f = 1/n
x = 10*sin(2*pi*f*tme) + 3*cos(2*pi*f*tme)
plot(tme, x, type = "o")
plot(1:50, abs(fft(x)[2:51])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)

n = 100
tme = 0:(n-1)
f = 6/n
x = 10*sin(2*pi*f*tme) + 3*cos(2*pi*f*tme)
plot(tme, x, type = "o")
plot(1:50, abs(fft(x)[2:51])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)

#Two Frequencies
n = 100
tme = 0:(n-1)
f1 = 0.25
f2 = 0.1
x = 10*sin(2*pi*f1*tme) + 3*cos(2*pi*f1*tme) + 5*cos(2*pi*f2*tme)
plot(tme, x, type = "o")
plot(1:50, abs(fft(x)[2:51])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)
#Once frequencies are figured out, one can do a linear regression to figure out amplitudes and phases. 
x1 = cos(2*pi*f1*tme)
x2 = sin(2*pi*f1*tme)
x3 = cos(2*pi*f2*tme)
x4 = sin(2*pi*f2*tme)
reg = lm(x ~ 1 + x1 + x2 + x3 + x4)
summary(reg)

#Five Frequencies
n = 100
tme = 0:(n-1)
f1 = 0.25
f2 = 0.1
f3 = 0.33
f4 = 0.01
f5 = 0.15
x = 10*sin(2*pi*f1*tme) + 3*cos(2*pi*f1*tme) + 5*cos(2*pi*f2*tme) + 8*sin(2*pi*f3*tme) + 4*cos(2*pi*f4*tme) + 7*sin(2*pi*f5*tme) + 3*cos(2*pi*f5*tme)
plot(tme, x, type = "o")
plot(1:50, abs(fft(x)[2:51])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)

#Structure + White Noise
n = 100
tme = 0:(n-1)
f1 = 0.25
f2 = 0.1
x = 10*sin(2*pi*f1*tme) + 3*cos(2*pi*f1*tme) + 5*cos(2*pi*f2*tme)
sig = 10
x = x + sig*rnorm(n)
plot(tme, x, type = "o")
plot(1:50, abs(fft(x)[2:51])^2/n, ylab = "Periodogram", type = "h")
abline(h = 0)
#Plot on log scale
plot(1:50, log(abs(fft(x)[2:51])^2/n), ylab = "Periodogram on Log Scale", type = "h")
abline(h = 0)

#Noisy Data:
n = 1000
tme = 0:(n-1)
f1 = 0.25
f2 = 0.1
f3 = 0.33
f4 = 0.01
f5 = 0.15
x = 10*sin(2*pi*f1*tme) + 3*cos(2*pi*f1*tme) + 5*cos(2*pi*f2*tme) + 8*sin(2*pi*f3*tme) + 4*cos(2*pi*f4*tme) + 7*sin(2*pi*f5*tme) + 3*cos(2*pi*f5*tme)
plot(tme, x, type = "l")
sig = 40
yy = x + sig*rnorm(n)
plot(tme, yy, type = "l")
points(tme, x, type = "l", col = "red")

plot((1:(n/2))/n, abs(fft(yy)[2:((n/2)+1)])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)

#Single Sinusoid at Non-Fourier Frequency: Leakage
n = 100
tme = 0:(n-1)
f = 0.06
x = 10*sin(2*pi*f*tme) + 3*cos(2*pi*f*tme)
plot(tme, x, type = "o")
plot(1:50, abs(fft(x)[2:51])^2/n, type = "h", ylab = "Periodogram") 

plot(1:50, log(abs(fft(x)[2:51])^2/n), type = "h", ylab = "Log(Periodogram)") #logarithmic scale
abline(h = 0)

#Plotting the periodogram on a more dense scale: 
grid.res = 0.001
f.val = seq(0.001, 0.5, grid.res)
pgram = rep(-1, length(f.val))
for(i in 1:length(f.val))
{
    pgram[i] = (((sum(x * cos(2*pi*f.val[i]*tme)))^2) + ((sum(x * sin(2*pi*f.val[i]*tme)))^2))/n
}
plot(f.val, pgram, type = "l")
f.val[which.max(pgram)]

#Leakage is eliminated if we increase n so as to obtain a Fourier Frequency
n = 1000
tme = 0:(n-1)
f = 0.065 
x = 10*sin(2*pi*f*tme) + 3*cos(2*pi*f*tme)
plot(tme, x, type = "o")
plot(1:500, abs(fft(x)[2:501])^2/n, type = "h", ylab = "Periodogram") #no more leakage as we now have a Fourier Frequency
abline(h = 0)

#More Leakage
n = 100
tme = 0:(n-1)
f1 = 0.253
f2 = 0.10
x = 10*sin(2*pi*f1*tme) + 3*cos(2*pi*f1*tme) + 5*cos(2*pi*f2*tme)
plot(tme, x, type = "o")
plot(1:50, abs(fft(x)[2:51])^2/n, type = "h", ylab = 'Periodogram')
abline(h = 0)

#White Noise
x = rnorm(n)
plot(x, type = "o", main = "White Noise")
plot(1:50, abs(fft(x)[2:51])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)

#Single sinusoid at Fourier frequency and periodogram
n = 100
tme = 0:(n-1)
f = 5/n
x = 10*sin(2*pi*f*tme) + 3*cos(2*pi*f*tme)
plot(tme, x, type = "o")
plot(1:50, abs(fft(x)[2:51])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)

#Plotting a denser periodogram at other frequencies as well:  
grid.res = 0.001
f.val = seq(1/n, 0.5, grid.res)
pgram = rep(-1, length(f.val))
for(i in 1:length(f.val))
{
    pgram[i] = (((sum(x * cos(2*pi*f.val[i]*tme)))^2) + ((sum(x * sin(2*pi*f.val[i]*tme)))^2))/n
}
plot(f.val, pgram, type = "l") #this is also related to leakage
abline(v = 5/n, col = "red")
abline(v = 6/n, col = "red")
abline(v = 7/n, col = "red")
abline(v = 4/n, col = "red")
abline(v = 3/n, col = "red")
plot(f.val, log(pgram), type = "l")
f.val[which.max(pgram)]

#Bayesian Posterior for Single Sinusoid at Non-Fourier Frequency
n = 100
tme = 0:(n-1)
f = 0.065
x = 10*sin(2*pi*f*tme) + 3*cos(2*pi*f*tme)

n = length(tme)
grid.res = 0.001
f.val = seq(1/n, 0.5, grid.res)
X = matrix(1, nrow = n, ncol = 3)
expos = rep(-1, length(f.val)) #exact marginal posterior for omega
log.values = rep(-1, length(f.val))
log.det.term = rep(-1, length(f.val))
for(i in 1:length(f.val))
{
    X[,2] = cos(2*pi*f.val[i]*tme)
    X[,3] = sin(2*pi*f.val[i]*tme)
    mod = lm(x ~ X)
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    log.det.term[i] = (0.5*(log(det(t(X) %*% X))))
    log.values[i] = log.value
}
log.values = log.values - max(log.values) #scaling to remove large values 
expos = exp(log.values)
expos = (expos/sum(expos))/grid.res
plot(f.val, expos, type = "l")
#The posterior is quite peaked. We can get a point estimate by the frequency which maximizes the posterior: 
ind.max = which.max(expos)
est.freq = f.val[ind.max] #posterior maximizer
est.freq
#Posterior mean: 
pm = (sum(f.val*expos))*grid.res
pm
#Posterior standard deviation:
psd = sqrt((sum(((f.val - pm)^2)*expos))*grid.res)
psd
c(pm - 2*psd, pm + 2*psd)

#Add noise:
n = 999 #n = 999
tme = 0:(n-1)
f = 0.065
x = 10*sin(2*pi*f*tme) + 3*cos(2*pi*f*tme)
sig = 50
y = x + sig*rnorm(n)

plot(y, type = "l")
points(x, type = "l", col = "red")

n = length(tme)
grid.res = 0.001
f.val = seq(1/n, 0.49, grid.res)
X = matrix(1, nrow = n, ncol = 3)
expos = rep(-1, length(f.val)) #exact marginal posterior for omega
log.values = rep(-1, length(f.val))
log.det.term = rep(-1, length(f.val))
for(i in 1:length(f.val))
{
    X[,2] = cos(2*pi*f.val[i]*tme)
    X[,3] = sin(2*pi*f.val[i]*tme)
    mod = lm(y ~ X)
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    log.det.term[i] = (0.5*(log(det(t(X) %*% X))))
    log.values[i] = log.value
}
log.values = log.values - max(log.values) #scaling to remove large values 
expos = exp(log.values)
expos = (expos/sum(expos))/grid.res
plot(f.val, expos, type = "l")
#The posterior is quite peaked. We can get a point estimate by the frequency which maximizes the posterior: 
ind.max = which.max(expos)
est.freq = f.val[ind.max] #posterior maximizer
est.freq
#Posterior mean: 
pm = (sum(f.val*expos))*grid.res
pm
#Posterior standard deviation:
psd = sqrt((sum(((f.val - pm)^2)*expos))*grid.res)
psd
c(pm - 2*psd, pm + 2*psd)
