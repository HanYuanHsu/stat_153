#Simulated Data for Illustrating Model Selection
#The problem involves figuring out the right order k for fitting k sinusoids (at frequencies 1/n,...,k/n) to the observed time series
#Generating the data:
n = 5000
ktrue = 20
sigtrue = 75
betatrue = c(0, rep(5, ktrue), rep(5, ktrue))
tme = 1:n
dt = rep(-1, n)
Xtrue = matrix(1, n, 1+(2*ktrue))
for(j in 1:ktrue)
{
    Xtrue[,(2*j)] = cos(2*pi*(j/n)*tme)
    Xtrue[,((2*j)+1)] = sin(2*pi*(j/n)*tme)
}
dt = Xtrue %*% betatrue + sigtrue*rnorm(n)
plot(dt, type = "l")

#Periodogram
plot(1:(n/2), abs(fft(dt)[2:((n/2)+1)])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)
#Note that the periodogram ordinates for the initial frequencies are quite high. 

#Let us calculate the Evidences for each model (for k = 1,...,50) using the formula from the class
#The prior here is Unif(-C, C) for each of the parameters
#The key quantity is C and the answer crucially depends on C
K = 50
C = 25  #play around with various values of C and convince yourself that the answers change quite a bit depending on the specific value of C
logevid = rep(-1, K)
Xin = matrix(1, n, 1)
for(kin in 1:K)
{
    Xin = cbind(Xin, cos(2*pi*(kin/n)*tme), sin(2*pi*(kin/n)*tme))
    mod = lm(dt ~ -1 + Xin)
    p = ncol(Xin)
    log.value = ((p - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(determinant(t(Xin) %*% Xin, logarithm = T)$modulus)) + ((p/2)*(log(pi))) + lgamma((n-p)/2)
    log.value = log.value - p*(log(2*C))
    logevid[kin] = log.value
}
logevid = logevid - max(logevid)
evid = exp(logevid)
evid = evid/(sum(evid))
plot(1:K, evid, type = "h")
which.max(evid)
evid[which.max(evid)]

#Another way of doing model selection is via the standard criteria AIC and BIC
#BIC and AIC Calculation (using formulae given in Lecture 10)
bic = rep(-9999, K)
aic = rep(-9999, K)
biclogevid = rep(-9999, K)
Xin = matrix(1, n, 1)
for(kin in 1:K)
{
    Xin = cbind(Xin, cos(2*pi*(kin/n)*tme), sin(2*pi*(kin/n)*tme))
    mod = lm(dt ~ -1 + Xin)
    p = ncol(Xin)
    sse = sum(mod$residuals^2) #this is also known as the residual sum of squares
    bic[kin] = n + n*log(2*pi*sse/n) + p*(log(n))
    aic[kin] = n + n*log(2*pi*sse/n) + p*(2)
    biclogevid[kin] = (-0.5)*bic[kin]
}
biclogevid = biclogevid - max(biclogevid)
bicevid = exp(biclogevid)
bicevid = bicevid/(sum(bicevid))
plot(1:K, bicevid, type = "h")
which.max(bicevid)
bicevid[which.max(bicevid)]

which.min(aic)
plot(c(1:K, 1:K), c(aic, bic), type = "n")
points(1:K, aic, type = "l")
points(1:K, bic, type = "l", col = "red")
#Generally aic gives a slightly larger model compared to the bic


#We can calculate the Evidence using a slightly different formula (where some of the p's are changed to (p-1)). This formula was given in class. The result of this calculation will generally be almost the same as the previous formula. This formula depends on the prior only through the prior density at the MLE.  
K = 50
C = 25
logevidalt = rep(-1, K)
Xin = matrix(1, n, 1)
for(kin in 1:K)
{
    Xin = cbind(Xin, cos(2*pi*(kin/n)*tme), sin(2*pi*(kin/n)*tme))
    mod = lm(dt ~ -1 + Xin)
    p = ncol(Xin)
    log.value = ((p - 1 - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(determinant(t(Xin) %*% Xin, logarithm = T)$modulus)) + ((p/2)*(log(pi))) + lgamma((n-p-1)/2)
    log.value = log.value - p*(log(2*C))
    logevidalt[kin] = log.value
}
logevidalt = logevidalt - max(logevidalt)
evidalt = exp(logevidalt)
evidalt = evidalt/(sum(evidalt))
plot(1:K, evidalt, type = "h")
which.max(evidalt)
evidalt[which.max(evidalt)]

plot(c(1:K, 1:K), c(evid, evidalt), type = "n")
points(1:K, evid, type = "l")
points(1:K, evidalt, type = "l", col = "blue")
cbind(evid, evidalt) #The two formulae for the Evidences give almost the same answers. 

#The most important question now is: How to choose C? 
#The alternative evidence formula depends on the prior only through the value of the prior density at the MLE. 
#Since the value of the density of Unif(-C, C) is the same as that of the density Unif(a - C, a + C) for any number a,
#we can think of the prior density of beta_j as being Unif(betahat_j - C, betahat_j + C) 
#where betahat is the least squares estimate in the model. 
#Now the rule of thumb for choosing C is as follows: 
#it should be large enough to cover the region of concentration of the likelihood. 
#It should not be too large compared to the region of concentration of the likelihood 
#(otherwise, we will be paying attention to irrelevant values of the parameters). 
#To illustrate this in the current example, consider, say, the model with k = K = 50 sinusoids that we fitted above:
summary(mod)
#The standard errors of the beta parameters in this fitted model are about 1.5. 
#This means that the likelihood will be very small beyond the range betahat + [-C, C]^p 
#when C is a small constant (such as 10) multiple of 1.5. 
#A reasonable choice of C is therefore 15. 
#One can check that for a fairly large set of values of C near 15, 
#the evidence for the true model (k = 20) remains quite high. 
K = 50
C = 15
logevidalt = rep(-1, K)
Xin = matrix(1, n, 1)
for(kin in 1:K)
{
    Xin = cbind(Xin, cos(2*pi*(kin/n)*tme), sin(2*pi*(kin/n)*tme))
    mod = lm(dt ~ -1 + Xin)
    p = ncol(Xin)
    log.value = ((p - 1 - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(determinant(t(Xin) %*% Xin, logarithm = T)$modulus)) + ((p/2)*(log(pi))) + lgamma((n-p-1)/2)
    log.value = log.value - p*(log(2*C))
    logevidalt[kin] = log.value
}
logevidalt = logevidalt - max(logevidalt)
evidalt = exp(logevidalt)
evidalt = evidalt/(sum(evidalt))
plot(1:K, evidalt, type = "h")
which.max(evidalt)
evidalt[which.max(evidalt)]
#If one bumps up C to be much larger than 15, then the evidence will be high for smaller models; such choices of C are not ideal as they focus on irrelevant regions of the parameter space. On the other hand, if C is much smaller than 15, then larger models will be preferred but this is less of an issue as small C will not be used as it indicates quite an informative prior. 

#Let us now look at a real dataset:
#US Accidental Deaths:
dt = USAccDeaths
n = length(dt)
tme = 1:n
plot(tme, dt, type = "l")
#This is a monthly dataset with a seasonal pattern. To model it, we can use sinusoids with a period of d=12: cos(2*pi*(j/d)*t) and sin(2*pi*(j/d)*t) for j = 1, 2, 3 etc. Before trying this however, let us first fit sinusoids at Fourier frequencies j/n for j = 1, 2, etc. Note that here n = 72 and j = 6 corresponds to j/n = 1/d

#To figure out a suitable value of C, let us first fit one of the candidate models 
#and look at the standard errors of the coefficients:
mod = lm(dt ~ cos(2*pi*(1/n)*tme) + sin(2*pi*(1/n)*tme) + cos(2*pi*(2/n)*tme) + sin(2*pi*(2/n)*tme) + cos(2*pi*(6/n)*tme) + sin(2*pi*(6/n)*tme))
summary(mod)
#the standard errors are around 100. A reasonable value of C is thus: 
C = 600
K = 20
logevidalt = rep(-1, K)
Xin = matrix(1, n, 1)
for(kin in 1:K)
{
    Xin = cbind(Xin, cos(2*pi*(kin/n)*tme), sin(2*pi*(kin/n)*tme))
    mod = lm(dt ~ -1 + Xin)
    p = ncol(Xin)
    log.value = ((p - 1 - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(determinant(t(Xin) %*% Xin, logarithm = T)$modulus)) + ((p/2)*(log(pi))) + lgamma((n-p-1)/2)
    log.value = log.value - p*(log(2*C))
    logevidalt[kin] = log.value
}
logevidalt = logevidalt - max(logevidalt)
evidalt = exp(logevidalt)
evidalt = evidalt/(sum(evidalt))
plot(1:K, evidalt, type = "h")
which.max(evidalt)
evidalt[which.max(evidalt)] #the best model is the model with k = 6.  

#BIC and AIC Calculation
bic = rep(-9999, K)
aic = rep(-9999, K)
biclogevid = rep(-9999, K)
Xin = matrix(1, n, 1)
for(kin in 1:K)
{
    Xin = cbind(Xin, cos(2*pi*(kin/n)*tme), sin(2*pi*(kin/n)*tme))
    mod = lm(dt ~ -1 + Xin)
    p = ncol(Xin)
    sse = sum(mod$residuals^2)
    bic[kin] = n + n*log(2*pi*sse/n) + p*(log(n))
    aic[kin] = n + n*log(2*pi*sse/n) + p*(2)
    biclogevid[kin] = (-0.5)*bic[kin]
}
biclogevid = biclogevid - max(biclogevid)
bicevid = exp(biclogevid)
bicevid = bicevid/(sum(bicevid))
plot(1:K, bicevid, type = "h")
which.max(bicevid)
bicevid[which.max(bicevid)]

which.min(aic)
plot(c(1:K, 1:K), c(aic, bic), type = "n")
points(1:K, aic, type = "l")
points(1:K, bic, type = "l", col = "red")

#Now instead of looking at all Fourier frequencies 1/n, 2/n,.., k/n, let us just consider models with sinusoidal components at frequencies 1/d, 2/d,.., k/d where d = 12. We shall consider up to k = 5 (note that, at k = 6, the sin function will be zero).  
K = 5
C = 1000
logevidalt = rep(-1, K)
rse = rep(-1, K) #this has the residual standard errors
betamin = rep(-1, K)
betamax = rep(-1, K)
betahatmat = matrix(-99999, K, 2+((2*K) + 1))
Xin = matrix(1, n, 1)
d = 12
for(kin in 1:K)
{
    Xin = cbind(Xin, cos(2*pi*(kin/d)*tme), sin(2*pi*(kin/d)*tme))
    mod = lm(dt ~ -1 + Xin)
    p = ncol(Xin)
    rse[kin] = sqrt((sum(mod$residuals^2))/(n-p))
    betamin[kin] = min(mod$coefficients)
    betamax[kin] = max(mod$coefficients)
    betahatmat[kin, 1:p] = mod$coefficients
    log.value = ((p - 1 - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(determinant(t(Xin) %*% Xin, logarithm = T)$modulus)) + ((p/2)*(log(pi))) + lgamma((n-p-1)/2)
    log.value = log.value - p*(log(2*C))
    logevidalt[kin] = log.value
}
logevidalt = logevidalt - max(logevidalt)
evidalt = exp(logevidalt)
evidalt = evidalt/(sum(evidalt))
plot(1:K, evidalt, type = "h")
which.max(evidalt)
evidalt[which.max(evidalt)]
evidalt #depending on the value of C (try C = 500 and C = 1000), the evidence prefers model 3 or model 2.
#Let us see what BIC and AIC give us here
#BIC and AIC Calculation
K = 5
d = 12
bic = rep(-9999, K)
aic = rep(-9999, K)
biclogevid = rep(-9999, K)
Xin = matrix(1, n, 1)
for(kin in 1:K)
{
    Xin = cbind(Xin, cos(2*pi*(kin/d)*tme), sin(2*pi*(kin/d)*tme))
    mod = lm(dt ~ -1 + Xin)
    p = ncol(Xin)
    sse = sum(mod$residuals^2)
    bic[kin] = n + n*log(2*pi*sse/n) + p*(log(n))
    aic[kin] = n + n*log(2*pi*sse/n) + p*(2)
}
cbind(aic, bic)
which.min(bic)
which.min(aic)

#Another real dataset
#Linear/quadratic plus seasonal trend : Google Trends (monthly data) for the query "Amazon" (we used this dataset previously in Lecture 3)
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

#We looked at four models for this dataset previously (see code for Lecture 3)
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

#Model Selection among the four models: 
#Model One:
n = length(amazon)
dt = amazon
tme = 1:n
X1 = matrix(1, n, 1)
X1 = cbind(X1, tme)
f1 = 1
d = 12
X1 = cbind(X1, cos(2*pi*f1*tme/d), sin(2*pi*f1*tme/d))
#Evidence for Model One:
mod1 = lm(dt ~ -1 + X1)
summary(mod1)
C1 = 10
p = ncol(X1)
log.value = ((p - 1 - n)/2)*(log(sum(mod1$residuals^2))) - (0.5*(determinant(t(X1) %*% X1, logarithm = T)$modulus)) + ((p/2)*(log(pi))) + lgamma((n-p-1)/2)
logev1 = log.value - p*(log(2*C1))
sse1 = sum(mod1$residuals^2)
bic1 = n + n*log(2*pi*sse1/n) + p*(log(n))
aic1 = n + n*log(2*pi*sse1/n) + p*(2)

#Model Two:
f2 = 2
X2 = cbind(X1, cos(2*pi*f2*tme/d), sin(2*pi*f2*tme/d))
#Evidence for Model One:
mod2 = lm(dt ~ -1 + X2)
summary(mod2)
C2 = 10
p = ncol(X2)
log.value = ((p - 1 - n)/2)*(log(sum(mod2$residuals^2))) - (0.5*(determinant(t(X2) %*% X2, logarithm = T)$modulus)) + ((p/2)*(log(pi))) + lgamma((n-p-1)/2)
logev2 = log.value - p*(log(2*C2))
sse2 = sum(mod2$residuals^2)
bic2 = n + n*log(2*pi*sse2/n) + p*(log(n))
aic2 = n + n*log(2*pi*sse2/n) + p*(2)


#Model Three: quadratic trend plus two seasonal terms:
X3 = cbind(X2, tme^2)
#Evidence for Model One:
mod3 = lm(dt ~ -1 + X3)
summary(mod3)
C3 = 10
p = ncol(X3)
log.value = ((p - 1 - n)/2)*(log(sum(mod3$residuals^2))) - (0.5*(determinant(t(X3) %*% X3, logarithm = T)$modulus)) + ((p/2)*(log(pi))) + lgamma((n-p-1)/2)
logev3 = log.value - p*(log(2*C3))
sse3 = sum(mod3$residuals^2)
bic3 = n + n*log(2*pi*sse3/n) + p*(log(n))
aic3 = n + n*log(2*pi*sse3/n) + p*(2)


#Model Four: quadratic trend plus four seasonal terms:
f3 = 3
X4 = cbind(X3, cos(2*pi*f3*tme/d), sin(2*pi*f3*tme/d))
mod4 = lm(dt ~ -1 + X4)
summary(mod4)
C4 = 10
p = ncol(X4)
log.value = ((p - 1 - n)/2)*(log(sum(mod4$residuals^2))) - (0.5*(determinant(t(X4) %*% X4, logarithm = T)$modulus)) + ((p/2)*(log(pi))) + lgamma((n-p-1)/2)
logev4 = log.value - p*(log(2*C4))
sse4 = sum(mod4$residuals^2)
bic4 = n + n*log(2*pi*sse4/n) + p*(log(n))
aic4 = n + n*log(2*pi*sse4/n) + p*(2)

logevec = c(logev1, logev2, logev3, logev4)
logevec = logevec - max(logevec)
evec = exp(logevec)
evec = evec/(sum(evec))
evec

c(sse1, sse2, sse3, sse4)
c(bic1, bic2, bic3, bic4)
c(aic1, aic2, aic3, aic4)



