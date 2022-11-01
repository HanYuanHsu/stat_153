#Fitting MA(1) model to data:
dt = arima.sim(n = 200, list(ma = 0.75))
n = length(dt)
plot(dt, type = "l")
acf(dt, lag.max = 30, type = "correlation", plot = T, main = "Sample Autocorrelation")

mamod = arima(dt, order = c(0, 0, 1))
mamod

library(emulator) #this library has the function quad.form for evaluating quadratic forms efficiently 
loglikfun = function(thet)
{
    n = length(dt)
    vec1 = c(1+(thet^2), thet, rep(0, (n-2)))
    Gamma = toeplitz(vec1)
    invGamma = solve(Gamma)
    wts = rowSums(invGamma)
    wts = wts/(sum(wts))
    muval = sum(dt * wts)
    muvec = rep(muval, n)
    sig = sqrt((quad.form(invGamma, (dt - muvec)))/n)
    logdetval = unname(unlist(determinant(Gamma, logarithm = T))[1])
    ans = -(n/2)*log(2*pi) - (n/2)*log(sig^2) - ((0.5)*logdetval) - (n/2)
    return(ans)
}

thetvals = seq(0, 2, 0.001)
M = length(thetvals)
loglikvals = rep(-9999, M)
for(i in 1:M)
{
   loglikvals[i] = loglikfun(thetvals[i])
}
plot(thetvals, loglikvals, type = "l")
thetvals[which.max(loglikvals)]

thethat = (optim(0, loglikfun, control = list(fnscale = -1)))$par
vec1 = c(1+(thethat^2), thethat, rep(0, (n-2)))
Gamma = toeplitz(vec1)
invGamma = solve(Gamma)
wts = rowSums(invGamma)
wts = wts/(sum(wts))
muhat = sum(dt * wts)
muvec = rep(muhat, n)
sig = sqrt((quad.form(invGamma, (dt - muvec)))/n)
c(muhat, sig^2)

#The Sum of Squares approach and Uncertainty Quantification
Sfunc = function(alpha)
{
    mu = alpha[1]
    thet = alpha[-1]
    n = length(dt)
    cm2 = mu
    cm1 = 0
    cm = cm1 + cm2 #this is the current conditional mean
    ans = (dt[1] - cm)^2
    for(t in 2:n)
    {
        cm1 = (-thet)*cm1 + thet*dt[t-1]
        cm2 = cm2 + (mu*((-thet)^(t-1)))
        cm = cm1 + cm2
        ans = ans + ((dt[t] - cm)^2)
    }
    return(ans)
}

alphaest = optim(c(0, 0), Sfunc)$par
alphaest

#Hessian calculation
library(numDeriv)
H = hessian(Sfunc, alphaest)
sighat = sqrt(Sfunc(alphaest)/(n-2))
sighat
covmat = (sighat^2)*(solve(0.5*H))
covmat
stderrs = sqrt(diag(covmat))
stderrs









