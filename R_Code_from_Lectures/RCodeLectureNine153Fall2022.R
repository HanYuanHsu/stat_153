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

#Periodogram:
plot(1:(n/2), abs(fft(dt)[2:((n/2)+1)])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)

K = 50
C = 100
logevid.noprior = rep(-1, K)
logevid = rep(-1, K)
Xin = matrix(1, n, 1)
for(kin in 1:K)
{
    Xin = cbind(Xin, cos(2*pi*(kin/n)*tme), sin(2*pi*(kin/n)*tme))
    mod = lm(dt ~ -1 + Xin)
    p = ncol(Xin)
    log.value = ((p - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(determinant(t(Xin) %*% Xin, logarithm = T)$modulus)) + ((p/2)*(log(pi))) + lgamma((n-p)/2)
    logevid.noprior[kin] = log.value
    log.value = log.value - p*(log(2*C))
    logevid[kin] = log.value
}

cbind(logevid.noprior, logevid)
plot(logevid.noprior, type = "l")
abline(v = 20)

logevid = logevid - max(logevid)
evid = exp(logevid)
evid = evid/(sum(evid))
plot(1:K, evid, type = "h")
which.max(evid)











