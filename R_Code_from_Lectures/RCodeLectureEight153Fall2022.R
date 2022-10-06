#Dataset
dt = c(26.6, 38.5, 34.4, 34, 31, 23.6, 100)
n = length(dt)
probdensnorm = function(thet, sigm)
{
    exp(sum(dnorm(dt, mean = thet, sd = sigm, log = T)))   
}
probdensnorm(31.35, 5)

th.gr = seq(-50, 150, 0.5)
logsig.gr = seq(-50, 50, 0.5)
g = expand.grid(x = th.gr, y = logsig.gr)
for(i in 1:nrow(g))
{
   g$normdens[i] = probdensnorm(g$x[i], exp(g$y[i]))
}
mean(g$normdens)
g[which.max(g$normdens),]

probdenslap = function(thet, sigm)
{
    logans = -sum((abs(dt - thet))/sigm) - n*(log(2*sigm))
    return(exp(logans))
}
th.gr = seq(-50, 150, 0.5)
logsig.gr = seq(-50, 50, 0.5)
g = expand.grid(x = th.gr, y = logsig.gr)
for(i in 1:nrow(g))
{
    g$normdens[i] = probdensnorm(g$x[i], exp(g$y[i]))
    g$lapdens[i] = probdenslap(g$x[i], exp(g$y[i]))
}
normev = mean(g$normdens) #ev stands for "Evidence" (I will introduce this term in the next class)
lapev = mean(g$lapdens)
c(normev, lapev)
c(normev/(normev + lapev), lapev/(normev + lapev)) #normalized to sum to one
g[which.max(g$normdens),]
g[which.max(g$lapdens),]

#Integration:
normal.integrand = function(t)
{
    ans = (sum((dt - t)^2))^(-n/2)
    return(ans)
}
res = 0.01
th = seq(-50, 150, res)
normal.integral = rep(-1, length(th))
for(i in 1:length(th))
{
    normal.integral[i] = normal.integrand(th[i])
}
normans = (gamma(n/2))*(sum(normal.integral*res))*(pi^(-n/2))*(0.5) #ignoring the (1/(4*M1*M2)) factor

laplace.integrand = function(t)
{
    ans = (sum(abs(dt - t)))^(-n)
    return(ans)
}
res = 0.01
th = seq(-50, 150, res)
laplace.integral = rep(-1, length(th))
for(i in 1:length(th))
{
    laplace.integral[i] = laplace.integrand(th[i])
}
lapans = (gamma(n))*(sum(laplace.integral*res))*(2^(-n)) #ignoring the (1/(4*M1*M2)) factor
c(normans, lapans)
c(normans/(normans + lapans), lapans/(normans + lapans)) #normalized to sum to one






