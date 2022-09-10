#FREQUENCY DOMAIN ANALYSIS
#SINUSOIDS
R = 3
f = 5
Ph = 0
bcf = function(t){R*cos(2*pi*f*t + Ph)}
plot(bcf, ylab = "Sinusoid")

Ph = pi/2
bcf = function(t){R*cos(2*pi*f*t + Ph)}
plot(bcf, ylab = "Sinuosoid")

Ph = -pi/2
bcf = function(t){R*cos(2*pi*f*t + Ph)}
plot(bcf, ylab = "Sinusoid")


#Guess sinusoid from data
t = 0:99; cos1 = cos(2*pi*t*(5/100))
plot(t, cos1, ylab = "Data")
plot(t, cos1, ylab = "Data", type = "o")

t = 0:99; cos2 = cos(2*pi*t*(10/100) + (pi/3))
plot(t, cos2, ylab = "Data", type = "o")

t = 0:99; cos2 = cos(2*pi*t*(5.5/100))
plot(t, cos2, ylab = "Data", type = "o")

t = 0:99; cos2 = cos(2*pi*t*(5.7/100))
plot(t, cos2, ylab = "Data", type = "o")
#Sinusoids with frequencies of the form (j/n) are easy to guess (here n is the data size and j is an integer). These frequencies are called Fourier frequencies. 

#Understanding the DFT
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

#Guessing Sinusoids
n = 11
tme = 0:(n-1)
f = 1/n
x = cos(2*pi*f*tme)
plot(tme, x, type = "o")

#DFT
fft(x)[1:6]
abs(fft(x)[1:6])

x = 100 + x
fft(x)[1:6]
abs(fft(x)[1:6])

#Periodogram
plot(1:5, abs(fft(x)[2:6])^2/n)
#Better plot:
plot(1:5, abs(fft(x)[2:6])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)
#Thus for the cosine of frequency 1/n, I(j/n) equals zero for all j = 2, 3, 4, 5.

f = 1/n
x = 100*sin(2*pi*f*tme)
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

