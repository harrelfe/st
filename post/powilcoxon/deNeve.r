f <- function(z) {
  r <- exp(z)
  r * (r - z - 1)/((r - 1)^2)
}

f(.01)
f(.001)
f(.0001)
f(.00001)
f(1e-6)
f(1e-7)

r <- seq(.5, 2, length=1000)
c <- r ^ 0.6453 / (1 + r ^ 0.6453)
plot(r, c, type='l', xlab='OR')
lines(r, f(log(r)), col='red')

r <- seq(.1, 10, length=1000)
g <- function(a) max(abs(r ^ a / (1 + r ^ a) - f(log(r))))
optimize(g, c(.6, .7), tol=1e-9). # max error 0.0013
c <- r ^ 0.675 / (1 + r ^ 0.675)
plot(r, c, type='l', xlab='OR')
lines(r, f(log(r)), col='red')
