# Find a random permutation that makes beta for x have nearly zero coefficient
require(rms)
y <- 1:100
beta <- numeric(1000)
for(i in 1:1000) {
  set.seed(i)
  x <- sample(1:100)
  co <- coef(orm(y ~ x))
  beta[i] <- co[length(co)]
}
i <- which(abs(beta) == min(abs(beta)))
i
beta[i]
set.seed(i)
x <- (sample(1:100) - 1) / 99 - 0.5
x <- x - mean(x)
f <- orm(y ~ x)
v <- diag(vcov(f, intercepts='all'))[1:99]
intsx <- coef(f)[1:99]

# Fit model with no predictor
g <- orm.fit(y=y)
w <- diag(vcov(g, intercepts='all'))
ints <- coef(g)
plot(ints, intsx)
hist(ints - intsx, nclass=20)
lm(intsx ~ ints)

plot(w, v); abline(a=0, b=1)
hist(v / w, nclass=20)
median(v / w)  # 2.5 if don't center x, 1.000018 if subtract 0.5 or subtract mean(x)

