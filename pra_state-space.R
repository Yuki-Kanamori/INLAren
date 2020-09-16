
# Red INLA 8.5 Space-state model --------------------------------
require(KFAS)
data("alcohol")

# make two column matrix; the first is observation (y) and the second is latent (x)
n = nrow(alcohol) -1 
y = matrix(NA, ncol = 2, nrow = n+(n-1))
y[1:n, 1] = alcohol[1:n, 1]
y[-c(1:n), 2] = 0

# offset
oset = c(alcohol[1:n, 5], rep(NA, n-1))

# x_t
i <- c(1:n, 2:n)
# x_(t-1) 2:n
j <- c(rep(NA, n), 1:(n - 1))
# Weight to have -1 * x_(t-1)
w1 <- c(rep(NA, n), rep(-1, n - 1)) 
# x_(t-1), 2:n
l <- c(rep(NA, n), 2:n)
# Weight to have * omega_(t-1)
w2 <- c(rep(NA, n), rep(-1, n - 1))

prec.prior <- list(prec = list(param = c(0.001, 0.001))) 
alc.inla <- inla(Y ~ 0 + offset(log(oset)) + 
                   f(i, model = "iid", hyper = list(prec = list(initial = -10, fixed = TRUE))) +
                   f(j, w1, copy = "i") + 
                   f(l, w2, model = "iid"), 
                   data = list(Y = Y, oset = oset), 
                   family = c("poisson", "gaussian"), 
                   control.family = list(list(), list(hyper = list(prec = list(initial = 10, fixed = TRUE)))), 
                   control.predictor = list(compute = TRUE))

rep(1:10, 2)


# binomial ------------------------------------------------------
n = 100
z = runif(n) #ノイズ
eta = 1 + 0.1*z
N = 20

## logit
p = exp(eta)/(1+exp(eta))
y = rbinom(n, size = N, prob = p) #n = データ数, size = 各乱数におけるベルヌーイ試行の回数, prob = 各ベルヌーイ試行における成功確率
r = inla(y ~ 1 + z, data = data.frame(y, z), 
         family = "binomial", 
         Ntrials = rep(N,n), 
         control.family = list(link = "logit"), 
         control.predictor = list(compute = TRUE))
