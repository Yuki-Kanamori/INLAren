df = rpois(100, 1)
summary(df)
z = (df>0)+0
y = ifelse(df>0, df, NA)
data = list(z = z, y = cbind(z, NA))



# ---------------------------------------------------------------
require(tidyverse)

set.seed(100)
s1 = data.frame(runif(100, 0, 1000))
s2 = runif(10, 0, 10) %>% data.frame()
s3 = runif(110, 0, 100) %>% data.frame()

n1 = dim(s1)[1]
n2 = dim(s2)[1]
n3 = dim(s3)[1]

nt1 = rep(NA, n1)
nt2 = rep(NA, n2)
nt3 = rep(NA, n3)

y = as.vector(s1)
yNA = as.vector(c(y, nt2, nt3))

z = as.vector(s2)
zNA = as.vector(c(nt1, z, nt3))

w = as.vector(s3)
wNA = as.vector(c(nt1, nt2, w))

outcome.matrix = matrix(c(yNA, zNA, wNA), ncol = 3)