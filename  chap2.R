## ----opts, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
setwd('/Users/Yuki/FRA/INLAren/spde-book-files')
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/spde-intro-'
)
library(fields)

## ----label = "nc", echo = FALSE, fig = TRUE, fig.cap = "Counties in North Carolina and their neighborhood structure.", results = 'hide'----
library(spdep)
library(rgdal)

# #Code from spData::nc.sids
# nc.sids <- readOGR(system.file("shapes/sids.shp", package="spData")[1])
# row.names(nc.sids) <- as.character(nc.sids$FIPS)
# rn <- row.names(nc.sids)
# ncCR85_nb <- read.gal(system.file("weights/ncCR85.gal", package="spData")[1],
#                       region.id=rn)
# xx <- coordinates(nc.sids)
# 
# plot(nc.sids, border="grey")
# plot(ncCR85_nb, xx, add = TRUE, pch = 19)


# Mattern correlation -----------------------------------------------------
cMatern = function(h, nu, kappa){
  ifelse(h > 0, besselK(h*kappa, nu)*(h*kappa)^nu/(gamma(nu)*2^(nu-1)), 1)
}


# sample from zero mean multivariate normal -------------------------------
rmvnorm0 = function(n, cov, R = NULL){
  if(is.null(R))
    R = chol(cov)
  
  return(crossprod(R, matrix(rnorm(n*ncol(R)), ncol(R))))
}



# make data ---------------------------------------------------------------
loc = 0:249/25
mdist = as.matrix(dist(loc))


# smoothness and range parameters -----------------------------------------
nu = c(0.5, 1.5, 2.5, 5.5) #smoothness
range = c(1, 4) #range; distance scale of spatial correlation

params = cbind(nu = rep(nu, length(range)), range = rep(range, each = length(nu)))


# noize -------------------------------------------------------------------
set.seed(123)
z = matrix(rnorm(nrow(mdist) * 5), ncol = 5)


## ----samples-------------------------------------------------------------
# Compute the correlated samples
# Scenarios (i.e., different set of parameters)
yy <- lapply(1:nrow(params), function(j) { 
  param <- c(params[j, 1], sqrt(8 * params[j, 1]) / params[j, 2], 
             params[j, 2])
  v <- cMatern(mdist, param[1], param[2])
  
  # fix the diagonal to avoid numerical issues
  diag(v) <- 1 + 1e-9 
  
  # Parameter scenario and computed sample
  return(list(params = param, y = crossprod(chol(v), z)))
})

## ----label = "maternsamples", echo = FALSE, results = 'hide', fig = TRUE, fig.width = 5, fig.cap =  "Five samples from the one-dimensional Matérn correlation function for two different range values (each column of plots) and four different values for the smoothness parameter (each line of plots)."----

par(mfcol = c(4, 2), mar = c(2, 2, 1, 0.1), 
    mgp = c(1.5, 0.5, 0), las = 1)

ry <- range(unlist(lapply(yy, tail, 1)))

# Each scenario
for (i in 1:length(yy)) { 
  plot(loc, yy[[i]]$y[, 1], ylim = ry, type = 'n', 
       xlab = '', ylab = '', 
       main = as.expression(bquote(paste(
         nu==.(yy[[i]]$params[1]), ', ',
         kappa==.(round(yy[[i]]$params[2], 2)), ', ',
         r==.(yy[[i]]$params[3])))))
  
  #Set colours
  cols <- book.color.d(5)
  
  # Each sample
  for (k in 1:5) 
    lines(loc, yy[[i]]$y[, k], col = cols[k])
}



# 2.1.4 make toy data -----------------------------------------------------
## ----rpts----------------------------------------------------------------
n <- 200
set.seed(123) 
pts <- cbind(s1 = sample(1:n / n - 0.5 / n)^2,
             s2 = sample(1:n / n - 0.5 / n)^2)

## ----distpts-------------------------------------------------------------
dmat <- as.matrix(dist(pts))

## ----params--------------------------------------------------------------
beta0 <- 10 #mean of obserbation
sigma2e <- 0.3 # nugget
sigma2u <- 5 # marginal variance of the process
kappa <- 7 # bessel
nu <- 1 # smoothness

## ----covMatm-------------------------------------------------------------
mcor <- cMatern(dmat, nu, kappa) 
mcov <- sigma2e * diag(nrow(mcor)) + sigma2u * mcor # observation error depending on distance, and marginal variance*Matern corr. (i.e., Matern covariance function)

## ----chol1mvnorm---------------------------------------------------------
R <- chol(mcov)
set.seed(234) 
y1 <- beta0 + drop(crossprod(R, rnorm(n))) # mean + random effect for space

## ----grf1, fig = TRUE, fig.height = 5, echo = FALSE, fig.cap = "The simulated toy example data."----
par(mar=c(3,3,1,1), mgp=c(1.5, 0.5, 0), las=1)
plot(pts, asp = 1, xlim = c(0, 1.2), cex = y1 / 10)
q <- quantile(y1, 0:5 / 5)
legend('topright', format(q, dig = 2), 
       pch = 1, pt.cex = q / 10, bty = "n")


## ----datatoy-------------------------------------------------------------
data(SPDEtoy)

## ----label = "strdata"---------------------------------------------------
str(SPDEtoy)

## ----label = "buildmesh5"------------------------------------------------
pl.dom <- cbind(c(0, 1, 1, 0.7, 0), c(0, 0, 0.7, 1, 1))
mesh5 <- inla.mesh.2d(loc.domain = pl.dom, max.e = c(0.092, 0.2))
plot(mesh5)

## ----label = "spde5def"--------------------------------------------------
spde5 <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh5, alpha = 2,
  # P(practic.range < 0.3) = 0.5
  prior.range = c(0.3, 0.5),
  # P(sigma > 1) = 0.01
  prior.sigma = c(10, 0.01)) 


## ----label = "proj2"-----------------------------------------------------
coords <- as.matrix(SPDEtoy[, 1:2])
A5 <- inla.spde.make.A(mesh5, loc = coords)

## ----label = "dima1"-----------------------------------------------------
dim(A5)

## ----label = "a5lines"---------------------------------------------------
table(rowSums(A5 > 0))

## ----label = "rsum"------------------------------------------------------
table(rowSums(A5))

## ----label = "colsA"-----------------------------------------------------
table(colSums(A5) > 0)

## ----label ="eacha1"-----------------------------------------------------
A1 <- inla.spde.make.A(mesh1, loc = coords)

## ----label = "summarya1"-------------------------------------------------
table(as.numeric(A1))

## ----label = "stackdata1b"-----------------------------------------------
stk5 <- inla.stack(
  data = list(resp = SPDEtoy$y),
  A = list(A5, 1), 
  effects = list(i = 1:spde5$n.spde,
                 beta0 = rep(1, nrow(SPDEtoy))),
  tag = 'est')

## ----label = "dimA"------------------------------------------------------
dim(inla.stack.A(stk5))

