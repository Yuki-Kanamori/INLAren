# spatio-temporal model -----------------------------------------
# year for fixed effect and spatio-temporal for random e --------
require(INLA)
require(tidyverse)

setwd("/Users/Yuki/Dropbox/Network2020")
data = read.csv("VASTdata.csv")
head(data, 2)
# mako = data #データ数が多くてmesh作ってプロットするのに時間がかかる
mako = data %>% filter(FISH == "makogarei", between(Y, 2015, 2018)) %>% arrange(Y)
summary(mako)
z = (mako$CPUE > 0) + 0
y = mako$CPUE


# lonlat data
lonlat = as.matrix(cbind(mako$Lon, mako$Lat))

# create a mesh with boundary
# maxedgeを0.03にすると重くて動かない
bound11 = inla.nonconvex.hull(lonlat, convex = 0.05, concave = -0.15)
mesh11 = inla.mesh.2d(boundary = bound11, max.edge = c(0.04, 0.04), cutoff = 0.08/5)
plot(mesh11)
points(lonlat, col = "red", pch = 16, cex = .5)
mesh11$n #529

# spde
cpue_spde = inla.spde2.pcmatern(mesh = mesh11, alpha = 2, prior.range = c(0.01, 0.05), prior.sigma = c(1, 0.01))

# index
# iset = inla.spde.make.index("i", n.spde = cpue_spde$n.spde, n.group = k)

# projection matrix
A_cpue_mako = inla.spde.make.A(mesh11, loc = cbind(mako$Lon, mako$Lat))
dim(A_cpue_mako) # of data times # of vertices in the mesh
table(rowSums(A_cpue_mako > 0))
table(rowSums(A_cpue_mako))
table(colSums(A_cpue_mako) > 0)

# stack for the occurance
stk_z = inla.stack(
  data = list(z = z, y = cbind(z, NA)),
  A = list(A_cpue_mako, 1),
  effects = list(i.z = 1:cpue_spde$n.spde, z.b0 = rep(1, length(z))),
  tag = "est.z"
)

# stack for the density
stk_y = inla.stack(
  data = list(r = y, y = cbind(NA, y)),
  A = list(A_cpue_mako, 1),
  effects = list(i.y = 1:cpue_spde$n.spde, y.b0 = rep(1, length(y))),
  tag = "est.y"
)

# fitting for the occurance
res_z = inla(z ~ 0 + z.b0 + f(i.z, model = cpue_spde), family = "binomial", data = inla.stack.data(stk_z), control.compute = list(dic = TRUE), control.predictor = list(A = inla.stack.A(stk_z), compute = TRUE))

res_y = inla(r ~ 0 + y.b0 + f(i.y, model = cpue_spde), data = inla.stack.data(stk_y), control.compute = list(dic = TRUE), control.predictor = list(A = inla.stack.A(stk_y), compute = TRUE))

# full data stack for the hurdle model
stk_zy = inla.stack(stk_z, stk_y)

# fitting the hurdle model
res_h = inla(y ~ 0 + z.b0 + y.b0 + f(i.z, model = cpue_spde) + f(i.y, copy = "i.z", fixed = FALSE), family = c("binomial", "gaussian"), data = inla.stack.data(stk_zy), control.compute = list(dic = TRUE), control.predictor = list(A = inla.stack.A(stk_zy), compute = TRUE))

