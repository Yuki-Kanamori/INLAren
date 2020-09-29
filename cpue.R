require(INLA)
require(tidyverse)

setwd("/Users/Yuki/Dropbox/Network2020")
data = read.csv("VASTdata.csv")
head(data, 2)
# mako = data #データ数が多くてmesh作ってプロットするのに時間がかかる
mako = data %>% filter(FISH == "makogarei", Y == 2018)

summary(mako)

dom_tok = cbind(c(139.7, 139.5, 139.7, 140.1, 140.3, 139.9), # x-axis 
                c(35.2,  35.4,  35.8,  35.8,  35.4,  35.2)) # matrix data
plot(dom_tok)
# dom_tok = SpatialPolygons(list(Polygons(list(Polygon(
#   cbind(c(139.7, 139.5, 139.7, 140.1, 140.3, 139.9), # x-axis 
#         c(35.2,  35.4,  35.8,  35.8,  35.4,  35.2)), # y-axis
#   FALSE)), '0')), proj4string = CRS(proj4string(map.sp))) # Formal class SpatialPolygons

cpue_mako_lonlat = as.matrix(cbind(mako$Lon, mako$Lat)) # matrix data
# cpue_mako_lonlat = SpatialPolygons(list(Polygons(list(Polygon(
#   cbind(mako$Lon, # x-axis 
#         mako$Lat), # y-axis
#   FALSE)), '0')), proj4string = CRS(proj4string(map.sp))) # Formal class SpatialPolygons
# 
# kmproj = CRS("+proj=utm +zone=53 ellps=WGS84 +units=km") #with warning
# dom_tok2 = spTransform(dom_tok, kmproj) # Formal class SpatialPolygons
# cpue_mako_lonlat2 = spTransform(cpue_mako_lonlat, kmproj) # Formal class SpatialPolygons
# plot(dom_tok2)

mesh1 = inla.mesh.2d(cpue_mako_lonlat, max.edge = c(0.2, 0.2))
plot(mesh1)
mesh2 = inla.mesh.2d(cpue_mako_lonlat, max.edge = c(1, 1), cutoff = 0.6)
plot(mesh2)
mesh3 = inla.mesh.2d(loc.domain = dom_tok, cpue_mako_lonlat, max.edge = c(0.2, 0.2), cutoff = 0.1)
plot(mesh3)
mesh4 = inla.mesh.2d(loc.domain = dom_tok, cpue_mako_lonlat, max.edge = c(0.2, 0.2), cutoff = 0.1, offset = c(0.5, 0.3))
plot(mesh4)
mesh5 = inla.mesh.2d(loc.domain = dom_tok, cpue_mako_lonlat, max.edge = c(1, 1), cutoff = 0.5, offset = c(0.3, 0.2))
plot(mesh5)
mesh6 = inla.mesh.2d(loc.domain = dom_tok, cpue_mako_lonlat, max.edge = c(0.5, 0.5), cutoff = 0.3, offset = c(0.6, 0.6))
plot(mesh6)
mesh7 = inla.mesh.2d(loc.domain = dom_tok, cpue_mako_lonlat, max.edge = c(0.8, 0.8), cutoff = 0.4, offset = c(0.8, 0.8))
plot(mesh7)
mesh8 = inla.mesh.2d(loc.domain = dom_tok, cpue_mako_lonlat, max.edge = c(0.05, 0.05), cutoff = 0.2, offset = c(0.8, 0.05))
plot(mesh8, asp=1, main = "")

# 暴走する
# mesh9 = inla.mesh.2d(cpue_mako_lonlat, max.edge = c(0.01, 0.01), cutoff = 0.005, offset = c(0.05, 0.05))
# plot(mesh9)

bound9 = inla.nonconvex.hull(cpue_mako_lonlat)
mesh9 = inla.mesh.2d(boundary = bound9, max.edge = c(0.5, 0.5))
plot(mesh9)
points(cpue_mako_lonlat, col = "red", pch = 16, cex = .5)


bound10 = inla.nonconvex.hull(cpue_mako_lonlat, convex = 0.05, concave = -0.15)
mesh10 = inla.mesh.2d(boundary = bound10, max.edge = c(0.08, 0.08), cutoff = 0.02)
plot(mesh10)
points(cpue_mako_lonlat, col = "red", pch = 16, cex = .5)
mesh10$n


# https://haakonbakkagit.github.io/btopic104.htmlの基準を参考に作ってみる
# ポイント1: 境界を1層以上作る
# ポイント2: maxrangeはrange（空間相関の範囲）の1/5以下に（でも空間相関の範囲は解析前に分からないので，サイトではstudy areaの1/3を使っていた）
# ポイント3: いびつな形のメッシュがない（cutoffを使って整える）
# maxedgeは適当（0.08, 0.05, 0.03あたりはどんどんメッシュが細かくなりつつ綺麗な三角形を保っていた．0.02は重くて動かない．バックでfmesherが暴走する）
# cutoffはmaxedgeの1/5に
# 2018年だけでなく全てのデータを使う時には，mesh11の条件を少し帰る必要がある（いびつな三角形が含まれるため）
bound11 = inla.nonconvex.hull(cpue_mako_lonlat, convex = 0.05, concave = -0.15)
mesh11 = inla.mesh.2d(boundary = bound11, max.edge = c(0.03, 0.03), cutoff = 0.08/5)
plot(mesh11)
points(cpue_mako_lonlat, col = "red", pch = 16, cex = .5)
mesh11$n #768

# mesh9 <- inla.mesh.2d(boundary = dom_tok2,
#                       max.edge = c(0.05, 0.05), 
#                       cutoff = 0.2, 
#                       offset = c(0.8, 0.05))
# plot(mesh9)



## PC-priorでrangeとmarginal varianceの範囲がどれくらいか分からない
cpue_spde = inla.spde2.pcmatern(mesh = mesh11, alpha = 2, prior.range = c(0.01, 0.05), prior.sigma = c(1, 0.01))

A_cpue_mako = inla.spde.make.A(mesh11, loc = cpue_mako_lonlat)

dim(A_cpue_mako) # of data times # of vertices in the mesh
table(rowSums(A_cpue_mako > 0))
table(rowSums(A_cpue_mako))
table(colSums(A_cpue_mako) > 0)
# table(as.numeric(A_cpue_mako))

i_index = inla.spde.make.index("i", n.spde = cpue_spde$n.spde, n.group = 1)

# stack_cpue = inla.stack(
#   data = mako$CPUE,
#   A = list(A_cpue_mako, 1),
#   effects = list(list(i = 1:cpue_spde$n.spde), #spatial random effect
#                  list(m = rep(1, nrow(mako)))), #intercept?
#   tag = "mako_cpue"
# )

stack_cpue = inla.stack(
  data = list(cpue = mako$CPUE),
  A = list(A_cpue_mako, 1),
  effects = list(i = i_index, #spatial random effect
                 m = rep(1, nrow(mako))), #intercept?
  tag = "mako_cpue"
)

eq = cpue ~ 0 + m + f(i, model = spde)

test = inla(eq, data = inla.stack.data(stack_cpue), control.predictor = list(A = inla.stack.A(stack_cpue), compute = TRUE))


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
mako$time = mako$Y-2014
mako$w = factor(mako$time)
table(mako$w)
class(mako$w)

# lonlat data
cpue_mako_lonlat = as.matrix(cbind(mako$Lon, mako$Lat))

# time step
# k = length(unique(mako$Y))
k = 4

# make a mesh based on a criteria written by Bakka
# maxedgeを0.03にすると重くて動かない
bound11 = inla.nonconvex.hull(cpue_mako_lonlat, convex = 0.05, concave = -0.15)
mesh11 = inla.mesh.2d(boundary = bound11, max.edge = c(0.04, 0.04), cutoff = 0.08/5)
plot(mesh11)
points(cpue_mako_lonlat, col = "red", pch = 16, cex = .5)
mesh11$n #529

# spde
cpue_spde = inla.spde2.pcmatern(mesh = mesh11, alpha = 2, prior.range = c(0.01, 0.05), prior.sigma = c(1, 0.01))

# index
# n.spde = mesh$n*time = 529*4 = 2116 ??
iset = inla.spde.make.index("i", n.spde = cpue_spde$n.spde, n.group = k)

# projection matrix
A_cpue_mako = inla.spde.make.A(mesh11, loc = cbind(mako$Lon, mako$Lat), group = mako$time)
dim(A_cpue_mako) # of data times # of vertices in the mesh
table(rowSums(A_cpue_mako > 0))
table(rowSums(A_cpue_mako))
table(colSums(A_cpue_mako) > 0)

# stack
sdat = inla.stack(
  data = list(y = mako$CPUE),
  A = list(A_cpue_mako, 1),
  effects = list(iset, w = mako$w),
  tag = "stdata"
)

# AR1
h.spec = list(theta = list(prior = "pccor1", param = c(0, 0.9)))

# foluma
formulae = y ~ 0 + w + f(i, model = cpue_spde, group = i.group, 
                         control.group = list(model = 'ar1', hyper = h.spec))
prec.prior = list(prior = "pc.prec", param = c(1, 0.01))

# fitting
res = inla(formulae, 
           data = inla.stack.data(sdat), 
           control.predictor = list(compute = TRUE, A = inla.stack.A(sdat)),
           control.family = list(hyper = list(theta = prec.prior)),
           control.fixed = list(expand.factor.strategy = "inla"))

res$summary.fixed
res$summary.linear.predictor
