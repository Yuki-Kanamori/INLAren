setwd("/Users/Yuki/Dropbox/eDNA_INLA")

require(INLA)
require(tidyverse)
require(openxlsx)

# 2019 data -----------------------------------------------------
# not have CPUE data --------------------------------------------
t = read.table("2019/COUNT.mifish1.nrOTU.id97.unoise.txt", header = T)
time1 = read.table ("times1.txt")


# 2018 ----------------------------------------------------------
mifish = read.table("2018/countOTUv2.fish.txt", header = T) %>% gather(key = OTU, value = count, -(1:3)) %>% mutate(year = 2018, month = as.numeric(str_sub(Date, 3, 4)), day = str_sub(Date, -2, -1))
head(mifish, 2)

sp_name = openxlsx::read.xlsx("2018/OTU_Species.xlsx")
head(sp_name, 2)
mifish = left_join(mifish, sp_name, by = "OTU")

lonlat = read.table("2018/sampling_points.txt", header = T)
head(lonlat, 2)
lonlat = lonlat %>% dplyr::rename(Point = pop)
mifish = left_join(mifish, lonlat, by = "Point")

df = mifish %>% filter(CPUE == 1) 
df = df %>% mutate(pa = ifelse(df$count > 0, 1, 0))
unique(df$name_j)


# env data ------------------------------------------------------
# site名がenvとednaで違うため，確認後修正が必要
env = read.csv("2018/env_data.csv", fileEncoding = "CP932")
summary(env)
env2018 = env %>% filter(year == 2018)
summary(env2018)
unique(env2018$site)


# INLA ----------------------------------------------------------
# state-space model -----------------------------------------------
# df_kono = df %>% filter(name_j == "コノシロ")
# summary(df_kono)
# 
# rep = df_kono %>% group_by(Point, Date, Depth) %>% summarize(repcount = n())
# summary(rep)


df_mako = df %>% filter(name_j == "マコガレイ")
# df_mako = df %>% filter(name_j == "コノシロ")
summary(df_mako)

rep = df_mako %>% group_by(Point, Date, Depth) %>% summarize(repcount = n())
summary(rep)
df_mako_b = df_mako %>% filter(Depth == "b") %>% mutate(rep = 1)


# data("SPDEtoy")
# coords = as.matrix(SPDEtoy[, 1:2])
# p5 = coords[1:5, ]

dom_tok = cbind(c(139.7, 139.5, 139.7, 140.1, 140.3, 139.9), # x-axis 
            c(35.2,  35.4,  35.8,  35.8,  35.4,  35.2))
mako_lonlat = cbind(df_mako_b$lng, df_mako_b$lat)

mesh1 = inla.mesh.2d(mako_lonlat, max.edge = c(0.2, 0.2))
plot(mesh1)
mesh2 = inla.mesh.2d(mako_lonlat, max.edge = c(0.2, 0.2), cutoff = 0.1)
plot(mesh2)
mesh3 = inla.mesh.2d(loc.domain = dom_tok, mako_lonlat, max.edge = c(0.2, 0.2), cutoff = 0.1)
plot(mesh3)
mesh4 = inla.mesh.2d(loc.domain = dom_tok, mako_lonlat, max.edge = c(0.2, 0.2), cutoff = 0.1, offset = c(0.5, 0.3))
plot(mesh4)
mesh5 = inla.mesh.2d(loc.domain = dom_tok, mako_lonlat, max.edge = c(1, 1), cutoff = 0.5, offset = c(0.3, 0.2))
plot(mesh5)
# mesh6 = inla.mesh.2d(loc.domain = dom_tok, mako_lonlat, max.edge = c(0.9, 0.9), cutoff = 0.5, offset = c(0.3, 0.3))
# plot(mesh6)


## PC-priorでrangeとmarginal varianceの範囲がどれくらいか分からない
spde = inla.spde2.pcmatern(mesh = mesh5, alpha = 2, prior.range = c(0.01, 0.05), prior.sigma = c(1, 0.01))

A_mako = inla.spde.make.A(mesh5, loc = mako_lonlat)

dim(A_mako) #199, 35; # of data times # of vertices in the mesh
table(rowSums(A_mako > 0))
table(rowSums(A_mako))
table(colSums(A_mako) > 0)
table(as.numeric(A_mako))

## make data for inla.stack()
n = nrow(df_mako_b)
y = matrix(NA, ncol = 2, nrow = 2*n)
y[1:n, 1] = df_mako_b$pa
y[-c(1:n), 2] = 0

# index for z
m1 = c(1:n, 1:n)
# index for env
m2 = c(rep(NA, n), 1:n)
# weight for env
w1 = c(rep(NA, n), rep(1, n))
# weight for spatial random effect
# w2 = 

## inla.stack()
# error; length(A)=2 should be equal to length(effects)=1 
# -> effectにはspatial random effectと固定効果の両方が必要？
# error; Row count mismatch for A: 199,1
# -> 観察データとfaked zeroデータの両方に対してstackが必要？(黑INLA 3.3 copying part of or the entire linear predictorを参照)

# stack = inla.stack(
#   data = y,
#   A = list(A_mako, 1),
#   effects = list(i = 1:spde$n.spde,
#                  ratio = 1),
#   tag = "makogarei"
# )

# observation data
# y_i = m_i*z_i + s_i
# l[[k]] <- rep(l[[k]], n.A) でエラー: 
# 置き換えられる個数よりも多くの要素が与えられました
# 
# parse.input.list(effects[[k]], input.ncol(A[[k]]), paste("Effect block ",  でエラー: 
# Effect block 2:
# Mismatching row sizes: 1440, n.A=35
# 
# (function (data, A, effects, tag = "", compress = TRUE, remove.unused = TRUE)  でエラー: 
# Row count mismatch for A: 35,199
stack_obs = inla.stack(
  data = df_mako_b$pa,
  A = list(A_mako, 1),
  effects = list(list(i = 1:spde$n.spde), #spatial random effect
                 list(m = rep(1, nrow(df_mako_b)))), #weight for latent variable, z_i
  tag = "mako_obs"
)

# faked zero data (system model)
# z_i = alpha + beta*enc_i + s_i

# l[[k]] <- rep(l[[k]], n.A) でエラー: 
# 置き換えられる個数よりも多くの要素が与えられました

# parse.input.list(effects[[k]], input.ncol(A[[k]]), paste("Effect block ",  でエラー: 
# Effect block 2:
# Mismatching row sizes: 394, n.A=284
stack_system = inla.stack(
  data = rep(0, n),
  A = list(A_mako, 1),
  effects = list(i = 1:spde$n.spde,
                 alpha = rep(1, n)),
  tag = "mako_fake"
)

stack = inla.stack(stack_obs, stack_system)

# fomula = y ~ 
  


data("SPDEtoy")

pl.dom <- cbind(c(0, 1, 1, 0.7, 0), c(0, 0, 0.7, 1, 1))
mesh5 <- inla.mesh.2d(loc.domain = pl.dom, max.e = c(0.092, 0.2))
plot(mesh5)

spde5 <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter 
  mesh = mesh5, alpha = 2,
  # P(practic.range < 0.3) = 0.5 
  prior.range = c(0.3, 0.5),
  # P(sigma > 1) = 0.01 
  prior.sigma = c(10, 0.01))

coords <- as.matrix(SPDEtoy[, 1:2])
A5 <- inla.spde.make.A(mesh5, loc = coords)
dim(A5)
table(rowSums(A5 > 0))
table(rowSums(A5))
table(colSums(A5) > 0)
table(as.numeric(A5))
A1 <- inla.spde.make.A(mesh1, loc = coords)
dim(A1)
table(rowSums(A1 > 0))
table(rowSums(A1))
table(colSums(A1) > 0)
table(as.numeric(A1))

stk5 <- inla.stack(
  data = list(resp = SPDEtoy$y),
  A = list(A5, 1),
  effects = list(i = 1:spde5$n.spde,
                 beta0 = rep(1, nrow(SPDEtoy))), 
  tag = 'est')




  
  
  
  
  
# n = 200
# loc = matrix(runif(n*2),n,2)
# mesh = inla.mesh.create.helper(points.domain=loc,
#                                max.edge=c(0.05, 0.2))
# proj.obs = inla.mesh.projector(mesh, loc=loc)
# proj.pred = inla.mesh.projector(mesh, loc=mesh$loc)
# spde = inla.spde2.matern(mesh,
#                          B.tau=cbind(log(1), 1, 0),
#                          B.kappa=matrix(c(log(sqrt(8)/0.2), 0, 1), 1, 3))
# 
# covar = rnorm(n)
# field = inla.qsample(n=1, Q=inla.spde2.precision(spde, theta=c(0,0)))[,1]
# y = 2*covar + inla.mesh.project(proj.obs, field)
# 
# A.obs = inla.spde.make.A(mesh, loc=loc)
# A.pred = inla.spde.make.A(mesh, loc=proj.pred$loc)
# dim(A.obs)
# stack.obs =
#   inla.stack(data=list(y=y),
#              A=list(A.obs, 1),
#              effects=list(c(list(intercept=rep(1, mesh$n)),
#                             inla.spde.make.index("spatial", spde$n.spde)),
#                           covar=covar),
#              tag="obs")