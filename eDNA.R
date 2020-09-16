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

df = mifish %>% filter(CPUE == 1) %>% mutate(pa = ifelse(df$count > 0, 1, 0))
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

## PC-priorでrangeとmarginal varianceの範囲がどれくらいか分からない
spde = inla.spde2.pcmatern(mesh = mesh4, alpha = 2, prior.range = c(0.01, 0.05), prior.sigma = c(1, 0.01))

A_mako = inla.spde.make.A(mesh4, loc = mako_lonlat)

dim(A_mako) #199, 284
table(rowSums(A_mako > 0))
table(rowSums(A_mako))
table(colSums(A_mako) > 0)

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
stack_obs = inla.stack(
  data = df_mako_b$pa,
  A = list(A_mako, 1),
  effects = list(i = 1:spde$n.spde,
                 m = 1:n),
  tag = "mako_obs"
)

# faked zero data (system model)
# z_i = alpha + beta*enc_i + s_i
# l[[k]] <- rep(l[[k]], n.A) でエラー: 
# 置き換えられる個数よりも多くの要素が与えられました
stack_system = inla.stack(
  data = rep(0, n),
  A = list(A_mako, 1),
  effects = list(i = 1:spde$n.spde,
                 alpha = rep(1, n)),
  tag = "mako_fake"
)

stack = inla.stack()