setwd("/Users/Yuki/Dropbox/Network2020")
data = read.csv("VASTdata.csv")
head(data, 2)
mako = data %>% filter(FISH == "makogarei")

dom_tok = cbind(c(139.7, 139.5, 139.7, 140.1, 140.3, 139.9), # x-axis 
                c(35.2,  35.4,  35.8,  35.8,  35.4,  35.2))
cpue_mako_lonlat = cbind(mako$Lon, mako$Lat)
mesh1 = inla.mesh.2d(cpue_mako_lonlat, max.edge = c(0.2, 0.2))
plot(mesh1)
mesh2 = inla.mesh.2d(cpue_mako_lonlat, max.edge = c(0.2, 0.2), cutoff = 0.1)
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

## PC-priorでrangeとmarginal varianceの範囲がどれくらいか分からない
cpue_spde = inla.spde2.pcmatern(mesh = mesh4, alpha = 2, prior.range = c(0.01, 0.05), prior.sigma = c(1, 0.01))

A_cpue_mako = inla.spde.make.A(mesh4, loc = cpue_mako_lonlat)

dim(A_cpue_mako) #199, 35; # of data times # of vertices in the mesh
table(rowSums(A_cpue_mako > 0))
table(rowSums(A_cpue_mako))
table(colSums(A_cpue_mako) > 0)
# table(as.numeric(A_cpue_mako))

stack_cpue = inla.stack(
  data = mako$CPUE,
  A = list(A_cpue_mako, 1),
  effects = list(list(i = 1:cpue_spde$n.spde), #spatial random effect
                 list(m = rep(1, nrow(mako)))), #intercept?
  tag = "mako_cpue"
)
