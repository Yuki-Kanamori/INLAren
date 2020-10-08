require(INLA)
require(tidyverse)
require(openxlsx)

# 2018 ----------------------------------------------------------
# eDNA ----------------------------------------------------------
setwd("/Users/Yuki/Dropbox/eDNA_INLA")
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

e_mako = df %>% filter(name_j == "マコガレイ", Depth == "b")
summary(e_mako)
edna = (e_mako$count > 0) + 0

# CPUE ----------------------------------------------------------
setwd("/Users/Yuki/Dropbox/Network2020")
data = read.csv("VASTdata.csv")
head(data, 2)
c_mako = data %>% filter(FISH == "makogarei", Y == 2018)
# c_mako = data %>% filter(FISH == "makogarei", Y == 2018, M > 2) #eDNAと月を揃えた方が良い？
summary(c_mako)
catch = (c_mako$CATCH > 0) + 0


# INLA ----------------------------------------------------------
#lonlat
e_loc = as.matrix(cbind(e_mako$lng, e_mako$lat))
c_loc = as.matrix(cbind(c_mako$Lon, c_mako$Lat))

# create a mesh with boundary
# maxedgeを0.03にすると重くて動かない
bound1 = inla.nonconvex.hull(c_loc, convex = 0.05, concave = -0.15)
mesh1 = inla.mesh.2d(boundary = bound1, max.edge = c(0.04, 0.04), cutoff = 0.08/5)
plot(mesh1)
points(c_loc, col = "red", pch = 16, cex = .5)
points(e_loc, col = "green", pch = 16, cex = .5)
mesh1$n #461

# eDNAがはみ出るので，両方のデータの緯度経度からメッシュを作成
loc = rbind(e_loc, c_loc)
bound2 = inla.nonconvex.hull(loc, convex = 0.05, concave = -0.15)
mesh2 = inla.mesh.2d(boundary = bound2, max.edge = c(0.04, 0.04), cutoff = 0.08/5)
plot(mesh2)
points(c_loc, col = "red", pch = 16, cex = .5)
points(e_loc, col = "green", pch = 16, cex = .5)
mesh2$n #618

# projector matricies
e_A = inla.spde.make.A(mesh2, loc = e_loc)
dim(e_A) # 199, 618
table(rowSums(e_A > 0))
table(rowSums(e_A))
table(colSums(e_A) > 0)
c_A = inla.spde.make.A(mesh2, loc = c_loc)
dim(c_A) # 2321, 618
table(rowSums(c_A > 0))
table(rowSums(c_A))
table(colSums(c_A) > 0)


# for prediction ------------------------------------------------
setwd('/Users/Yuki/FRA/INLAren/spde-book-files')

## ----opts, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/barrier-'
)
library(scales)
library(rgeos)
## High resolution maps when using map()
library(mapdata) 
## Map features, map2SpatialPolygons()
library(maptools)


# Tokyo Bay ---------------------------------------------------------------
# Select region 
map <- map("world", "Japan", fill = TRUE,
           col = "transparent", plot = TRUE)
IDs <- sapply(strsplit(map$names, ":"), function(x) x[1])
map.sp <- map2SpatialPolygons(
  map, IDs = IDs,
  proj4string = CRS("+proj=longlat +datum=WGS84")) #緯度経度データ
summary(map.sp)

# make a polygon ------------------------------------------------------------
pl.sel <- SpatialPolygons(list(Polygons(list(Polygon(
  cbind(c(139.7, 139.5, 139.7, 140.1, 140.3, 139.9), # x-axis 
        c(35.2,  35.4,  35.8,  35.8,  35.4,  35.2)), # y-axis
  FALSE)), '0')), proj4string = CRS(proj4string(map.sp))) #緯度経度データ
summary(pl.sel)
poly.water <- gDifference(pl.sel, map.sp)
plot(pl.sel)
plot(map.sp)
plot(poly.water)

tok_bor = poly.water@polygons[[1]]@Polygons[[1]]@coords

bb_tok = poly.water@bbox
x = seq(bb_tok[1, "min"] - 1, bb_tok[1, "max"] + 1, length.out = 150)
y = seq(bb_tok[2, "min"] - 1, bb_tok[2, "max"] + 1, length.out = 150)
coop = as.matrix(expand.grid(x, y))
ind = point.in.polygon(coop[, 1], coop[, 2],
                       tok_bor[, 1], tok_bor[, 2])
coop = coop[which(ind == 1), ]
plot(coop, asp = 1)

Ap = inla.spde.make.A(mesh = mesh2, loc = coop)
dim(Ap) #398, 618


# spde
spde = inla.spde2.pcmatern(mesh = mesh2, alpha = 2, prior.range = c(0.01, 0.05), prior.sigma = c(1, 0.01))

# stack for eDNA
# e_index = inla.spde.make.index("i.e", spde$n.spde)
e_stk = inla.stack(data = list(y = cbind(edna, NA)),
                   A = list(e_A, 1),
                   effects = list(i.e = 1:mesh2$n, eb.0 = rep(1, length(edna))),
                   tag = "e_dat")
na = as.matrix(cbind(rep(NA, nrow(coop)), rep(NA, nrow(coop))))
ep_stk = inla.stack(data = list(y = cbind(na[, 1], na[, 2])),
                   A = list(Ap, 1),
                   effects = list(i.e = 1:mesh2$n, eb.0 = rep(1, nrow(coop))),
                   tag = "ep_dat")
stk_edna = inla.stack(e_stk, ep_stk)

c_stk = inla.stack(data = list(y = cbind(NA, catch)),
                   A = list(c_A, 1),
                   effects = list(list(i.c = 1:mesh2$n, x = 1:mesh2$n), cb.0 = rep(1, length(catch))),
                   tag = "c_dat")
cp_stk = inla.stack(data = list(y = cbind(na[, 1], na[, 2])),
                   A = list(Ap, 1),
                   effects = list(list(i.c = 1:mesh2$n, x = 1:mesh2$n), cb.0 = rep(1, nrow(coop))),
                   tag = "cp_dat")
stk_catch = inla.stack(c_stk, cp_stk)
# stk = inla.stack(e_stk, c_stk)
stk = inla.stack(stk_edna, stk_catch)

# formula
formula = y ~ 0 + eb.0 + cb.0 + f(i.e, model = spde) + f(x, model = spde) + f(i.c, copy = "i.e", fixed = FALSE)

# fitting the joint model
# res_joint = inla(formula, data = inla.stack.data(stk), family = c("binomial", "binomial"), control.predictor = list(compute = TRUE, A = inla.stack.A(stk)))

# not estimate the posterior marginal distribution
res_joint = inla(formula, data = inla.stack.data(stk), family = c("binomial", "binomial"), control.predictor = list(compute = TRUE, A = inla.stack.A(stk)), control.results = list(return.marginals.random = FALSE, return.marginals.predictor = FALSE))



# plot the fitted values on a map -------------------------------
index_ep = inla.stack.index(stk, tag = "ep_dat")$data
index_cp = inla.stack.index(stk, tag = "cp_dat")$data

pred_mean_e = res_joint$summary.fitted.values[index_ep, "mean"]
pred_mean_c = res_joint$summary.fitted.values[index_cp, "mean"]
pred_ll_e = res_joint$summary.fitted.values[index_ep, "0.025quant"]
pred_ul_e = res_joint$summary.fitted.values[index_ep, "0.975quant"]
pred_ll_c = res_joint$summary.fitted.values[index_cp, "0.025quant"]
pred_ul_c = res_joint$summary.fitted.values[index_cp, "0.975quant"]

# test1 = data.frame(east = coop[, 1], north = coop[, 2],
#                    value = pred_mean_e, variable = "pred_mean_eDNA")
# test2 = data.frame(east = coop[, 1], north = coop[, 2],
#                    value = pred_ll_e, variable = "pred_ll_eDNA")
# test3 = data.frame(east = coop[, 1], north = coop[, 2],
#                    value = pred_ul_e, variable = "pred_ul_eDNA")
# colnames(test1)
# colnames(test2)
# colnames(test3)
# dpm_e = rbind(test1, test2, test3)

dpm_e = rbind(data.frame(east = coop[, 1], north = coop[, 2],
                         value = pred_mean_e, variable = "pred_mean_eDNA"),
              data.frame(east = coop[, 1], north = coop[, 2],
                         value = pred_ll_e, variable = "pred_ll_eDNA"),
              data.frame(east = coop[, 1], north = coop[, 2],
                         value = pred_ul_e, variable = "pred_ul_eDNA"))
dpm_c = rbind(data.frame(east = coop[, 1], north = coop[, 2],
                         value = pred_mean_c, variable = "pred_mean_catch"),
              data.frame(east = coop[, 1], north = coop[, 2],
                         value = pred_ll_c, variable = "pred_ll_catch"),
              data.frame(east = coop[, 1], north = coop[, 2],
                         value = pred_ul_c, variable = "pred_ul_catch"))

dpm_e$variable = as.factor(dpm_e$variable)
dpm_c$variable = as.factor(dpm_c$variable)

g1 = ggplot(data = dpm_e, aes(east, north, fill = value))
g2 = ggplot(data = dpm_c, aes(east, north, fill = value))
t = geom_tile()
f = facet_wrap(~ variable)
c = coord_fixed(ratio = 1)
s = scale_fill_gradient(name = "encounter probability", low = "blue", high = "orange")
g1+t+f+c+s+theme_bw()
g2+t+f+c+s+theme_bw()

# outputs -------------------------------------------------------
summary(res_joint)
# fixed effects
res_joint$summary.fix # intercepts

# ???
post.se <- inla.tmarginal(function(x) sqrt(1/x), res_joint$marginals.hy[[1]])
inla.emarginal(function(x) x, post.se)
inla.hpdmarginal(0.95, post.se)
data.inla.field <- inla.spde2.result(res_joint, "i.e", spde, do.transf = TRUE)
inla.emarginal(function(x) x, data.inla.field$marginals.kappa[[1]])
inla.emarginal(function(x) x, data.inla.field$marginals.variance.nominal[[1]])


# prediction ----------------------------------------------------
summary(loc)
# dom_tok = cbind(c(139.7, 139.5, 139.7, 140.1, 140.3, 139.9), # x-axis 
#                 c(35.2,  35.4,  35.8,  35.8,  35.4,  35.2)) # matrix data
coords.grid <- as.matrix(expand.grid(long = seq(139.5, 140.5, len = 100), lat = seq(35, 36, len = 100)))
head(coords.grid)
# data.inla.projector <- inla.mesh.projector(mesh2, loc = loc)
data.inla.projector <- inla.mesh.projector(mesh2, loc = coords.grid)
book.plot.field(field = ((res_joint$summary.random$i.e$mean + res_joint$summary.random$x$mean + res_joint$summary.random$i.c$mean) + res_joint$summary.fixed$mean[1]), projector = coords.grid)


index = inla.stack.index()





#?
newdata <- data.frame(loc, mean = inla.mesh.project(data.inla.projector, res_joint$summary.random$i.e$mean + res_joint$summary.random$x$mean + res_joint$summary.random$i.c$mean) + res_joint$summary.fixed$mean[1])
str(newdata$mean)

# ggplot(newdata, aes(y = X2, x = X1)) + geom_tile(aes(fill = mean))
# ggplot(loc, aes(y = V1, x = V2)) + geom_point(aes(color = y), size = 2)

p_test = inla.mesh.projector(mesh2, xlim = 139.5:140.5, ylim = 35:36)
pred_mean = inla.mesh.project(p_test, res_joint$summary.random$i.e$mean)
plot(pred_mean)



# from Dina -----------------------------------------------------
ggplot() +  geom_point(data = data.frame(loc), aes(x = X1, y = X2, color = res_joint$summary.random$i.c$mean), size=1.5) + theme(text=element_text(size=5),axis.text.x=element_text(size=5),axis.text.y=element_text(size=5), panel.background=element_rect(fill="light grey"),panel.grid.minor=element_blank(),panel.grid.major=element_blank())+ scale_color_gradient2("",low="blue", high="red")
# ggplot() +  geom_point(data = data0, aes(x = Lon, y = Lat, color = res_joint$summary.random$i.c$mean), size=1.5) + coord_cartesian(xlim = c(-10, 10), ylim = c(48, 62))  + theme(text=element_text(size=5),axis.text.x=element_text(size=5),axis.text.y=element_text(size=5), panel.background=element_rect(fill="light grey"),panel.grid.minor=element_blank(),panel.grid.major=element_blank())+ scale_color_gradient2("",low="blue", high="red")

test = data.frame(as.matrix(e_A))
test = test %>% mutate(data = 1:nrow(test))
ncol(test)
test2 = test %>% gather(key = projector, value = sparse, 1:(ncol(test)-1))
test3 = test2 %>% filter(sparse > 0)



# 何か描けた ---------------------------------------------------------
source("/Users/Yuki/FRA/INLAren/spde-book-files/R/spde-book-functions.R")
prj = inla.mesh.projector(mesh2, xlim = c(min(loc[,1]), max(loc[,1])), ylim = c(min(loc[,2]), max(loc[,2])))
book.plot.field(field = res_joint$summary.random$i.e$mean, projector = prj)
book.plot.field(field = res_joint$summary.random$i.e$sd, projector = prj)
book.plot.field(field = res_joint$summary.random$i.c$mean, projector = prj)
book.plot.field(field = res_joint$summary.random$i.c$sd, projector = prj)
book.plot.field(field = res_joint$summary.random$x$mean, projector = prj)
book.plot.field(field = res_joint$summary.random$x$sd, projector = prj)



# test for prediction (from inlabru) ----------------------------
library(INLA) 
library(inlabru) 
library(RColorBrewer) 
data(gorillas)
gnests <- gorillas$nests
mesh <- gorillas$mesh
boundary <- gorillas$boundary
gcounts <- gorillas$plotsample$counts 
plots <- gorillas$plotsample$plots 
plotnests <- gorillas$plotsample$nests

fig.plots <- ggplot() +
  gg(boundary) +
  gg(plots) +
  gg(gnests, size=0.07,color = "blue") +
  gg(plotnests,size=0.07,color = "red") + geom_label(aes(label=gcounts$count, x=gcounts$x+0.35, y=gcounts$y-0.3), fill = "white", size = 1.5, label.padding = unit(0.1, "lines")) + coord_fixed() + xlab("x") + ylab("y")
fig.plots

mdl <- count ~ spat(map = coordinates, model = inla.spde2.matern(mesh, constr = FALSE)) + Intercept
fit = bru(mdl, family = "poisson", data = gcounts, options = list(E = gcounts$exposure))

pxl <- pixels(mesh, mask = boundary, nx = 500, ny = 400)
dens <- predict(fit, pxl, ~ exp(spat+Intercept))

pxl2 <- pixels(mesh2, mask = NULL, nx = 50, ny = 50)
dens2 <- predict(fit, pxl, ~ exp(spat+Intercept))

x = 150; y = 150
lattice <- INLA::inla.mesh.lattice(x = x, y = y)


pxl <- pixels(, nx = 50, ny = 50)
