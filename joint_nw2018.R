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

mifish = mifish %>% filter(CPUE == 1)
mifish = mifish %>% mutate(Depth = ifelse(mifish$Depth == "b", "B", "S"), Day = as.numeric(mifish$day)) %>% dplyr::rename(lon = lng)
summary(mifish)
mifish$tag = paste(mifish$month, mifish$Day, mifish$Point, mifish$Depth, sep = "_")


# env -----------------------------------------------------------
env = read.csv("/Users/Yuki/Dropbox/eDNA_INLA/env2018.csv")
env = env %>% mutate(tag = paste(Month, Day, Point, SorB, sep = "_"))
temp = env %>% select(temp, tag)
sal = env %>% select(salinity, tag)
do = env %>% select(DO_mg, tag)
ph = env %>% select(pH, tag)

mifish = merge(mifish, temp, by = "tag", all = T)
mifish = merge(mifish, sal, by = "tag", all = T)
mifish = merge(mifish, do, by = "tag", all = T)
mifish = merge(mifish, ph, by = "tag", all = T)


# CPUE ----------------------------------------------------------
setwd("/Users/Yuki/Dropbox/Network2020")
data = read.csv("VASTdata.csv")
head(data, 2)



# マアナゴ ------------------------------------
e_fish = mifish %>% filter(name_j == "マアナゴ", Depth == "B")
summary(e_fish)
edna = (e_fish$count > 0) + 0

c_fish = data %>% filter(FISH == "maanago", Y == 2018)
# c_mako = data %>% filter(FISH == "makogarei", Y == 2018, M > 2) #eDNAと月を揃えた方が良い？
summary(c_fish)
catch = (c_fish$CATCH > 0) + 0
summary(catch)


# INLA ----------------------------------------------------------
#lonlat
e_loc = as.matrix(cbind(e_fish$lon, e_fish$lat))
c_loc = as.matrix(cbind(c_fish$Lon, c_fish$Lat))

loc = rbind(e_loc, c_loc)
bound2 = inla.nonconvex.hull(loc, convex = 0.05, concave = -0.15)
mesh2 = inla.mesh.2d(boundary = bound2, max.edge = c(0.04, 0.04), cutoff = 0.08/5)
plot(mesh2)
points(c_loc, col = "red", pch = 16, cex = .5) #pdfにする時はcex = 1にする
points(e_loc, col = "green", pch = 16, cex = .5)
mesh2$n #618

# projector matricies
e_A = inla.spde.make.A(mesh2, loc = e_loc)
dim(e_A) # 246, 611
table(rowSums(e_A > 0))
table(rowSums(e_A))
table(colSums(e_A) > 0)
c_A = inla.spde.make.A(mesh2, loc = c_loc)
dim(c_A) # 36, 611
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
dim(Ap) #398, 611

# spde
spde = inla.spde2.pcmatern(mesh = mesh2, alpha = 2, prior.range = c(0.01, 0.05), prior.sigma = c(1, 0.01))

# stack for eDNA
# e_index = inla.spde.make.index("i.e", spde$n.spde)
e_stk = inla.stack(data = list(y = cbind(edna, NA)),
                   A = list(e_A, 1),
                   effects = list(i.e = 1:mesh2$n, list(eb.0 = rep(1, length(edna)), temp = e_fish$temp, salinity = e_fish$salinity, DO = e_fish$DO, pH = e_fish$pH)),
                   tag = "e_dat")
na = as.matrix(cbind(rep(NA, nrow(coop)), rep(NA, nrow(coop))))
ep_stk = inla.stack(data = list(y = cbind(na[, 1], na[, 2])),
                    A = list(Ap, 1),
                    effects = list(i.e = 1:mesh2$n, list(eb.0 = rep(1, nrow(coop)), temp = rep(1, nrow(coop)), salinity = rep(1, nrow(coop)), DO = rep(1, nrow(coop)), pH = rep(1, nrow(coop)))),
                    tag = "ep_dat")
# ep_stk = inla.stack(data = list(y = cbind(na[, 1], na[, 2])),
#                     A = list(Ap, 1),
#                     effects = list(i.e = 1:mesh2$n, list(eb.0 = rep(1, nrow(coop)), temp = e_ana$temp, salinity = e_ana$salinity, DO = e_ana$DO, pH = e_ana$pH)),
#                     tag = "ep_dat")
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
formula = y ~ 0 + eb.0 + cb.0 + f(inla.group(temp), model = "rw1") + f(inla.group(salinity), model = "rw1") + f(inla.group(DO), model = "rw1") + f(inla.group(pH), model = "rw1") + f(i.e, model = spde) + f(x, model = spde) + f(i.c, copy = "i.e", fixed = FALSE)

# fitting
res_ana = inla(formula, data = inla.stack.data(stk), family = c("binomial", "binomial"), control.predictor = list(compute = TRUE, A = inla.stack.A(stk)), control.results = list(return.marginals.random = FALSE, return.marginals.predictor = FALSE))
summary(res_ana)

# plot the fitted values on a map -------------------------------
index_ep = inla.stack.index(stk, tag = "ep_dat")$data
index_cp = inla.stack.index(stk, tag = "cp_dat")$data

pred_mean_e = res_ana$summary.fitted.values[index_ep, "mean"]
pred_mean_c = res_ana$summary.fitted.values[index_cp, "mean"]
pred_ll_e = res_ana$summary.fitted.values[index_ep, "0.025quant"]
pred_ul_e = res_ana$summary.fitted.values[index_ep, "0.975quant"]
pred_ll_c = res_ana$summary.fitted.values[index_cp, "0.025quant"]
pred_ul_c = res_ana$summary.fitted.values[index_cp, "0.975quant"]

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
s = scale_fill_gradient(name = "encounter prob. (logit)", low = "blue", high = "orange")
g1+t+f+c+s+theme_bw()
g2+t+f+c+s+theme_bw()


# with map
world_map <- map_data("world")
jap <- subset(world_map, world_map$region == "Japan")
jap_cog <- jap[jap$lat > 35 & jap$lat < 38 & jap$long > 139 & jap$long < 141, ]
pol = geom_polygon(data = jap_cog, aes(x=long, y=lat, group=group), colour="gray 50", fill="gray 50")
c_map = coord_map(xlim = c(139.5, 140.3), ylim = c(35, 35.75))

dpm = rbind(dpm_e, dpm_c)
m_dpm = dpm %>% filter(str_detect(variable, "mean"))
unique(m_dpm$variable)

g = ggplot(data = m_dpm, aes(east, north, fill = value))
t = geom_tile()
f = facet_wrap(~ variable)
c = coord_fixed(ratio = 1)
s = scale_fill_gradient(name = "encounter prob. (logit)", low = "blue", high = "orange")
g+t+f+c+s+pol+c_map+theme_bw()+labs(title = "maanago")


# projecting the spatial field ----------------------------------
range_e = apply(mesh2$loc[, c(1, 2)], 2, range)
# range_e = apply(coop, 2, range)
proj_e = inla.mesh.projector(mesh2, xlim = range_e[, 1], ylim = range_e[, 2], dims = c(50, 50))
mean_s_ie = inla.mesh.project(proj_e, res_ana$summary.random$i.e$mean)
sd_s_ie = inla.mesh.project(proj_e, res_ana$summary.random$i.e$sd)

df_ie = expand.grid(x = proj_e$x, y = proj_e$y)
df_ie$mean_s = as.vector(mean_s_ie)
df_ie$sd_s = as.vector(sd_s_ie)

require(viridis)
require(cowplot)
require(gridExtra)

g = ggplot(df_ie, aes(x = x, y = y, fill = mean_s_ie))
r = geom_raster()
v = scale_fill_viridis(na.value = "transparent")
c = coord_fixed(ratio = 1)
g+r+v+c+theme_bw()

# with map
world_map <- map_data("world")
jap <- subset(world_map, world_map$region == "Japan")
jap_cog <- jap[jap$lat > 35 & jap$lat < 38 & jap$long > 139 & jap$long < 141, ]
pol = geom_polygon(data = jap_cog, aes(x=long, y=lat, group=group), colour="gray 50", fill="gray 50")
c_map = coord_map(xlim = c(139.5, 140.3), ylim = c(35, 35.75))

g1 = ggplot(df_ie, aes(x = x, y = y, fill = mean_s_ie))
g2 = ggplot(df_ie, aes(x = x, y = y, fill = sd_s_ie))
# r = geom_raster()
t = geom_tile()
v = scale_fill_viridis(na.value = "transparent")
c = coord_fixed(ratio = 1)
labs1 = labs(x = "Longitude", y = "Latitude", title = "")
labs2 = labs(x = "Longitude", y = "Latitude", title = "")
m = g1+t+v+c+pol+c_map+labs1+theme_bw()
sd = g2+t+v+c+pol+c_map+labs2+theme_bw()
grid.arrange(m, sd, ncol = 2)







