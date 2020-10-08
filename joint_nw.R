require(INLA)
require(tidyverse)
require(openxlsx)


# eDNA ----------------------------------------------------------
setwd("/Users/Yuki/Dropbox/Data")
edf = read.csv("d_2019_merge_df.csv")
colnames(edf)
edf = edf %>% gather(key = sp, value = leads, 8:ncol(edf))
sp_j = data.frame(sp = c("Conger.myriaster", "Pseudopleuronectes.yokohamae", "Kareius.bicoloratus", "Lateolabrax.japonicus", "Konosirus.punctatus"), sp_j = c("maanago", "makogarei", "isigarei", "suzuki", "konosiro"))
edf = merge(edf, sp_j, by = "sp")
unique(edf$sp_j)
unique(edf$Point)


lonlat = read.table("/Users/Yuki/Dropbox/eDNA_INLA/2018/sampling_points.txt", header = T)
head(lonlat, 2)
lonlat = lonlat %>% dplyr::rename(Point = pop) 
lonlat = lonlat %>% mutate(Point = ifelse(lonlat$Point == "Ft", "FT", lonlat$Point))
unique(lonlat$Point)

edf = left_join(edf, lonlat, by = "Point")
summary(edf)
edf[is.na(edf$leads), ] = 0

edf_b = edf %>% filter(Depth == "B") %>% dplyr::rename(lon = lng) %>% mutate(site_id = paste(Month, Point, Depth, Sampling, sep = "_")) #nrow = 1113


# env -----------------------------------------------------------
env = read.csv("/Users/Yuki/Dropbox/Data/Env_data_merged_unique.csv")
summary(env)
unique(env$Variable)
need = c("Temp", "Salinity", "DO", "pH")
env_b = env %>% filter(SorB == "B", Variable %in% need) %>% mutate(site_id = paste(site_id, number, sep = "_"))
unique(env_b$Variable) #nrow = 683

for(i in need){
  assign(paste0(i),
         env_b %>% filter(Variable == i) %>% select(value, site_id)
  )
}

summary(Temp)
summary(Salinity)
summary(DO)
summary(pH)

edf_b = merge(edf_b, Temp, by = "site_id", all = T) 
edf_b = edf_b %>% dplyr::rename(temp = value)
edf_b = merge(edf_b, Salinity, by = "site_id", all = T)
edf_b = edf_b %>% dplyr::rename(salinity = value)
edf_b = merge(edf_b, DO, by = "site_id", all = T)
edf_b = edf_b %>% dplyr::rename(DO = value)
edf_b = merge(edf_b, pH, by = "site_id", all = T)
edf_b = edf_b %>% dplyr::rename(pH = value)



# CPUE ----------------------------------------------------------
cdf = read.csv("/Users/Yuki/Dropbox/eDNA_INLA/added_data2019_2.csv")



# マアナゴ ----------------------------------------------------------
e_ana = edf_b %>% filter(sp_j == "maanago")
summary(e_ana)
edna = (e_ana$leads > 0) + 0
summary(edna)

c_ana = cdf %>% filter(sp == "maanago")
summary(c_ana)
catch = (c_ana$catch > 0) + 0
summary(catch) #1しかない！


# INLA ----------------------------------------------------------
#lonlat
e_loc = as.matrix(cbind(e_ana$lon, e_ana$lat))
c_loc = as.matrix(cbind(c_ana$lon, c_ana$lat))

loc = rbind(e_loc, c_loc)
bound2 = inla.nonconvex.hull(loc, convex = 0.05, concave = -0.15)
mesh2 = inla.mesh.2d(boundary = bound2, max.edge = c(0.04, 0.04), cutoff = 0.08/5)
plot(mesh2)
points(c_loc, col = "red", pch = 16, cex = .5) #pdfにする時はcex = 1にする
points(e_loc, col = "green", pch = 16, cex = .5)
mesh2$n #611

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
                   effects = list(i.e = 1:mesh2$n, list(eb.0 = rep(1, length(edna)), temp = e_ana$temp, salinity = e_ana$salinity, DO = e_ana$DO, pH = e_ana$pH)),
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

# projecting the spatial field ----------------------------------
range_e = apply(mesh2$loc[, c(1, 2)], 2, range)
# range_e = apply(coop, 2, range)
proj_e = inla.mesh.projector(mesh2, xlim = range_e[, 1], ylim = range_e[, 2], dims = c(50, 50))
mean_s_ie = inla.mesh.project(proj_e, res_ana$summary.random$i.e$mean)
sd_s_ie = inla.mesh.project(proj_e, res_ana$summary.random$i.e$sd)

df_e = expand.grid(x = proj_e$x, y = proj_e$y)
df_e$mean_s = as.vector(mean_s_ie)
df_e$sd_s = as.vector(sd_s_ie)

require(viridis)
require(cowplot)

g = ggplot(df_e, aes(x = x, y = y, fill = mean_s_ie))
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
g1 = ggplot(df_e, aes(x = x, y = y, fill = mean_s_ie))
g2 = ggplot(df_e, aes(x = x, y = y, fill = sd_s_ie))
# r = geom_raster()
t = geom_tile()
v = scale_fill_viridis(na.value = "transparent")
c = coord_fixed(ratio = 1)
labs1 = labs(x = "Longitude", y = "Latitude", title = "Mean")
labs2 = labs(x = "Longitude", y = "Latitude", title = "SD")
g1+t+v+c+pol+c_map+labs1+theme_bw()
g2+t+v+c+pol+c_map+labs2+theme_bw()



# マコガレイ ---------------------------------------------------------
e_fish = edf_b %>% filter(sp_j == "makogarei")
summary(e_fish)
edna = (e_fish$leads > 0) + 0
summary(edna)

c_fish = cdf %>% filter(sp == "makogarei")
summary(c_fish)
catch = (c_fish$catch > 0) + 0
summary(catch) #1しかない！
