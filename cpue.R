require(INLA)
require(tidyverse)

setwd("/Users/Yuki/Dropbox/Network2020")
data = read.csv("VASTdata.csv")
head(data, 2)
mako = data %>% filter(FISH == "makogarei", Y == 2018)

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
mesh8 = inla.mesh.2d(loc.domain = dom_tok, cpue_mako_lonlat, max.edge = c(0.05, 0.05), cutoff = 0.2, offset = c(0.8, 0.05))
plot(mesh8, asp=1, main = "")


mesh9 = inla.mesh.2d(cpue_mako_lonlat, max.edge = c(0.01, 0.01), cutoff = 0.005, offset = c(0.05, 0.05))
plot(mesh9)



points(cpue_mako_lonlat, col = "red", pch = 16, cex = .5)
mesh8$n

# mesh9 <- inla.mesh.2d(boundary = dom_tok2,
#                       max.edge = c(0.05, 0.05), 
#                       cutoff = 0.2, 
#                       offset = c(0.8, 0.05))
# plot(mesh9)



## PC-priorでrangeとmarginal varianceの範囲がどれくらいか分からない
cpue_spde = inla.spde2.pcmatern(mesh = mesh8, alpha = 2, prior.range = c(0.01, 0.05), prior.sigma = c(1, 0.01))

A_cpue_mako = inla.spde.make.A(mesh8, loc = cpue_mako_lonlat)

dim(A_cpue_mako) #199, 35; # of data times # of vertices in the mesh
table(rowSums(A_cpue_mako > 0))
table(rowSums(A_cpue_mako))
table(colSums(A_cpue_mako) > 0)
# table(as.numeric(A_cpue_mako))

i_index = inla.spde.make.index("i", n.spde = cpue_spde$n.spde)

# stack_cpue = inla.stack(
#   data = mako$CPUE,
#   A = list(A_cpue_mako, 1),
#   effects = list(list(i = 1:cpue_spde$n.spde), #spatial random effect
#                  list(m = rep(1, nrow(mako)))), #intercept?
#   tag = "mako_cpue"
# )

stack_cpue = inla.stack(
  data = mako$CPUE,
  A = list(A_cpue_mako, 1),
  effects = list(i = i_index, #spatial random effect
                 m = rep(1, nrow(mako))), #intercept?
  tag = "mako_cpue"
)




# re-try --------------------------------------------------------
## ----opts, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
# source('R/initial_setup.R')
# opts_chunk$set(
#   fig.path = 'figs/barrier-'
# )
library(scales)
library(rgeos)
## High resolution maps when using map()
library(mapdata) 
## Map features, map2SpatialPolygons()
library(maptools)
require(tidyverse)
require(INLA)

setwd("/Users/Yuki/Dropbox/Network2020")
data = read.csv("VASTdata.csv")
head(data, 2)
mako = data %>% filter(FISH == "makogarei", Y == 2018)

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

## ------------------------------------------------------------------------
# Define UTM projection
kmproj <- CRS("+proj=utm +zone=53 ellps=WGS84 +units=km")
# Project data
poly.water = spTransform(poly.water, kmproj)
pl.sel = spTransform(pl.sel, kmproj)
map.sp = spTransform(map.sp, kmproj)

## ------------------------------------------------------------------------
max.edge = 3
bound.outer = 5
mesh <- inla.mesh.2d(boundary = poly.water,
                     max.edge = c(1, 2) * max.edge,
                     cutoff = 2,
                     offset = c(max.edge, bound.outer))
plot(mesh)

water.tri = inla.over_sp_mesh(poly.water, y = mesh, 
                              type = "centroid", ignore.CRS = TRUE)
num.tri = length(mesh$graph$tv[, 1])
barrier.tri = setdiff(1:num.tri, water.tri)
poly.barrier = inla.barrier.polygon(mesh, 
                                    barrier.triangles = barrier.tri)
plot(poly.barrier)

plot(mesh, lwd = 0.5, add = FALSE)
plot(pl.sel, add = TRUE)
plot(map.sp, add = TRUE, col = gray(.9))
plot(poly.barrier, border = "red", add = TRUE)


plot(mesh, asp = 1)
points(x = mako$Lon,
       y = mako$Lat,
       col = 2,
       pch = 16, 
       cex = 0.5) #not run


range <- 200
barrier.model <- inla.barrier.pcmatern(mesh, 
                                       barrier.triangles = barrier.tri)
Q <- inla.rgeneric.q(barrier.model, "Q", theta = c(0, log(range)))

loc.corr <- c(10, 10)
corr <- book.spatial.correlation(Q, loc = loc.corr, mesh)

par(mfrow = c(1, 2), mar = c(0, 0, 0, 2), mgp = c(1, 0.5, 0), las = 1)
book.plot.field(corr, mesh = mesh, poly = poly.barrier) 


# zone = 53
# LongLatToUTM<-function(x,y,zone){
#   require(sp)
#   xy <- data.frame(ID = 1:length(x), X = x, Y = y)
#   coordinates(xy) <- c("X", "Y")
#   proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
#   res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
#   return(as.data.frame(res))
# }
# test = LongLatToUTM(mako$Lon, mako$Lat, zone) #合ってるのか．．．？

cpue_mako_lonlat = SpatialPolygons(list(Polygons(list(Polygon(
  cbind(mako$Lon, # x-axis 
        mako$Lat), # y-axis
  FALSE)), '0')), proj4string = CRS(proj4string(map.sp))) # Formal class SpatialPolygons
cpue_mako_lonlat2 = spTransform(cpue_mako_lonlat, kmproj) # Formal class SpatialPolygons
lonlat = cpue_mako_lonlat2@polygons[[1]]@Polygons[[1]]@coords

A.data <- inla.spde.make.A(mesh, lonlat)
