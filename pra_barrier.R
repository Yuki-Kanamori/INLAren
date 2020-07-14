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
  proj4string = CRS("+proj=longlat +datum=WGS84"))
summary(map.sp)

# make a polygon ------------------------------------------------------------
pl.sel <- SpatialPolygons(list(Polygons(list(Polygon(
  cbind(c(139.7, 139.5, 139.7, 140.1, 140.3, 139.9), # x-axis 
        c(35.2,  35.4,  35.8,  35.8,  35.4,  35.2)), # y-axis
  FALSE)), '0')), proj4string = CRS(proj4string(map.sp)))
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


## ------------------------------------------------------------------------
water.tri = inla.over_sp_mesh(poly.water, y = mesh, 
                              type = "centroid", ignore.CRS = TRUE)
num.tri = length(mesh$graph$tv[, 1])
barrier.tri = setdiff(1:num.tri, water.tri)
poly.barrier = inla.barrier.polygon(mesh, 
                                    barrier.triangles = barrier.tri)
plot(poly.barrier)


## ----label = "plot-barr-mesh2", fig = TRUE, echo = FALSE, fig.align = "center", fig.width = 6, heigh = 4.5, width = '97%', fig.cap = "The mesh constructed both over water and land. The grey region is the original land map. The inner red outline marks the coastline barrier."----

plot(mesh, lwd = 0.5, add = FALSE)
plot(pl.sel, add = TRUE)
plot(map.sp, add = TRUE, col = gray(.9))
plot(poly.barrier, border = "red", add = TRUE)



# North Japan -----------------------------------------------------------
# Select region 
map <- map("world", "Japan", fill = TRUE,
           col = "transparent", plot = TRUE)
IDs <- sapply(strsplit(map$names, ":"), function(x) x[1])
map.sp <- map2SpatialPolygons(
  map, IDs = IDs,
  proj4string = CRS("+proj=longlat +datum=WGS84"))
summary(map.sp)




# make a polygon ------------------------------------------------------------
pl.sel <- SpatialPolygons(list(Polygons(list(Polygon(
  cbind(c(135, 137, 139, 143, 144, 144), # x-axis 
        c(35,  39,  42,  42,  38,  36.2)), # y-axis
  FALSE)), '0')), proj4string = CRS(proj4string(map.sp)))
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
mesh.not <- inla.mesh.2d(boundary = poly.water, max.edge = 30,
                         cutoff = 2)

## ----label = "plot-barr-mesh1", fig = TRUE, echo = FALSE, fig.align = "center", fig.width = 10, heigh = 4.5, width = '97%', fig.cap = "The left plot shows the polygon for land in grey and the manually constructed polygon for our study area in light blue. The right plot shows the simple mesh, constructed only in the water."----
par(mfrow = c(1, 2), mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.7, 0), las = 1)
par(mar = c(0, 0, 0, 0))

plot(pl.sel, col = alpha("skyblue", 0.5), asp = 1)
plot(map.sp, add = TRUE, col = alpha(gray(0.9), 0.5))

plot(pl.sel, asp = 1)
plot(map.sp, add = TRUE, col = gray(0.9))
plot(mesh.not, add = TRUE)

## ------------------------------------------------------------------------
max.edge = 30
bound.outer = 80
mesh <- inla.mesh.2d(boundary = poly.water,
                     max.edge = c(1, 5) * max.edge,
                     cutoff = 1,
                     offset = c(max.edge, bound.outer))
plot(mesh)


## ------------------------------------------------------------------------
water.tri = inla.over_sp_mesh(poly.water, y = mesh, 
                              type = "centroid", ignore.CRS = TRUE)
num.tri = length(mesh$graph$tv[, 1])
barrier.tri = setdiff(1:num.tri, water.tri)
poly.barrier = inla.barrier.polygon(mesh, 
                                    barrier.triangles = barrier.tri)
plot(poly.barrier)


## ----label = "plot-barr-mesh2", fig = TRUE, echo = FALSE, fig.align = "center", fig.width = 6, heigh = 4.5, width = '97%', fig.cap = "The mesh constructed both over water and land. The grey region is the original land map. The inner red outline marks the coastline barrier."----

plot(mesh, lwd = 0.5, add = FALSE)
plot(pl.sel, add = TRUE)
plot(map.sp, add = TRUE, col = gray(.9))
plot(poly.barrier, border = "red", add = TRUE)

## ---- warning = FALSE, message = FALSE-----------------------------------
range <- 200
barrier.model <- inla.barrier.pcmatern(mesh, 
                                       barrier.triangles = barrier.tri)
Q <- inla.rgeneric.q(barrier.model, "Q", theta = c(0, log(range)))

## ---- warning = FALSE, message = FALSE-----------------------------------
stationary.model <- inla.spde2.pcmatern(mesh, 
                                        prior.range = c(1, 0.1), prior.sigma = c(1, 0.1))
Q.stat <- inla.spde2.precision(stationary.model, 
                               theta = c(log(range), 0))

## ---- warning = FALSE, message = FALSE-----------------------------------
# The location we find the correlation with respect to
summary(map.sp)
summary(pl.sel)
loc.corr <- c(1065, 4510)
corr <- book.spatial.correlation(Q, loc = loc.corr, mesh)
corr.stat <- book.spatial.correlation(Q.stat, loc = loc.corr,
                                      mesh)

## ----label = "plot-canada-b-corr", echo = FALSE, fig.width = 14, fig.heigh = 3.5, fig.cap = "The left plot shows the correlation structure of the Barrier model, with respect to the black point, while the right plot shows the correlation structure of the stationary model."----

par(mfrow = c(1, 2), mar = c(0, 0, 0, 2), mgp = c(1, 0.5, 0), las = 1)
book.plot.field(corr, mesh = mesh, poly = poly.barrier, 
                xlim = c(500, 1300), ylim = c(3873, 4700), zlim = c(0.1, 1)) 
points(loc.corr[1], loc.corr[2], pch = 19)
book.plot.field(corr.stat, mesh = mesh, poly = poly.barrier, 
                xlim = c(500, 1300), ylim = c(3873, 4700), zlim = c(0.1, 1)) 
points(loc.corr[1], loc.corr[2], pch = 19)

