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



# higashi nihon -----------------------------------------------------------
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
  cbind(c(135, 137, 139, 142, 143, 141), # x-axis 
        c(35,  39,  42,  42,  38,  35.85)), # y-axis
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
max.edge = 30
bound.outer = 90
mesh <- inla.mesh.2d(boundary = poly.water,
                     max.edge = c(1, 3) * max.edge,
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

