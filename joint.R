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

# spde
spde = inla.spde2.pcmatern(mesh = mesh2, alpha = 2, prior.range = c(0.01, 0.05), prior.sigma = c(1, 0.01))

# stack for eDNA
e_stk = inla.stack(data = list(y = cbind(edna, NA)),
                   A = list(e_A, 1),
                   effects = list(i.e = 1:mesh2$n, eb.0 = rep(1, length(edna))),
                   tag = "e_dat")
c_stk = inla.stack(data = list(y = cbind(NA, catch)),
                   A = list(c_A, 1),
                   effects = list(list(i.c = 1:mesh2$n, x = 1:mesh2$n), cb.0 = rep(1, length(catch))),
                   tag = "c_dat")

stk = inla.stack(e_stk, c_stk)

# formula
formula = y ~ 0 + eb.0 + cb.0 + f(i.e, model = spde) + f(x, model = spde) + f(i.c, copy = "i.e", fixed = FALSE)

# fitting the joint model
res_joint = inla(formula, data = inla.stack.data(stk), family = c("binomial", "binomial"), control.predictor = list(compute = TRUE, A = inla.stack.A(stk)))



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
data.inla.projector <- inla.mesh.projector(mesh2, loc = loc)
newdata <- data.frame(loc, mean = inla.mesh.project(data.inla.projector, res_joint$summary.random$i.e$mean + res_joint$summary.random$x$mean + res_joint$summary.random$i.c$mean) + res_joint$summary.fixed$mean[1])
str(newdata$mean)
ggplot(newdata, aes(y = X2, x = X1)) + geom_tile(aes(fill = mean))
ggplot(loc, aes(y = V1, x = V2)) + geom_point(aes(color = y), size = 2)


p_test = inla.mesh.projector(mesh2, xlim = 139.5:140.5, ylim = 35:36)
pred_mean = inla.mesh.project(p_test, res_joint$summary.random$i.e$mean)
plot(pred_mean)
