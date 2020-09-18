library(INLA)
#library(INLAutils)
library(tidyverse)
library(RColorBrewer)


#Make data set*************************************************************************************
load("/Users/00007269/保存/■2020/06_IOTC_SWO/2_figs/d_glm.Rdata")

#d0 = d_glm %>% filter(area=="NW",yr>=1994) %>%
#              mutate(yr=as.factor(yr),
#                     qtr=as.factor(qtr),
#                     latlon=as.factor(latlon),
#                    jp_name=as.factor(jp_name))

#Aggregate data-set to reduce calculation time
d = d_glm %>% filter(area == "NW",
                     yr >= 1994) %>%
  group_by(yr,
           month,
           jp_name,
           lat,
           lon,
           hpb) %>%
  summarise(swo = sum(swo),
            hooks = sum(hooks)) %>%
  ungroup() %>%
  mutate(cpue = swo/(hooks/1000),
         hpb = as.factor(hpb),
         year = as.numeric(yr),
         yr = as.factor(yr),
         qtr = as.factor(floor((as.numeric(month)-1)/3)+1),
         jp_name = as.factor(jp_name),
         lat5 = 5*floor(lat/5),
         lon5 = 5*floor(lon/5),
         latlon = as.factor(paste(lat5,lon5,sep=""))) 

#g = ggplot(data=d, aes(x=factor(yr),y=cpue))
#g = g + geom_boxplot()
#g

#Set a prior for random effect (idd) model*********************************************************
halfcauchy = "expression:
              lambda = 0.022;
              precision = exp(log_precision);
              logdens = -1.5*log_precision-log(pi*lambda)-log(1+1/(precision*lambda^2));
              log_jacobian = log_precision;
              return(logdens+log_jacobian);"

hcprior = list(prec = list(prior = halfcauchy))

#Null model****************************************************************************************
m_null = inla(swo ~ 1, 
              data = d,
              offset = log(d$hooks/1000),
              family = "poisson", 
              control.family = list(link="log"),
              control.predictor = list(link=1, compute=TRUE),
              control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

#summary(m_null)

#Traditional GLM***********************************************************************************

m_glm = inla(swo ~ yr + qtr + latlon, 
             data = d,
             offset = log(d$hooks/1000),
             family = "poisson", 
             control.family = list(link="log"),
             control.predictor = list(link=1, compute=TRUE),
             control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

#summary(m_glm)

#Simple Poisson GLMM*******************************************************************************
#When d0 was used, it needs four hours!! Result was not different between d0 and d (includes big catch).
#prec.prior = list(prec = list(param = c(0.001, 0.001)))
#swo ~ yr + qtr + f(latlon, model="iid", hyper=prec.prior) + f(jp_name, model="iid", hyper=prec.prior) #Use prior

m_glmm = inla(swo ~ yr + qtr + f(latlon, model="iid", hyper=hcprior) + 
                f(jp_name, model="iid") + f(hpb, model="iid"), 
              data = d,
              offset = log(d$hooks/1000),
              family = "poisson", 
              control.family = list(link="log"),
              control.predictor = list(compute=TRUE),
              control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

#summary(m_glmm)
#m_glmm$summary.random
#m_glmm$summary.fixed


#ZIP GLMM******************************************************************************************
m_zip_glmm = inla(swo ~ yr + qtr + f(latlon, model="iid") + f(jp_name, model="iid"), 
                  data = d,
                  offset = log(d$hooks/1000),
                  family = "zeroinflatedpoisson1", 
                  control.family = list(link="log"),
                  control.predictor = list(compute=TRUE),
                  control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

#summary(m_zip_glmm)

#Predict standardized CPUE

#Spatial Poisson GLMM******************************************************************************
#Step 1: Make a mesh
loc = as.matrix(bind_cols(lon=d$lon,lat=d$lat))
mesh = inla.mesh.2d(loc=loc, max.edge=c(5,5), cutoff=1.5, offset=c(2,2))
plot(mesh, asp=1, main = "")
points(loc, col = "red", pch = 16, cex = .5)
mesh$n

#Step 2: Define the  weighting factors a_ik
A = inla.spde.make.A(mesh=mesh,loc=loc)
dim(A)

#Step 3: Define spatial field (Make spde matrix)
spde = inla.spde2.matern(mesh, alpha=2)

#Step 4: Make a stack
w.index = inla.spde.make.index(name="w",
                               n.spde=spde$n.spde,
                               n.group=1)

X = data.frame(intercept=rep(1, nrow(d)), yr=d$yr, qtr=d$qtr)

StackFit = inla.stack(tag = "Fit",
                      data = list(swo=d$swo),  
                      A = list(1, 1, 1, A),                  
                      effects = list(X=X, hpb=d$hpb, jp_name=d$jp_name, w=w.index))

#str(StackFit)

#Step 5: Estimates model parameters
m_spde = inla(swo ~ 0 + intercept + yr + qtr + f(hpb, model="iid") + f(jp_name, model="iid") + f(w, model=spde), 
              data = inla.stack.data(StackFit),
              offset = log(d$hooks/1000),
              family = "poisson", 
              control.family = list(link="log"),
              control.predictor = list(A=inla.stack.A(StackFit), compute=TRUE),
              control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

#summary(m_spde)

#Spatial Poisson GLMM-2****************************************************************************

StackFit2 = inla.stack(tag = "Fit",
                       data = list(swo=d$swo),  
                       A = list(1, 1, 1, 1, 1, A),                  
                       effects = list(intercept=rep(1, nrow(d)), yr=d$year, month=factor(d$month), 
                                      hpb=d$hpb, jp_name=d$jp_name, w=w.index))

m_spde2 = inla(swo ~ 0 + intercept + f(yr, model="ar1") + f(month, model="iid", hyper=hcprior) + f(hpb, model="iid", hyper=hcprior) +
                 f(hpb, model="iid", hyper=hcprior) + f(jp_name, model="iid", hyper=hcprior) + f(w, model=spde), 
               data = inla.stack.data(StackFit2),
               offset = log(d$hooks/1000),
               family = "poisson", 
               control.family = list(link="log"),
               control.predictor = list(A=inla.stack.A(StackFit2), compute=TRUE),
               control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

summary(m_spde2)
m_spde2$summary.random
m_spde2$summary.fixed

#Randomized quantile residuals using posterior mean
#res_resid = data.frame(fitted=m_spde2$summary.fitted.values[c(1:length(d$swo)),], obs=d$swo) %>% mutate(res = obs-fitted.mean)

#Fx = ppois(res_resid$obs, res_resid$fitted.mean) 
#px = dpois(res_resid$obs, res_resid$fitted.mean)
#u = runif(length(Fx))
#pvalue = Fx-u*px
#pvalue = pmin(pmax(pvalue,10^{-10}),1-10^{-10})
#rqr = qnorm(pvalue)

#res_resid = res_resid %>% mutate(rqr=rqr) #%>% filter(rqr!=Inf, rqr!=-Inf)

#g = ggplot()
#g = g + theme_light()
#g = g + geom_point(data=res_resid, aes(x=fitted.mean, y=rqr),alpha=0.2)
#g = g + xlab("Fitted mean value") + ylab("Randomized quantile residuals")
#g = g + labs(title = "Randomized quantile residuals")
#g

#Plot  posterior probability distribution
#fixed_par = data.frame(m_spde2$marginals.fixed[[1]]) %>% mutate(par=names(m_spde2$marginals.fixed[1]))

#tmp = 0

#for (i in 1:7) {
#  tmp = rbind(tmp,(data.frame(m_spde2$marginals.hyperpar[[i]]) %>% mutate(par=names(m_spde2$marginals.hyperpar)[i])))
#}

#hyp_par = tmp[-1,]
#rm(tmp)

#pars=rbind(fixed_par, hyp_par)

#g = ggplot()
#g = g + theme_light()
#g = g + geom_line(data=pars,aes(x=x,y=y))
#g = g + facet_wrap(~ par, scales = "free", ncol=2)
#g = g + xlab("Estimated value") + ylab("Posterior density")
#g

#Predict for latent spatial field

#lon_range = diff(range(d$lon))
#lat_range = diff(range(d$lat))
#step_size = 0.5
#nlonlat = round(c(lon_range, lat_range)/step_size)

#projgrid = inla.mesh.projector(mesh, xlim=range(d$lon), ylim=range(d$lat), dims=nlonlat)

#xmean = inla.mesh.project(projgrid, m_spde2$summary.random$w$mean)
#xsd = inla.mesh.project(projgrid, m_spde2$summary.random$w$sd)

#locpred = as.data.frame(projgrid$lattice$loc) %>% rename(Longitude=V1, Latitude=V2) %>%
#          mutate(xmean = as.vector(xmean)) %>%
#          filter(!is.na(xmean))

#col = rev(brewer.pal(11, "Spectral"))
#world.map = map_data("world")

#g = ggplot()
#g = g + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank())
#g = g + theme(panel.border=element_rect(fill=NA,colour="gray",size=0.75))
#g = g + theme(axis.ticks=element_line(colour="black"),axis.text=element_text(colour = "black"))
#g = g + geom_raster(data=locpred,aes(x=Longitude,y=Latitude,fill=xmean))
#g = g + geom_polygon(data=world.map,aes(x=long,y=lat,group=group),fill="black")
#g = g + scale_fill_gradientn(colours=col,breaks=c(-1,-0.5,0,0.5,1), labels=c("-1","-0.5","0","0.5","1"))
#g = g + coord_fixed(xlim=c(25,110),ylim=c(-40,20))
#g = g + scale_x_continuous(breaks=seq(25,110,25),labels=c("25E","50E","75E","100E"))
#g = g + xlab("Longitude") + ylab("Latitude")
#g = g + labs(fill="Latent \n Spatial \n Feild")
#g

#Calculate standardized CPUE

#pred.month = m_spde2$summary.random$month %>% mutate(mean2=mean^2) %>% filter(mean2==min(mean2))
#pred.hpb = m_spde2$summary.random$hpb %>% mutate(mean2=mean^2) %>% filter(mean2==min(mean2)) 
#pred.jpname = m_spde2$summary.random$jp_name %>% mutate(mean2=mean^2) %>% filter(mean2==min(mean2))
#pred.loc = m_spde2$summary.random$w %>% mutate(mean2=mean^2) %>% filter(mean2==min(mean2))

#d.pred = expand.grid(lon=mesh$loc[pred.loc$ID,1], lat=mesh$loc[pred.loc$ID,2], year=c(rep(1994:2010),rep(2012:2018)), month=pred.month$ID, hpb=pred.hpb$ID,
#                     jp_name=pred.jpname$ID, hooks=mean(d$hooks)) %>% 
#         mutate(lon=as.numeric(lon),
#                lat=as.numeric(lat),
#               year=as.numeric(year),
#               month=as.factor(month),
#               hpb=as.factor(hpb),
#               jp_name=as.factor(jp_name),
#               hooks=as.numeric(hooks))

#A.pred = inla.spde.make.A(mesh=mesh, loc=as.matrix(d.pred[,c(1:2)]))

#StackPred2 = inla.stack(tag = "Pred",
#                        data = list(swo=NA),  
#                        A = list(1,1,1,1,1,A.pred),                  
#                        effects = list(intercept=rep(1, nrow(d.pred)), yr=d.pred$year,
#                                       month=d.pred$month, hpb=d.pred$hpb, jp_name=d.pred$jp_name, w=w.index))

#join.stack = inla.stack(StackFit2, StackPred2)

#joint.output = inla(swo ~ 0 + intercept + f(yr, model="ar1") + f(month, model="iid") + f(hpb, model="iid") +
#                    f(jp_name, model="iid") + f(w, model=spde), 
#                   data = inla.stack.data(join.stack),
#                    offset = c(log(d$hooks/1000), log(d.pred$hooks/1000)),
#                    family = "poisson", 
#                    control.family = list(link="log"),
#                    control.predictor = list(A=inla.stack.A(join.stack), compute=TRUE),
#                    control.results = list(return.marginals.random = FALSE, return.marginals.predictor = FALSE), 
#                    control.mode = list(theta = m_spde2$mode$theta, restart = FALSE),
#                    verbose=FALSE)

#index.pred.response = inla.stack.index(join.stack,tag='Pred')$data
#nominal_cpue = d %>% group_by(year) %>% summarise(swo=sum(swo),hooks=sum(hooks)) %>% mutate(nom_cpue=swo/(hooks/1000))
#post.mean.pred.resonse = joint.output$summary.linear.predictor[index.pred.response,] 
#tmp = cbind(nominal_cpue,post.mean.pred.resonse) %>% mutate(std_cpue=exp(mean-log(mean(d$hooks/1000))),
#                                                            upper=exp(post.mean.pred.resonse[,5]-log(mean(d$hooks/1000))),
#                                                            lower=exp(post.mean.pred.resonse[,3]-log(mean(d$hooks/1000))))
#res_fig = rbind(tmp,c(2011,rep(NA,13)))

#g = ggplot(data=res_fig)
#g = g + geom_ribbon(aes(x=year, ymin = lower, ymax = upper), fill = "grey70")
#g = g + geom_point(aes(x=year,y=nom_cpue))
#g = g + geom_line(aes(x=year,y=std_cpue))
#g = g + ylim(c(0,1))
#g

#Spatial ZIP GLMM**********************************************************************************
#Spatial assumption and data sets are same as "Spatial GLMM"
m_zip_spde = inla(swo ~ 0 + intercept + yr + qtr + f(hpb, model="iid") + f(jp_name, model="iid") + f(w, model=spde), 
                  data = inla.stack.data(StackFit),
                  offset = log(d$hooks/1000),
                  family = "zeroinflatedpoisson1", 
                  control.family = list(link="log"),
                  control.predictor = list(A=inla.stack.A(StackFit), compute=TRUE),
                  control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

#summary(m_zip_spde)

#Spatial ZIP GLMM-2********************************************************************************

StackFit2 = inla.stack(tag = "Fit",
                       data = list(swo=d$swo),  
                       A = list(1, 1, 1, 1, 1, A),                  
                       effects = list(intercept=rep(1, nrow(d)), yr=d$year, month=factor(d$month), 
                                      hpb=d$hpb, jp_name=d$jp_name, w=w.index))

m_zip_spde2 = inla(swo ~ 0 + intercept + f(yr, model="ar1") + f(month, model="iid", hyper=hcprior) + f(hpb, model="iid", hyper=hcprior) +
                     f(hpb, model="iid", hyper=hcprior) + f(jp_name, model="iid", hyper=hcprior) + f(w, model=spde), 
                   data = inla.stack.data(StackFit2),
                   offset = log(d$hooks/1000),
                   family = "zeroinflatedpoisson1", 
                   control.family = list(link="log"),
                   control.predictor = list(A=inla.stack.A(StackFit2), compute=TRUE),
                   control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

summary(m_zip_spde2)
#m_zip_spde2$summary.random
#m_zip_spde2$summary.fixed

#Model selection using WAIC************************************************************************

waics = data.frame(model=c("m_null","m_glm","m_glmm","m_zip_glmm","m_spde","m_spde2","m_zip_spde","m_zip_spde2"),
                   waic=c(m_null$waic$waic,m_glm$waic$waic,m_glmm$waic$waic,m_zip_glmm$waic$waic,m_spde$waic$waic,
                          m_spde2$waic$waic,m_zip_spde$waic$waic,m_zip_spde2$waic$waic))

#Results of selected model (m_zip_spde2)***********************************************************
#Randomized quantile residuals using posterior mean
res_resid = data.frame(fitted=m_zip_spde2$summary.fitted.values[c(1:length(d$swo)),], obs=d$swo) %>% 
  mutate(res = obs-fitted.mean)

dzpois = function(x,lambda,p){return((1-p)*dpois(x,lambda)+p*(x==0))}
pzpois = function(x,lambda,p){return((1-p)*ppois(x,lambda)+p*(x>=0))} 
theta = m_zip_spde2$summary.hyperpar[1,1]
p =  rep(exp(theta)/(1+exp(theta)),length(res_resid$obs))

Fx = pzpois(res_resid$obs-1, res_resid$fitted.mean, p) 
px = dzpois(res_resid$obs, res_resid$fitted.mean, p)
u = runif(length(Fx))
pvalue = Fx+u*px
pvalue = pmin(pmax(pvalue,10^{-10}),1-10^{-10})
rqr = qnorm(pvalue)

res_resid = res_resid %>% mutate(rqr=rqr) %>% filter(rqr!=Inf, rqr!=-Inf)

g = ggplot()
g = g + theme_light()
g = g + geom_point(data=res_resid, aes(x=fitted.mean, y=rqr),alpha=0.2)
g = g + xlab("Fitted mean value") + ylab("Randomized quantile residuals")
g = g + labs(title = "Randomized quantile residuals")
g

#Plot  posterior probability distribution
fixed_par = data.frame(m_zip_spde2$marginals.fixed[[1]]) %>% mutate(par=names(m_zip_spde2$marginals.fixed[1]))

tmp = 0

for (i in 1:7) {
  tmp = rbind(tmp,(data.frame(m_zip_spde2$marginals.hyperpar[[i]]) %>% mutate(par=names(m_spde2$marginals.hyperpar)[i])))
}

hyp_par = tmp[-1,]
rm(tmp)

pars=rbind(fixed_par, hyp_par)

g = ggplot()
g = g + theme_light()
g = g + geom_line(data=pars,aes(x=x,y=y))
g = g + facet_wrap(~ par, scales = "free", ncol=2)
g = g + xlab("Estimated value") + ylab("Posterior density")
g

#Predict of latent spatial field

lon_range = diff(range(d$lon))
lat_range = diff(range(d$lat))
step_size = 0.5
nlonlat = round(c(lon_range, lat_range)/step_size)

projgrid = inla.mesh.projector(mesh, xlim=range(d$lon), ylim=range(d$lat), dims=nlonlat)

xmean = inla.mesh.project(projgrid, m_zip_spde2$summary.random$w$mean)
xsd = inla.mesh.project(projgrid, m_zip_spde2$summary.random$w$sd)

locpred = as.data.frame(projgrid$lattice$loc) %>% rename(Longitude=V1, Latitude=V2) %>%
  mutate(xmean = as.vector(xmean)) %>%
  filter(!is.na(xmean))

col = rev(brewer.pal(11, "Spectral"))
world.map = map_data("world")

g = ggplot()
g = g + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank())
g = g + theme(panel.border=element_rect(fill=NA,colour="gray",size=0.75))
g = g + theme(axis.ticks=element_line(colour="black"),axis.text=element_text(colour = "black"))
g = g + geom_raster(data=locpred,aes(x=Longitude,y=Latitude,fill=xmean))
g = g + geom_polygon(data=world.map,aes(x=long,y=lat,group=group),fill="black")
g = g + scale_fill_gradientn(colours=col,breaks=c(-1,-0.5,0,0.5,1), labels=c("-1","-0.5","0","0.5","1"))
g = g + coord_fixed(xlim=c(25,110),ylim=c(-40,20))
g = g + scale_x_continuous(breaks=seq(25,110,25),labels=c("25E","50E","75E","100E"))
g = g + xlab("Longitude") + ylab("Latitude")
g = g + labs(fill="Latent \n Spatial \n Feild")
g

#Calculate standardized CPUE
pred.month = m_zip_spde2$summary.random$month %>% mutate(mean2=mean^2) %>% filter(mean2==min(mean2))
pred.hpb = m_zip_spde2$summary.random$hpb %>% mutate(mean2=mean^2) %>% filter(mean2==min(mean2)) 
pred.jpname = m_zip_spde2$summary.random$jp_name %>% mutate(mean2=mean^2) %>% filter(mean2==min(mean2))
pred.loc = m_zip_spde2$summary.random$w %>% mutate(mean2=mean^2) %>% filter(mean2==min(mean2))

d.pred = expand.grid(lon=mesh$loc[pred.loc$ID,1], lat=mesh$loc[pred.loc$ID,2], year=rep(1994:2018), month=pred.month$ID, hpb=pred.hpb$ID,
                     jp_name=pred.jpname$ID, hooks=mean(d$hooks)) %>% 
  mutate(intercept=1,
         year=as.numeric(year),
         month=as.factor(month),
         hpb=as.factor(hpb),
         jp_name=as.factor(jp_name),
         hooks=as.numeric(hooks))

A.pred = inla.spde.make.A(mesh=mesh, loc=as.matrix(d.pred[,c(1:2)]))

StackPred2 = inla.stack(tag = "Pred",
                        data = list(swo=NA),  
                        A = list(1,1,1,1,A.pred),                  
                        effects = list(yr=d.pred$year, month=d.pred$month, hpb=d.pred$hpb, 
                                       jp_name=d.pred$jp_name, w=w.index))

join.stack = inla.stack(StackFit2, StackPred2)

joint.output = inla(swo ~ 0 + intercept + f(yr, model="ar1") + f(month, model="iid") + f(hpb, model="iid") +
                      f(jp_name, model="iid") + f(w, model=spde), 
                    data = inla.stack.data(join.stack),
                    offset = c(log(d$hooks/1000), log(d.pred$hooks/1000)),
                    family = "poisson", 
                    control.family = list(link="log"),
                    control.predictor = list(A=inla.stack.A(join.stack), compute=TRUE, link=1),
                    control.results = list(return.marginals.random = FALSE, return.marginals.predictor = FALSE))
#control.mode = list(theta = m_zip_spde2$mode$theta, restart = FALSE),
#verbose=FALSE)

index.pred.response = inla.stack.index(join.stack,tag='Pred')$data
nominal_cpue = d %>% group_by(year) %>% summarise(swo=sum(swo),hooks=sum(hooks)) %>% mutate(nom_cpue=swo/(hooks/1000))
post.mean.pred.resonse = joint.output$summary.linear.predictor[index.pred.response,]  %>% mutate(std_cpue=exp(mean-log(mean(d$hooks/1000))),
                                                                                                 upper=exp(post.mean.pred.resonse[,5]-log(mean(d$hooks/1000))),
                                                                                                 lower=exp(post.mean.pred.resonse[,3]-log(mean(d$hooks/1000))),
                                                                                                 year=rep(1994:2018))
res_fig = full_join(nominal_cpue,post.mean.pred.resonse)

g = ggplot(data=res_fig)
g = g + geom_ribbon(aes(x=year, ymin = lower, ymax = upper), fill = "grey70")
g = g + geom_point(aes(x=year,y=nom_cpue))
g = g + geom_line(aes(x=year,y=std_cpue))
g = g + xlab("Year") + ylab("CPUE (SWO/1000hooks)")
g = g + ylim(c(0,2.5))
g

#Spatio-temporal Poisson GLMM**********************************************************************
#tmp = data.frame("yr"="2011", "month"=NA, "jp_name"=NA, "lat"=NA, "lon"=NA, "hpb"=NA, "swo"=NA, 
#                 "hooks"=NA, "cpue"=NA, "year"=NA, "qtr"=NA, "lat5"=NA, "lon5"=NA, "latlon"=NA)
#d2 = rbind(d,tmp)
#Step 1: Make a mesh

#Step 2: Define the  weighting factors a_ik at time t
#At = inla.spde.make.A(mesh=mesh, loc=loc, group=as.numeric(d$yr))
#dim(At)

#Step 3: Define spatial field
#Same as spatial GLMM

#Step 4: Make a stack
#w.index2 = inla.spde.make.index(name="w",
#                               n.spde=spde$n.spde,
#                               n.group=24)

#Xm = model.matrix(~ -1  + yr + qtr, data=d) 
#N = nrow(d)
#colnames(Xm)

#X = as.data.frame(Xm[,-which(colnames(Xm)%in%c("yr1994"))]) 

#StackFit = inla.stack(tag = "Fit",
#                      data = list(swo=d$swo),  
#                      A = list(1, 1, 1, 1, 1, At),                  
#                      effects = list(intercept=rep(1,N), yr=d$yr, qtr=d$qtr, hpb=d$hpb, jp_name=d$jp_name, w=w.index2))

#Step 5: Estimates model parameters
#m_sp = inla(swo ~ 0 + intercept + yr + qtr + f(jp_name, model="iid") +
#                  f(hpb, model="iid") +
#                  f(w, model=spde, group=w.group, control.group=list(model="ar1")), 
#            data = inla.stack.data(StackFit),
#            offset = log(d$hooks/1000),
#            family = "poisson", 
#            control.family = list(link="log"),
#            control.predictor = list(A=inla.stack.A(StackFit), compute=TRUE),
#            control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

#summary(m_sp)

#Spatio-temporal ZIP GLMM**************************************************************************
#Spatio-temporal assumption is same as "Spatio-temporal Poisson GLMM"
#m_zip_sp = inla(swo ~ 0 + intercept + yr + qtr + f(jp_name, model="iid") +
#                f(hpb, model="iid") +
#                f(w, model=spde, group=w.group, control.group=list(model="ar1")), 
#                data = inla.stack.data(StackFit),
#                offset = log(d$hooks/1000),
#                family = "zeroinflatedpoisson1", 
#                control.family = list(link="log"),
#                control.predictor = list(A=inla.stack.A(StackFit), compute=TRUE),
#                control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

#summary(m_zip_sp)

#Spatio-temporal ZIP GLMM2*************************************************************************
#Spatio-temporal assumption is same as "Spatio-temporal Poisson GLMM"

#StackFit2 = inla.stack(tag = "Fit",
#                       data = list(swo=d$swo),  
#                       A = list(1, 1, 1, 1, 1, At),                  
#                       effects = list(intercept=rep(1, nrow(d)), yr=d$year, month=factor(d$month), 
#                                      hpb=d$hpb, jp_name=d$jp_name, w=w.index2))

#m_zip_sp2 = inla(swo ~ 0 + intercept + f(yr, model="ar1") + f(month, model="iid", hyper=hcprior) + f(hpb, model="iid", hyper=hcprior) +
#                                       f(hpb, model="iid", hyper=hcprior) + f(jp_name, model="iid", hyper=hcprior) + f(w, model=spde), 
#                data = inla.stack.data(StackFit2),
#                offset = log(d$hooks/1000),
#                family = "zeroinflatedpoisson1", 
#                control.family = list(link="log"),
#               control.predictor = list(A=inla.stack.A(StackFit2), compute=TRUE),
#               control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

