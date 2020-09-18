### Code for Gompertz state-space model ###

#The observational equation is defined as 
# y_{j,t} \sim NegBinom(\mu_{j,t}, \phi)
# And the process model is defined as 
# Z_{j,t}= theta Z_{j,t-1} + beta  X_{j,t} + epsilon_{t}
# where theta is an AR1 parameter for density dependence, 
# beta is the regression coefficient and intercept parameter 
# comprising the growth rate, and epsilon_t is latent error.  
# The parameter theta does not vary for each {j}, 
# and is listed as the parameter "Beta for ix1b" in the output 
# below. 

# Below is example code for a state-space model with snow cover 
# extent.

# Importing and arranging the count data and covariates to
# set the process model equal to 0, or
# 0=\theta Z_{j,t-1}-Z_{j,t}+\beta_{0,j}+\beta X_{t}+ \mathtt{offset_{j,t}}+ \epsilon.

# As shown in Ruiz-Cardenas et al. (2010), 
# we set up our data matrix into two columns.  
# The first column contains n*k rows of count data followed 
# by NAs, while the second column contains NAs for the 
# first n*k rows followed by 0s for n*k-n rows. 

require(INLA)

ss.data <- read.csv("ss_model_inla.csv",header=TRUE)
head(ss.data)
n=length(unique(ss.data$stratum)) #number of strata (study sites) in study
k=length(unique(ss.data$year)) #number of years
y=matrix(as.vector(ss.data$y),nrow=6,ncol=43) #matrix of counts

# calculate offset parameter for varying sampling intensity
N.offset1=ss.data$N
N.off.tmp=matrix(N.offset1,n,k)
N.off.1=matrix(NA,n,k)
for(i in 1:n){
N.off.1[i,]=N.off.tmp[i,]-min(N.off.tmp[i,])
}
N.offset2=N.off.1[,-1]
N.offset=c(rep(NA,n*k),t(N.offset2))
N.offset=log(N.offset+1)

## Set up data matrix
nd <- n*k
Y <- matrix(NA, nd*2-n, 2)
Y[1:nd             , 1] <- as.vector(t(y)) #counts in first column
Y[1:(nd-n) + nd    , 2] <- 0 #Z_t in second column

head(Y)
Y[(nd-1):(nd+2),]

id1  <- (1:nd)[-((1:n)*k)]
id2  <- (1:nd)[-c(1,((1:(n-1))*k)+1)]
ix1  <- c(1:nd, id2)                  ## indices for Z_t 
ix1b <- c(rep(NA,nd), id1)            ## indices for Z_{t-1} 
wx1b <- c(rep(1,nd), rep(-1,nd-n))   ## weights for Z_{t-1} 
iw22 <- c(rep(NA,nd),id2)     ## indices for w_t 
st.1=c(rep(NA,nd),rep(1,(k-1)),rep(0,(nd-6-(k-1)))) #intercepts 
#for each growth rate for each stratum 
st.2=c(rep(NA,nd),rep(0,(k-1)),rep(1,(k-1)),rep(0,(nd-6-((k-1)*2))))
st.3=c(rep(NA,nd),rep(0,((k-1)*2)),rep(1,(k-1)),rep(0,(nd-6-((k-1)*3))))
st.4=c(rep(NA,nd),rep(0,((k-1)*3)),rep(1,(k-1)),rep(0,(nd-6-((k-1)*4))))
st.5=c(rep(NA,nd),rep(0,((k-1)*4)),rep(1,(k-1)),rep(0,(nd-6-((k-1)*5))))
st.6=c(rep(NA,nd),rep(0,((k-1)*5)),rep(1,(k-1)),rep(0,(nd-6-((k-1)*6))))
snow.idx <- matrix(ss.data$stand.snow[7:258],6,42) #rearranging snow data
snow <- c(rep(NA,nd),c(snow.idx[1,],snow.idx[2,],snow.idx[3,],
snow.idx[4,],snow.idx[5,],snow.idx[6,]))

dat.snow=list(Y=Y, ix1=ix1,ix1b=ix1b,wx1b=wx1b,iw22=iw22,snow=snow,
st.1=st.1,st.2=st.2,st.3=st.3,st.4=st.4,
st.5=st.5,st.6=st.6,N.offset=N.offset) #collection of all the data

# Specifing the model for the state-space model and 
# calling the inla() funciton

formula.snow <- Y ~ f(iw22,model="iid") + 
  f(ix1,wx1b, model="iid",initial=-10, fixed=TRUE) + 
  f(ix1b, copy="ix1",fixed=FALSE) + N.offset +
  st.1 + st.2 + st.3 + st.4 + st.5 + st.6 + snow -1

r.snow <- inla(formula.snow, data = dat.snow,
  family = c("nbinomial","gaussian"),
  control.family = list(list(link="log"),
  list(initial=10, fixed=TRUE)),
  control.predictor=list(compute=TRUE),
  control.compute=list(dic=TRUE,cpo=TRUE))

#snow.recalc <- inla.cpo(r.snow,force=TRUE) #recalculate CPO values

summary(r.snow)
r.snow$summary.ran$ix1[43,]

# Plotting the abundance for a single spatial unit (e.g., Stratum 13, id=1).

rang <- range(c(exp(r.snow$summary.ran$ix1[1:k,4]),exp(r.snow$summary.ran$ix1[1:k,6])))
plot(exp(r.snow$summary.ran$ix1[1:k,2]), type="l", 
ylim=rang, col="red", xlim=c(1,k),ylab="Abundance of Scaup Pairs",xlab="time",
main="Stratum 13")
lines(exp(r.snow$summary.ran$ix1[1:k,4]), col="blue", lty=3,lwd=2)
lines(exp(r.snow$summary.ran$ix1[1:k,6]), col="blue", lty=3,lwd=2)
points(y[1,],col="black")
legend("topleft",legend=c("observed","posterior mean","95% CI"), 
col=c("black", "black","blue"),lty=c(NA,1,2),
pch=c(1,NA,NA),bty="n")

# Plotting the difference in abundance from Z_{j,t-1} to Z_{j,t} in Stratum 13 as snow changes. 

x.min=min(snow[259:300])
x.max=max(snow[259:300])
x.snow=seq(x.min,x.max,length=100)
#Intercept for stratum 13 + snow
y.pred=exp(.7968 +.1358*x.snow)
plot(y=y.pred,x=x.snow,cex.lab=1.75,cex.axis=1.5,
xlab="Snow Cover Extent",
ylab="Change in Scaup Pairs",type="l",
lwd=1.5,lty=1,mgp=c(2.5,1,0))

# To determine the change in the predicted value of y_{j,43} 
# given y_{j,42}=93 and the predicted value of 
# standardized snow cover extent is the mean for stratum 13 (0.55), 
# we can estimate y_{j,43} as

lc1 = inla.make.lincomb(snow=0.55,st.1=1,
  st.2=0,st.3=0,st.4=0,st.5=0,st.6=0,
  ix1=c(rep(NA,42),1,rep(NA,215)),
  ix1b=c(rep(NA,41),-1,rep(NA,216)),
  N.offset=1.09,iw22=c(rep(NA,42),1,rep(NA,215)))

r.snow.predict <- inla(formula.snow, data = dat.snow,
  lincomb=lc1,
  family = c("nbinomial","gaussian"),
  control.family = list(list(link="log"),
  list(initial=10, fixed=T)),
  control.predictor=list(compute=TRUE),
  control.compute=list(dic=TRUE,cpo=TRUE),
  control.inla=list(lincomb.derived.only=FALSE))

lc1.2= r.snow.predict$marginals.lincomb.derived$lc

r.snow.predict$summary.lincomb.derived

#Distribution of predicted value on log scale
plot(lc1.2)

#transform posterior of linear combination
plot(inla.tmarginal(exp,lc1.2),xlab="Change in Scaup Pairs",ylab="Density") 

#quantiles of linear combinations
inla.qmarginal(c(0.025,0.975),inla.tmarginal(exp,lc1.2)) 
