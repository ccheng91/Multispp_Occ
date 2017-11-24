################################################################################
# R and WinBUGS code for running occupancy model reported in:
# "Occupancy modeling of a management indicator species: black-backed woodpeckers
# on burned Sierra Nevada forests" by Saracco, Siegel, and Wilkerson
#
# Last edited by J. Saracco on 1 Sept. 2010
################################################################################

library(R2WinBUGS)

setwd("C:/Work/Manuscripts/BBWO/")

bbwo.dat <- read.table("bbwo_dat.txt", header = T, sep = ",")

#create index for fires
fireID <- as.vector(bbwo.dat[,"fire"])
nfire <- length(unique(fireID))
fire.ind <- rep(NA, nfire)
for (i in 1:nfire){
  fire.ind[fireID == unique(fireID)[i]] <- i
}

# standardized fire age covariate
fa <- (bbwo.dat$fa-mean(bbwo.dat$fa))/sd(bbwo.dat$fa)
# standardized snag basal area covariate
snag <- (bbwo.dat$snag.ct*0.2295687*10 - mean(bbwo.dat$snag.ct*0.2295687*10))/sd(bbwo.dat$snag.ct*0.2295687*10)  # multipliers convert stem counts ft^2/acre x10 (where 10 is BAF) to m^2/ha
# standardized canopy change covariate
cc <- (bbwo.dat$cc - mean(bbwo.dat$cc))/sd(bbwo.dat$cc)
# standardized latitude covariate
lat <- (bbwo.dat$lat - mean(bbwo.dat$lat))/sd(bbwo.dat$lat)
# standardized elevation residuals covariate
elev.res <- (bbwo.dat$elev.res - mean(bbwo.dat$elev.res))/sd(bbwo.dat$elev.res)

K = 5 # no. survey intervals
npoint = 899 # no. survey points

#create vector with effort indicator - first interval = 3 minutes (= 1), remaining intervals = 2 minutes (= 0)
ef <- rep(NA, K)
ef[1] <- 1
ef[2:5] <- 0

#create vector with survey type indicator(passive = 0; broadcast = 1)
itype <- rep(NA, K)
itype[1:2] <- 0
itype[3:5] <- 1

# create matrix with encounter histories
y = as.matrix(bbwo.dat[,3:7], ncol = 5, nrow = 899)

#write winbugs model to file ###################################################
modelFilename = 'bbwoEAmodel.txt'
cat('
model {

#prior for random fire-area effect
for (j in 1:nfire){
  fire[j] ~ dnorm(0, fire.tau)
  }

# prior for random fire effect variance
fire.tau ~ dgamma(0.001, 0.001)
fire.sd <- pow(fire.tau, -0.5)

# prior for occupancy intercept
psi0 ~ dunif(0,1)
b0 <- log(psi0/(1-psi0))

#priors for logit-linear model coefficients
a1 ~ dunif(-10,10)
a2 ~ dunif(-10,10)
b1 ~ dnorm(0,0.001)
b2 ~ dnorm(0,0.001)
b3 ~ dnorm(0,0.001)
b4 ~ dnorm(0,0.001)
b5 ~ dnorm(0,0.001)

# prior for p intercept
p0 ~ dunif(0,1)
a0 <- log(p0/(1-p0))

for (k in 1:K){
  logit(p[k]) <- a0 + a1*ef[k] + a2*itype[k]       # observation model
  }
for (i in 1:npoint){
  for (k in 1:K){
    muy[i,k] <- z[i]*p[k]
    y[i,k] ~ dbern(muy[i,k])
    }
  z[i] ~ dbern(psi[i])     # state model
  logit(psi[i]) <- b0 + fire[fire.ind[i]] + b1*fa[i] + b2*snag[i] + b3*lat[i] + b4*cc[i] + b5*elev.res[i]
}
propZ <- sum(z[])/npoint
}
', fill=TRUE, file=modelFilename)

# Run WinBUGS ##################################################################

# data needed for WinBUGS
data <- list ("y", "K", "npoint", "fire.ind", "nfire", "itype", "fa", "snag", "lat", "cc", "ef", "elev.res")

# initial values for WinBUGS
zst<-rbinom(npoint,1,0.5)
beta.st <- rnorm(1,0,1)
fire.sd.st <-runif(1,0,5)
pst <- runif(1,.05,.6)
psist <- runif(1, .05, .6)

inits <- function(){
list (z  = zst, p0 = pst, a1 = runif(1, -1, 1), a2 = runif(1, -1, 1), psi0 = psist, fire = rnorm(nfire,0,1), fire.tau = runif(1,.3,.6), b1 = beta.st, b2 = beta.st, b3 = beta.st, b4 = beta.st, b5 = beta.st)
}

#parameters to trace
parameters <- c("a0", "p0", "a1", "a2", "b0", "psi0", "fire", "fire.sd", "b1", "b2", "b3", "b4", "b5", "propZ")

#call WinBUGS, store results in 'out'
out <- bugs(data, inits, parameters, "bbwoEAmodel.txt", n.thin=2, n.chains=2, n.burnin=20000, n.iter=50000, bugs.directory = "c:/Program Files (x86)/WinBUGS14/", debug=F)

#calculate interval specific and type specific p's and overall p
# 3 min interval detection probability
p3min <- inv.logit(out$sims.list$a0 + out$sims.list$a1)
p3min.mn <- mean(p3min)
p3min.025 <- quantile(p3min, prob = 0.025)
p3min.975 <- quantile(p3min, prob = 0.975)

# overall passive detection probability
p.p <- 1 - (1-out$sims.list$p0)*(1-p3min)
p.p.mn <- mean(p.p)
p.p.025 <- quantile(p.p, prob = 0.025)
p.p.975 <- quantile(p.p, prob = 0.975)

# broadcast interval detection probability
p.b <- inv.logit(out$sims.list$a0 + out$sims.list$a2)
p.b.mn <- mean(p.b) 
p.b.025 <- quantile(p.b, prob = 0.025)
p.b.975 <- quantile(p.b, prob = 0.975)

# overall playback detection probability
p.b.tot <- 1 - (1-p.b)^3
p.b.tot.mn <- mean(p.b.tot)
p.b.tot.025 <- quantile(p.b.tot, prob = 0.025)
p.b.tot.975 <- quantile(p.b.tot, prob = 0.975)

# overall detection probability
poverall <- 1 - (1 - p.pa)*(1 - p.br.tot)
poverall.mn <- mean(poverall)
poverall.025 <- quantile(poverall, prob = 0.025)
poverall.975 <- quantile(poverall, prob = 0.975)
