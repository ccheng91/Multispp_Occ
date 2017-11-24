### MSMSOccbanna
library(R2jags)

### 3-demensinal data
####
rm(list=ls(all=TRUE))
list.filenames <- list.files(path = "data/All_widedata", pattern="*.csv")
list.filenames <- paste(c("data/All_widedata/"), list.filenames, sep="")
# create an empty list that will serve as a container to receive the incoming files
spp.temp <- c("blackbear","brushtailedporcupine", "chineseferretbadger",  "commonmacaque", "commonpalmcivet","crabeatingmongoose", "camhour","dhole",               
              "gaur", "goral","hogbadger","leopardcat","maskedpalmcivet", "muntjac", "pigtailedmacaque", "porcupine", "sambar", "serow", "smallindiancivet", 
              "spotedlinsang", "weasel", "wildboar",   "yellowthroatedmarten")
spp <- spp.temp[-7]
list.data<-list()
# create a loop to read in your data
for (i in 1:length(list.filenames)){
  a <- read.csv(list.filenames[i])
  a <- a[-1]
  a <- as.matrix(a)
  colnames(a) <- NULL
  list.data[[i]]<-a
}

names(list.data)<-spp.temp
camhour <- as.data.frame(list.data[7]) # save camhour data
a <- read.csv(list.filenames[1]) 
station <- a[,1] # GET colnames of camhour 
a <- a[-1] # remove "station"
colnames(camhour) <- colnames(a)
list.data <- list.data[-7] # remove camhour from spp list 

##############################
# use 5-day as a sampling period
sitecov <- read.csv("data/sitecov_temp_Cheng_reorded.csv")
days <- sitecov$days
days.5 <- days%/%5
daysd.5<- days%%5

fac <- rep(1:max(days.5), each = 5)       # the index of every 5 day
fac <- c(fac, rep((max(fac) + 1), max(daysd.5))) # add the last 4 days at the end
rep <- c(1:max(fac)) # number of replication 

list.data.5d <-list() # initiate a empty list to store data
for(i in 1:length(list.data)){
  occ.5d.row <-matrix(999, nrow= nrow(list.data[[1]]), ncol=max(fac))
  for(j in 1:nrow(list.data[[1]])){
    q <- tapply(list.data[[i]][j,], fac, FUN = max, na.rm=T)
    q[which(q == -Inf)] <- NA
    occ.5d.row[j,] <- q
    rownames(occ.5d.row) <- NULL
  }
  list.data.5d[[i]] <- occ.5d.row
}
# use 5-day as a sampling period
#####################################

X <- array(
  c(list.data.5d[[1]],list.data.5d[[2]],list.data.5d[[3]],list.data.5d[[4]],
    list.data.5d[[5]],list.data.5d[[6]],list.data.5d[[7]],list.data.5d[[8]],
    list.data.5d[[9]],list.data.5d[[10]],list.data.5d[[11]],list.data.5d[[12]],
    list.data.5d[[13]],list.data.5d[[14]],list.data.5d[[15]],list.data.5d[[16]],
    list.data.5d[[17]],list.data.5d[[18]],list.data.5d[[19]],list.data.5d[[20]],
    list.data.5d[[21]],list.data.5d[[22]]), dim=c(nrow(list.data[[1]]) , max(fac), length(list.data))
)

dimnames(X) <- list(station, rep, spp)

### convert count into present absent "1 or 0"
n <- which((X[] > 1))
X[n] <- 1            

##### DATA AUGMENTATION ######
# Create all zero encounter histories to add to the detection array X 
# as part of the data augmentation to account for additional 
# species (beyond the n observed species). 

# nzeroes is the number of all-zero encounter histories to be added
nzeroes = 20
X.zero <- list.data.5d[[1]]
X.zero[which(X.zero > 0)] <- 0
Xaug <- array(0, dim=c(dim(X)[1],dim(X)[2],dim(X)[3]+nzeroes))
Xaug[,,(dim(X)[3]+1):dim(Xaug)[3]] = rep(X.zero, nzeroes)
dimnames(X)=NULL # ??? why set dimname of X have to be Null??
Xaug[,,1:dim(X)[3]] <-  X

# K is a vector of length J indicating the number of reps at each point j  
KK <- X.zero
a=which(KK==0); KK[a] <- 1
K=apply(KK,1, sum, na.rm=TRUE)
K=as.vector(K)

################## SAMPLING COVARIATES #######################################
#Standardize Camhour 
camhour<-as.matrix(camhour)
camhour.5d <- matrix(999, nrow= nrow(camhour), ncol=max(fac))
for(i in 1:nrow(camhour)){
  q <- tapply(camhour[i,], fac, FUN = mean, na.rm=T)
  q[which(q == "NaN")] <- NA
  camhour.5d[i,]<- q
  rownames(camhour.5d) <- NULL
}
camhour.5d <- round(camhour.5d, digits=2) # only two decimal places

mcamhour <- mean(camhour.5d, na.rm = T)
sdcamhour <-sd(camhour.5d, na.rm=T)
camhour <- (camhour.5d-mcamhour)/sdcamhour

#cam_angle
cam_angle<-X[,,1] 
NA.n <- which(is.na(cam_angle) == TRUE) 
caman.1 <- which(sitecov$Cam_angle == 1)
cam_angle[caman.1,] <- 1
cam_angle[-caman.1,] <- 0
cam_angle[NA.n] <- NA

################## OCCUPANCY COVARIATES ########################
head(sitecov)
ele.ave <- mean(sitecov$ele)
ele.sd <- sd(sitecov$ele)
ele <- (sitecov$ele-ele.ave)/ele.sd

dis.ave <- mean(sitecov$dis)
dis.sd <- mean(sitecov$dis)
distance <- (sitecov$dis-dis.ave)/dis.sd

pop.ave <- mean(sitecov$pop3000m)
pop.sd <- sd(sitecov$pop3000m)
pop <- (sitecov$pop3000m-pop.ave)/pop.sd

size.ave <- mean(sitecov$size.m)
size.sd <- sd(sitecov$size.m)
pasize <- (sitecov$size.m-size.ave)/size.sd

punish.ave <- mean(sitecov$punish)
punish.sd <- mean(sitecov$punish)
punish <- (sitecov$punish - punish.ave)

reach.ave <- mean(sitecov$reach)
reach.sd <- sd(sitecov$reach)
reach<- (sitecov$reach - reach.ave)/reach.sd

relationship.ave <- mean(sitecov$relationship)
relationship.sd <- sd(sitecov$relationship)
relationship <- (sitecov$relationship - relationship.ave)/relationship.sd

score.P.ave <- mean(sitecov$score.P)
score.P.sd <- sd(sitecov$score.P)
score.P <- (sitecov$score.P - score.P.ave)/score.P.sd

score.H.ave <- mean(sitecov$score.H)
score.H.sd <- sd(sitecov$score.H)
score.H <- (sitecov$score.H - score.H.ave)/score.H.sd

edu.ave <- mean(sitecov$edu)
edu.sd <- sd(sitecov$edu)
edu <- (sitecov$edu - sitecov$edu)/edu.sd

income.ave <- mean(sitecov$income)
income.sd <- sd(sitecov$income)
income <- (sitecov$income - income.ave)/income.sd

######## BAYESIAN MODELING ########
n <- dim(X)[3] ## number of speices 
J <- dim(X)[1] ## numver of sites 

#Specify the data
occ.data <- list(n=n, nzeroes=nzeroes, J=J, K=K, X=Xaug, elev=ele, pop=pop,pasize=pasize, punish=punish, reach=reach,
                 income=income,camhour = camhour, cam_angle=cam_angle)

X <-Xaug
elev <- ele
dist<- distance
occ.data.qou <- list('n', 'nzeroes', 'J', 'K', 'X', 'elev', 'pop','pasize', 'punish', 'reach',
                     'income','camhour' , 'cam_angle')
#Specify the parameters to be monitored      
occ.params <- c('u', 'v', 'w', 
                'mu.u', 'mu.v', 'tau.u', 'tau.v','omega',
                'mu.a1','mu.a3','mu.a4', 'mu.a5', 'mu.a6','mu.a11', 
                'mu.b1', 'mu.b2',
                'a1', 'a3', 'a4', 'a5', 'a6', 'a11', 
                'b1','b2',
                'Z', 'Nsite',
                'p.fit','p.fitnew')         

modelFile='R/model_6_var.txt'

#Specify the initial values

# JAGS to provide initial values for the incompletely 
# observed occupancy state $z$ that are consistent with observed presences
# In other words if $X(j, i, k)=1$, provide an intial value of $1$ for $z(j, i) 
# Unlike WinBUGS and OpenBUGS, if you do not do this, you’ll often (but not always) encounter an error message such as:
# Observed node inconsistent with parents 

Z.int <- array(dim=c( dim(Xaug)[1], dim(Xaug)[3]))
for(i in 1:dim(Xaug)[1]){
  for(j in 1:dim(Xaug)[3]){
    Z.int[i, j] <- max(Xaug[i,,j], na.rm=TRUE)
  }
}

occ.inits = function() {
  omegaGuess = runif(1, n/(n+nzeroes), 1)
  psi.meanGuess = runif(1, .25,1)
  list(
    omega=omegaGuess,w=c(rep(1, n), rbinom(nzeroes, size=1, prob=omegaGuess)),
    u=runif(n+nzeroes,-1,1), v=runif(n+nzeroes,-1,1),
    Z = Z.int,   
    b1=rnorm(n+nzeroes), b2=rnorm(n+nzeroes),
    a1=rnorm(n+nzeroes), 
    a3=rnorm(n+nzeroes), a4=rnorm(n+nzeroes),
    a5=rnorm(n+nzeroes), a6=rnorm(n+nzeroes),
    a11=rnorm(n+nzeroes) 
  )           
}
##

#Z = matrix(rbinom((n+nzeroes)*J, size=1, prob=psi.meanGuess),
#           nrow=J, ncol=(n+nzeroes)) 
#psi.meanGuess = runif(1, .25,1)
## run mutiply chains at same time 
#library(snow)
#library(dclone)
#load.module("lecuyer")
#cl <- makeCluster(3,type = "SOCK")
#parLoadModule(cl, "lecuyer", quiet=T)
#re.int <- parallel.inits(occ.inits, n.chains = 3)
#parJagsModel(cl, name = "rep.fit", file = modelFile, data=occ.data, n.chains = 3,n.adapt = 1000)
#stopCluster(cl)
#parCodaSamples()
### Run the model and call the results ?fit?

fitppp <- jags.parallel(occ.data.qou, occ.inits, occ.params, modelFile,     ## 
                        n.chains=3, n.iter=70000, n.burnin=30000, n.thin=20, n.cluster = 3, export_obj_names =c('Z.int') )

summary(fitppp)
summary(fit$BUGSoutput$sims.array)
fit.mcmc <-as.mcmc(fit)
summary(fit.mcmc)
traceplot(fit)
fit <- jags(occ.data, occ.inits, occ.params, modelFile,     
            n.chains=3, n.iter=50, n.burnin=30, n.thin=1)

#fit <- jags(occ.data, occ.inits, occ.params, modelFile,     ## 3 hours run
#            n.chains=3, n.iter=8000, n.burnin=3000, n.thin=8)

#fit <- jags(occ.data, occ.inits, occ.params, modelFile,     ##  final qulity maybe
#            n.chains=3, n.iter=50000, n.burnin=5000, n.thin=2)
recompile(fit)
load.module("dic")
fit.up <- autojags(fit)

##################################################################################################################################################################
### analysis
library("mcmcplots")
mcmcplot(fit)
mcmcfit <- as.mcmc.rjags(fit)
mu.fit <- coda.samples(mcmcfit, c("a1"),thin = 1,n.iter = 1000)
heidel.diag(mcmcfit)
raftery.diag(mcmcfit)
gelman.diag(fit)
heidel.diag(mcmcfit)
summary(mcmcfit)
fit <- fit.par
fit$BUGSoutput$DIC
mu.a1 <- fit$BUGSoutput$sims.list$mu.a1
a1<-fit$BUGSoutput$sims.list$a1
a2<-fit$BUGSoutput$sims.list$a2
mean(fit$BUGSoutput$sims.list$a1)
quantile(fit$BUGSoutput$sims.list$mu.a4, probs = 0.025)
quantile(fit$BUGSoutput$sims.list$mu.a4, probs = 0.975)

fit$BUGSoutput$mean


quantile(a1, prob=0.025
         fit$BUGSoutput$DIC
         
         summary(fit)
         print(fit)
         p.fit <- fit$BUGSoutput$sims.list$p.fit
         p.fitnew <- fit$BUGSoutput$sims.list$p.fitnew
         modelfit <- length(which(p.fit-p.fitnew>0))/length(p.fit)
         
         fit$BUGSoutput$summary
         
         modelfit
         
         # Matrix of prop model runs > or < 0 for each species for the specified parameter (e.g. a3)
         psim =fit$BUGSoutput$sims.list$a3  #matrix with 1col/spp and 1row/run (n.burnin)
         propruns <-matrix(999,30,2,byrow=T)
         for (i in 1:(nrow(propruns))){
           propruns[i,1]=(length(which(psim[,i]<0)))/nrow(psim)
           propruns[i,2]=(length(which(psim[,i]>0)))/nrow(psim)
         }
         propruns
         
         u = fit$BUGSoutput$sims.list$u
         occ = plogis(u)
         occupancy <- matrix(rep(0), nrow=length(spp), ncol=4)
         for (i in 1:length(spp)) {
           occupancy[i,1] <- mean(occ[,i])
           occupancy[i,2] <- sd(occ[,i])
           occupancy[i,3] <- quantile(occ[,i], prob=0.025)
           occupancy[i,4] <- quantile(occ[,i], prob=0.975)
         }
         rownames(occupancy) <- spp
         colnames(occupancy) <- c("Mean", "SD", "95 lower", "95 higher")
         
         # Detection 
         v = fit$BUGSoutput$sims.list$v
         det = plogis(fit$BUGSoutput$sims.list$v)
         detection <- matrix(rep(0), nrow=length(spp), ncol=4)
         for (i in 1:length(spp)) {
           detection[i,1] <- mean(det[,i])
           detection[i,2] <- sd(det[,i])
           detection[i,3] <- quantile(det[,i], prob=0.025)
           detection[i,4] <- quantile(det[,i], prob=0.975)
         }
         rownames(detection) <- spp
         colnames(detection) <- c("Mean", "SD", "95 lower", "95 higher")
         
         # Posterior distribution of occupancy for each species in each habitat
         u  = fit$BUGSoutput$sims.list$u
         a1 = fit$BUGSoutput$sims.list$a1  #elev
         a2 = fit$BUGSoutput$sims.list$a2  #pop
         a3 = fit$BUGSoutput$sims.list$a3  #punish
         a4 = fit$BUGSoutput$sims.list$a4  #reach
         for (i in 1:length(spp)) {
           par(mfrow=c(2,4))
           hist(a1[,i], xlab = spp[i], main= "Elevation")
           hist(a2[,i], xlab = spp[i], main="distance")
           hist(a3[,i], xlab = spp[i], main="Pop")
           hist(a4[,i], xlab = spp[i], main="PAsize")
           browser()
         }
         Q # quit browser
         
         species.occ = fit$BUGSoutput$sims.list$u
         species.det = fit$BUGSoutput$sims.list$v
         
         
         
         ##################################################################################################################################################################
         library("ggplot2")
         #-------------------------- Line graph of hunting effects in each habitat
         r.dummy=seq(0,10,length.out=100)
         mreach <- mean(sitecov$reach)
         sdreach<- sd(sitecov$reach)
         hreach1=(r.dummy-mreach)/sdreach
         
         #reach effect
         psi=matrix(0,nrow=length(spp),ncol=dim(u)[1])
         richness = matrix(0, nrow=length(hreach1), ncol=dim(u)[1])
         for (j in 1:length(r.dummy)) {
           for(i in 1: length(spp)) {
             for (k in 1:dim(u)[1]) {
               psi[i,k] <- plogis(u[k,i] +a1[k,i]*ele.ave +a2[k,i]*pop.ave + a3[k,i]*punish.ave + a4[k,i]*hreach1[j] )  
             } 
           }
           richness[j,] <- apply(psi,2,sum)  
         }   
         richness1<-cbind(apply(richness,1,mean),apply(richness,1,quantile,0.975),apply(richness,1,quantile,0.025))
         Reachrichness<-cbind(r.dummy, richness1)
         Reachrichness <- as.data.frame(Reachrichness)
         names(Reachrichness) <- c("reach", "richness", "up", "low")
         
         plot(Reachrichness$richness ~ Reachrichness$reach )
         
         plot(richness ~ reach, data=Reachrichness, type="l", ylim=c(0,20), 
              ylab="Predicted species richness", xlab= "Frequncy of outreach by park staff(times/year)", col = "blue", lwd = 3)
         lines(low ~ reach, Reachrichness, type="l", col=gray(0.5), lty=2)
         lines(up ~ reach, Reachrichness,type="l", col=gray(0.5), lty=2)
         
         #punish effect
         
         p.dummy=seq(0.1, 1.5,length.out=50)
         
         mpunish <- mean(sitecov$punish)
         sdpunish<- sd(sitecov$punish)
         hpunish1=(p.dummy-mpunish)/sdpunish
         
         psi=matrix(0,nrow=length(spp),ncol=dim(u)[1])
         richness.p = matrix(0, nrow=length(hpunish1), ncol=dim(u)[1])
         for (j in 1:length(hpunish1)) {
           for(i in 1: length(spp)) {
             for (k in 1:dim(u)[1]) {
               psi[i,k] <- plogis(u[k,i] +a1[k,i]*mean(ele) +a2[k,i]*mean(pop) + a3[k,i]*p.dummy[j] + a4[k,i]*mreach )  
             } 
           }
           richness.p[j,] <- apply(psi,2,sum)  
         }   
         richness.p1<-cbind(apply(richness.p,1,mean),apply(richness.p,1,quantile,0.975),apply(richness.p,1,quantile,0.025))
         punrichness<-cbind(p.dummy, richness.p1)
         punrichness <- as.data.frame(punrichness)
         
         names(punrichness) <- c("punish", "richness", "up", "low")
         
         plot(richness ~ punish, data=punrichness, type="l", ylim=c(0,10), 
              ylab="Predicted species richness", xlab= "Punishment park staff meted out for law violator(times/year)", col = "blue", lwd = 3)
         lines(low ~ punish, punrichness, type="l", col=gray(0.5), lty=2)
         lines(up ~ punish, punrichness,type="l", col=gray(0.5), lty=2)
         
         #elevation
         e.dummy=seq(500, 2000,length.out=50)
         hele1=(e.dummy-ele.ave)/ele.sd
         
         psi=matrix(0,nrow=length(spp),ncol=dim(u)[1])
         richness.e = matrix(0, nrow=length(hele1), ncol=dim(u)[1])
         for (j in 1:length(hele1)) {
           for(i in 1: length(spp)) {
             for (k in 1:dim(u)[1]) {
               psi[i,k] <- plogis(u[k,i] +a1[k,i]*hele1[j] +a2[k,i]*mean(pop) + a3[k,i]*mean(punish) + a4[k,i]*mreach )  
             } 
           }
           richness.e[j,] <- apply(psi,2,sum)  
         }   
         richness.e1<-cbind(apply(richness.e,1,mean),apply(richness.e,1,quantile,0.975),apply(richness.e,1,quantile,0.025))
         elerichness<-cbind(e.dummy, richness.e1)
         elerichness <- as.data.frame(elerichness)
         names(elerichness) <- c("elevation", "richness", "up", "low")
         plot(richness ~ elevation, data=elerichness, type="l", ylim=c(5,20), 
              ylab="Predicted species richness", xlab= "Elevation(m)", col = "blue", lwd = 3)
         lines(low ~ richness, elerichness, type="l", col=gray(0.5), lty=2)
         lines(up ~ richness, elerichness,type="l", col=gray(0.5), lty=2)
         
         ## population 
         pop.dummy <- seq(min(sitecov$pop3000m), max(sitecov$pop3000m),length.out=50)
         hpop1 <- (pop.dummy-pop.ave)/pop.sd
         
         psi=matrix(0,nrow=length(spp),ncol=dim(u)[1])
         richness.pop = matrix(0, nrow=length(hpop1), ncol=dim(u)[1])
         for (j in 1:length(hpop1)) {
           for(i in 1: length(spp)) {
             for (k in 1:dim(u)[1]) {
               psi[i,k] <- plogis(u[k,i] +a1[k,i]*mean(ele) + a2[k,i]*hpop1[j] + a3[k,i]*mean(punish) + a4[k,i]*mreach )  
             } 
           }
           richness.pop[j,] <- apply(psi,2,sum)  
         }   
         richness.pop1<-cbind(apply(richness.pop,1,mean),apply(richness.pop,1,quantile,0.975),apply(richness.pop,1,quantile,0.025))
         poprichness<-cbind(pop.dummy, richness.pop1)
         poprichness <- as.data.frame(poprichness)
         names(poprichness) <- c("population", "richness", "up", "low")
         plot(richness ~ population, data=poprichness, type="l", ylim=c(0,10), 
              ylab="Predicted species richness", xlab= "Population index", col = "blue", lwd = 3)
         lines(low ~ population, poprichness, type="l", col=gray(0.5), lty=2)
         lines(up ~ population, poprichness,type="l", col=gray(0.5), lty=2)
         
         
         #-------------------------- Melt plot of posterior distributions of beta coefficients
         library(ggplot2)
         library(denstrip)
         library(lattice)
         library(reshape)
         
         mu.a1=fit$BUGSoutput$sims.list$mu.a1
         mu.a2=fit$BUGSoutput$sims.list$mu.a2
         mu.a3=fit$BUGSoutput$sims.list$mu.a3
         mu.a4=fit$BUGSoutput$sims.list$mu.a4
         
         toExport<-cbind(mu.a1,mu.a2,mu.a3,mu.a4)
         preds<-as.data.frame(toExport)
         
         idx<-sort(abs(apply(preds,2,median)),index.return=T,decreasing=F)
         idx2<-c(1,2,3,4)
         # apply and sort labels
         labs=c('mu.a1','mu.a2','mu.a3','mu.a4')[idx2]
         mp=melt(preds[,idx2])
         
         rpp=bwplot(variable~value,data=mp,xlab=list(label="relative effect",cex=1),
                    ,panel = function(x, y) { 
                      #grid.segments(1,0,0,0)
                      xlist <- split(x, factor(y))
                      for (i in seq(along=xlist))
                        panel.denstrip(x=xlist[[i]], at=i)
                    },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1),y=list(draw=T,labels=labs)))
         print(rpp)
         
         # draw line at 0 across
         trellis.focus("panel", 1, 1)
         panel.abline(v=0,col=1,lty=2)
         trellis.unfocus()
         dic.samples(fit,1000,"pD")
         
         
         