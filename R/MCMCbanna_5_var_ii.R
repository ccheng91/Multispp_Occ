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
warnings()
### convert count into present absent "1 or 0"
n <- which((X[] > 1))
X[n] <- 1            

##### DATA AUGMENTATION ######
# Create all zero encounter histories to add to the detection array X 
# as part of the data augmentation to account for additional 
# species (beyond the n observed species). 

# nzeroes is the number of all-zero encounter histories to be added
nzeroes = 25
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
occ.data <- list(n=n, nzeroes=nzeroes, J=J, K=K, X=Xaug, elev=ele, pop=pop, punish=punish, reach=reach,
                 income=income,camhour = camhour, cam_angle=cam_angle)

X <-Xaug
elev <- ele
dist<- distance
occ.data.qou <- list('n', 'nzeroes', 'J', 'K', 'X', 'elev', 'pop', 'punish', 'reach',
                     'income','camhour' , 'cam_angle')
#Specify the parameters to be monitored      
occ.params <- c('u', 'v', 'w', 
                'mu.u', 'mu.v', 'tau.u', 'tau.v','omega',
                'mu.a1','mu.a3', 'mu.a5', 'mu.a6', 'mu.a11',
                'mu.b1', 'mu.b2',
                'a1', 'a3', 'a5', 'a6', 'a11', 
                'b1','b2',
                'Z', 'Nsite',
                'p.fit','p.fitnew','N')         

modelFile='R/model_5_var_ii.txt'

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
    a3=rnorm(n+nzeroes), 
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

fit <- jags(occ.data.qou, occ.inits, occ.params, modelFile,     # test run
            n.chains=3, n.iter=50, n.burnin=30, n.thin=1)

#fit <- jags(occ.data, occ.inits, occ.params, modelFile,     
#            n.chains=3, n.iter=100000, n.burnin=70000, n.thin=20)

## run model parallely, each chain load to different core 'n.cluster = NO. of core' you using don't use all core of you computer
fitppp <- jags.parallel(occ.data.qou, occ.inits, occ.params, modelFile,     ## 10 hours
                        n.chains=3, n.iter=100000, n.burnin=30000, n.thin=20, n.cluster = 3, export_obj_names =c('Z.int') )

summary(fitppp)
summary(fit$BUGSoutput$sims.array)
#fit.mcmc <-as.mcmc(fit)
#summary(fit.mcmc)
traceplot(fit)
#fit <- jags(occ.data, occ.inits, occ.params, modelFile,     
#            n.chains=3, n.iter=50, n.burnin=30, n.thin=1)

#fit <- jags(occ.data, occ.inits, occ.params, modelFile,     ## 3 hours run
#            n.chains=3, n.iter=8000, n.burnin=3000, n.thin=8)

#fit <- jags(occ.data, occ.inits, occ.params, modelFile,     ##  final quality 
#            n.chains=3, n.iter=50000, n.burnin=5000, n.thin=2)

############################################################################################################################################
### analysis
fit <- fitppp
library("mcmcplots")
library("ggplot2")
mcmcplot(fit)

par <- c("mu.a1","mu.a3","mu.a5","mu.a6","mu.a11","mu.b1","mu.b2")
mcmcplot(fit, parms = par)


#mcmcfit <- as.mcmc.rjags(fit)
#mu.fit <- coda.samples(mcmcfit, c("a1"),thin = 1,n.iter = 1000)
#heidel.diag(mcmcfit)
#raftery.diag(mcmcfit)
#gelman.diag(mcmcfit)
#heidel.diag(mcmcfit)
#summary(mcmcfit)
#n <- fit$BUGSoutput$sims.list$N
#mean(n)
#quantile(n, probs = 0.025)
#quantile(n, probs = 0.975)

## use Gelman and Rubin's diagnosis diagnose convergence
galman.diag <- gelman.diag(mcmcfit, multivariate=FALSE)
galman.diag$psrf[354]
a <- as.data.frame(galman.diag$psrf)
b<-which(a[,1]>1.1)
a[b,1]
rname <- rownames(a)
rname[b]

mu.a1 <- fit$BUGSoutput$sims.list$mu.a1
mu.a3 <- fit$BUGSoutput$sims.list$mu.a3
mu.a5 <- fit$BUGSoutput$sims.list$mu.a5
mu.a6 <- fit$BUGSoutput$sims.list$mu.a6
mu.a11 <- fit$BUGSoutput$sims.list$mu.a11
mu.b1 <- fit$BUGSoutput$sims.list$mu.b1
mu.b2 <- fit$BUGSoutput$sims.list$mu.b2
Nsite <- fit$BUGSoutput$sims.list$Nsite
u <- fit$BUGSoutput$sims.list$u
v <- fit$BUGSoutput$sims.list$v

mu.beta <- matrix(000,7,5)
mu.beta[1,3]<-mean(mu.a1)
mu.beta[1,4]<-quantile(mu.a1, prob=0.025)
mu.beta[1,5]<-quantile(mu.a1, prob=0.975)

mu.beta[2,3]<-mean(mu.a3)
mu.beta[2,4]<-quantile(mu.a3, prob=0.025)
mu.beta[2,5]<-quantile(mu.a3, prob=0.975)

mu.beta[3,3]<-mean(mu.a5)
mu.beta[3,4]<-quantile(mu.a5, prob=0.025)
mu.beta[3,5]<-quantile(mu.a5, prob=0.975)

mu.beta[4,3]<-mean(mu.a6)
mu.beta[4,4]<-quantile(mu.a6, prob=0.025)
mu.beta[4,5]<-quantile(mu.a6, prob=0.975)

mu.beta[5,3]<-mean(mu.a11)
mu.beta[5,4]<-quantile(mu.a11, prob=0.025)
mu.beta[5,5]<-quantile(mu.a11, prob=0.975)

mu.beta[6,3]<-mean(mu.b1)
mu.beta[6,4]<-quantile(mu.b1, prob=0.025)
mu.beta[6,5]<-quantile(mu.b1, prob=0.975)

mu.beta[7,3]<-mean(mu.b2)
mu.beta[7,4]<-quantile(mu.b2, prob=0.025)
mu.beta[7,5]<-quantile(mu.b2, prob=0.975)

mu.beta <- as.data.frame(mu.beta)
mu.beta[,1] <- par
mu.beta[,2] <- c("elevation","population","punishment","outreach","income","camhour","camangle")
names(mu.beta) <- c("coef","name","mean","lower","upper")
write.csv(mu.beta, file = "result/beta_coef.csv")


#-------------------------- Melt plot of posterior distributions of beta coefficients
library(reshape)
library(denstrip)
library(lattice)

## fast way
caterplot(fit, c("mu.a1","mu.a3","mu.a5","mu.a6","mu.a11","mu.b1","mu.b2"), reorder=FALSE,
          labels=c("elevation","population","punishment","outreach","income","camhour","camangle"))
caterplot(fit, c("a1"), reorder=FALSE)

## fast way

preds<-as.data.frame(cbind(mu.a1,mu.a3,mu.a5,mu.a6,mu.a11,mu.b1,mu.b2))
# sort effects by median size
idx<-sort(abs(apply(preds,2,median)),index.return=T,decreasing=F)$ix
idx2<-c(1,2,3,4,5,6,7)

names(preds) <- c("elevation","population","punishment","outreach","income","camhour","camangle")

# apply and sort labels
labs=c("elevation","population","punishment","outreach","income","camhour","camangle")[idx2]

mp=melt(preds[,idx2])
rpp=bwplot(variable~value,data=mp,xlab=list(label="Standardized beta coefficients",cex=1),
           ,panel = function(x, y) { 
             #grid.segments(1,0,0,0)
             xlist <- split(x, factor(y))
             for (i in seq(along=xlist))
               panel.denstrip(x=xlist[[i]], at=i)
           },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1),y=list(draw=T,labels=labs)))
print(rpp)

mean(preds[,1],na.rm = TRUE)
mean(preds[,2],na.rm = TRUE)
mean(preds[,3],na.rm = TRUE)
mean(preds[,4],na.rm = TRUE)
mean(preds[,5],na.rm = TRUE)
mean(preds[,6],na.rm = TRUE)
mean(preds[,7],na.rm = TRUE)


# draw line at 0 across
trellis.focus("panel", 1, 1)
panel.abline(v=0,col=1,lty=2)
trellis.unfocus()


?panel.bwplot 
######

fit <- fit.par
fit$BUGSoutput$DIC
mu.a1 <- fit$BUGSoutput$sims.list$mu.a1
a1<-fit$BUGSoutput$sims.list$a1
a3<-fit$BUGSoutput$sims.list$a3
a6<-fit$BUGSoutput$sims.list$a6
a5<-fit$BUGSoutput$sims.list$a5
a11<-fit$BUGSoutput$sims.list$a11
u <- fit$BUGSoutput$sims.list$u
v <- fit$BUGSoutput$sims.list$v
mean(fit$BUGSoutput$sims.list$a1)
quantile(fit$BUGSoutput$sims.list$mu.a4, probs = 0.025)
quantile(fit$BUGSoutput$sims.list$mu.a5, probs = 0.975)


## beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices ##
## beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices ##
## beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices ##
#### a1 ###### 
betaa1 <- matrix(rep(0), nrow=length(spp), ncol=4)
for (i in 1:length(spp)) {
  betaa1[i,1] <- mean(a1[,i])
  betaa1[i,2] <- sd(a1[,i])
  betaa1[i,3] <- quantile(a1[,i], prob=0.025)
  betaa1[i,4] <- quantile(a1[,i], prob=0.975)
}

betaa1 <- as.data.frame(betaa1)
rownames(betaa1) <- spp
betaa1[,5] <- spp
colnames(betaa1) <- c("Mean", "SD", "lower", "higher","spp")
sig<- (betaa1$higher - abs(betaa1$lower))

betaa1 <- transform(betaa1, spp = reorder(spp,-Mean))
limit <- aes(ymax = higher, ymin=lower)
P <- ggplot(betaa1, aes(y=Mean, x=spp))
P + geom_point(stat="identity") +
  theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust = 1, size=18), axis.line.y=element_line()) + geom_hline(yintercept=0)+
  theme(text = element_text(size=15)) + geom_errorbar(limit, na.rm =F, width=0.25) + ylab("Beta coefficient of elevation parameters") + labs(x = "") +
  scale_y_continuous(breaks = seq(-3, 3, 0.5))

##### a3 #####
betaa3 <- matrix(rep(0), nrow=length(spp), ncol=4)
for (i in 1:length(spp)) {
  betaa3[i,1] <- mean(a3[,i])
  betaa3[i,2] <- sd(a3[,i])
  betaa3[i,3] <- quantile(a3[,i], prob=0.025)
  betaa3[i,4] <- quantile(a3[,i], prob=0.975)
}

betaa3 <- as.data.frame(betaa3)
rownames(betaa3) <- spp
betaa3[,5] <- spp
colnames(betaa3) <- c("Mean", "SD", "lower", "higher","spp")

betaa3 <- transform(betaa3, spp = reorder(spp, -Mean))
limit <- aes(ymax = higher, ymin=lower)
P <- ggplot(betaa3, aes(y=Mean, x=spp))
P + geom_point(stat="identity") + scale_y_continuous(breaks = seq(-8, 8, 1))+
  theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust = 1, size=18), axis.line.y=element_line()) + geom_hline(yintercept=0)+
  theme(text = element_text(size=15)) + geom_errorbar(limit, na.rm =F, width=0.25) + ylab("Beta coefficient of human density parameters") +labs(x = "") 

##### a5 #####
betaa5 <- matrix(rep(0), nrow=length(spp), ncol=4)
for (i in 1:length(spp)) {
  betaa5[i,1] <- mean(a5[,i])
  betaa5[i,2] <- sd(a5[,i])
  betaa5[i,3] <- quantile(a5[,i], prob=0.025)
  betaa5[i,4] <- quantile(a5[,i], prob=0.975)
}

betaa5 <- as.data.frame(betaa5)
rownames(betaa5) <- spp
betaa5[,5] <- spp
colnames(betaa5) <- c("Mean", "SD", "lower", "higher","spp")

betaa5 <- transform(betaa5, spp = reorder(spp, -lower))
#betaa5$spp <- reorder(betaa5$spp, -(betaa5$lower))
limit <- aes(ymax = higher, ymin=lower)
P <- ggplot(betaa5, aes(y=Mean, x=spp))
P + geom_point(stat="identity")  + 
  theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust = 1, size=18), axis.line.y=element_line()) + geom_hline(yintercept=0)+ 
  theme(text = element_text(size=15)) + geom_errorbar(limit, na.rm =F, width=0.25) + ylab("Beta coefficient of punishment parameters") + labs(x = "") +
  scale_y_continuous(breaks = seq(-8,30, 2))

##### a6 ######
betaa6 <- matrix(rep(0), nrow=length(spp), ncol=4)
for (i in 1:length(spp)) {
  betaa6[i,1] <- mean(a6[,i])
  betaa6[i,2] <- sd(a6[,i])
  betaa6[i,3] <- quantile(a6[,i], prob=0.025)
  betaa6[i,4] <- quantile(a6[,i], prob=0.975)
}

betaa6 <- as.data.frame(betaa6)
rownames(betaa6) <- spp
betaa6[,5] <- spp
colnames(betaa6) <- c("Mean", "SD", "lower", "higher","spp")

betaa6 <- transform(betaa6, spp = reorder(spp, -lower))
limit <- aes(ymax = higher, ymin=lower)
P <- ggplot(betaa6, aes(y=Mean, x=spp))
P + geom_point(stat="identity") + 
  theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust = 1, size=18), axis.line.y=element_line()) + geom_hline(yintercept=0) + 
  theme(text = element_text(size=15)) + geom_errorbar(limit, na.rm =F, width=0.25) + ylab("Beta coefficient of outreach parameters")  + labs(x = "") +
 scale_y_continuous(breaks = seq(-4,8, 1))

#### a11 #####
betaa11 <- matrix(rep(0), nrow=length(spp), ncol=4)
for (i in 1:length(spp)) {
  betaa11[i,1] <- mean(a11[,i])
  betaa11[i,2] <- sd(a11[,i])
  betaa11[i,3] <- quantile(a11[,i], prob=0.025)
  betaa11[i,4] <- quantile(a11[,i], prob=0.975)
}

betaa11 <- as.data.frame(betaa11)
rownames(betaa11) <- spp
betaa11[,5] <- spp
colnames(betaa11) <- c("Mean", "SD", "lower", "higher","spp")

betaa11 <- transform(betaa11, spp = reorder(spp, -lower))
limit <- aes(ymax = higher, ymin=lower)
P <- ggplot(betaa11, aes(y=Mean, x=spp))
P + geom_point(stat="identity")  + 
  theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust = 1, size=18), axis.line.y=element_line()) + geom_hline(yintercept=0)+  
  theme(text = element_text(size=15)) + geom_errorbar(limit, na.rm =F, width=0.25) + ylab("Beta coefficient of income parameters") + labs(x = "") +
  scale_y_continuous(breaks = seq(-1,2, 0.5))

## beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices ##
## beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices ##
## beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices ##
## beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices #### beta cof of each speices ##

quantile(a1, prob=0.025)
spp[5]
spp[1]
spp[2]
fit$BUGSoutput$DIC
         
summary(fit)
print(fit)
p.fit <- fit$BUGSoutput$sims.list$p.fit
p.fitnew <- fit$BUGSoutput$sims.list$p.fitnew
modelfit <- length(which(p.fit-p.fitnew>0))/length(p.fit)

modelfit
n0 <- w[(n+1):(n+nzeroes)]

summary(w)

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
         
### 
occupancy<-as.data.frame(occupancy)
spp <- row.names(occupancy)
occupancy <- cbind(spp,occupancy)
ggplot(occupancy, aes(y=Mean, x = spp)) + geom_bar(stat="identity") 

#####
#reach effect

####
r.dummy=seq(min(sitecov$reach),max(sitecov$reach),length.out=100)
mreach <- mean(sitecov$reach)
sdreach<- sd(sitecov$reach)
hreach1=(r.dummy-mreach)/sdreach

psi=matrix(0,nrow=length(spp),ncol=dim(u)[1])
richness = matrix(0, nrow=length(r.dummy), ncol=dim(u)[1])
for (j in 1:length(r.dummy)) {
  for(i in 1: length(spp)) {
    for (k in 1:dim(u)[1]) {
      #psi[i,k] <- plogis(u[k,i] +mu.a1[k]*mean(ele) + mu.a3[k]*pop.ave + mu.a5[k]*punish.ave + 
 #                          mu.a6[k]*mreach + mu.a11[k]*income.ave)  
            psi[i,k] <- plogis(u[k,i] +a1[k,i]*mean(ele) + a3[k,i]*pop.ave + a5[k,i]*punish.ave + 
                                a6[k,i]*hreach1[j] + a11[k,i]*income.ave  )  
    } 
  }
  richness[j,] <- apply(psi,2,sum)  
}   
richness1<-cbind(apply(richness,1,mean),apply(richness,1,quantile,0.975),apply(richness,1,quantile,0.025))
Reachrichness<-cbind(r.dummy, richness1)
Reachrichness <- as.data.frame(Reachrichness)
names(Reachrichness) <- c("reach", "richness", "up", "low")

plot(Reachrichness$richness ~ Reachrichness$reach )

plot(richness ~ reach, data=Reachrichness, type="l", ylim=c(0,15), 
     ylab="Predicted species richness", xlab= "Frequncy of outreach by park staff(times/year)", col = "blue", lwd = 3)
lines(low ~ reach, Reachrichness, type="l", col=gray(0.5), lty=2)
lines(up ~ reach, Reachrichness,type="l", col=gray(0.5), lty=2)


#punish effect
p.dummy=seq(0.1, 5,length.out=100)

mpunish <- mean(sitecov$punish)
sdpunish<- sd(sitecov$punish)
hpunish1=(p.dummy-mpunish)/sdpunish

psi=matrix(0,nrow=length(spp),ncol=dim(u)[1])
richness.p = matrix(0, nrow=length(hpunish1), ncol=dim(u)[1])
for (j in 1:length(hpunish1)) {
  for(i in 1: length(spp)) {
    for (k in 1:dim(u)[1]) {
      psi[i,k] <- plogis(u[k,i] +a1[k,i]*mean(ele) + a3[k,i]*pop.ave + a5[k,i]*hpunish1[j] + 
                           a6[k,i]*reach.ave + a11[k,i]*income.ave  )  
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

#punish effect when outreach is low 

p.dummy=seq(0.1, 5,length.out=100)
mpunish <- mean(sitecov$punish)
sdpunish<- sd(sitecov$punish)
hpunish1=(p.dummy-mpunish)/sdpunish
minreach <- min(reach)

psi=matrix(0,nrow=length(spp),ncol=dim(u)[1])
richness.p = matrix(0, nrow=length(hpunish1), ncol=dim(u)[1])
for (j in 1:length(hpunish1)) {
  for(i in 1: length(spp)) {
    for (k in 1:dim(u)[1]) {
      psi[i,k] <- plogis(u[k,i] +a1[k,i]*mean(ele) + a3[k,i]*pop.ave + a5[k,i]*hpunish1[j] + 
                           a6[k,i]*minreach + a11[k,i]*income.ave  )  
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

#punish effect when outreach is high 

p.dummy=seq(0.1, 5,length.out=100)
mpunish <- mean(sitecov$punish)
sdpunish<- sd(sitecov$punish)
hpunish1=(p.dummy-mpunish)/sdpunish
maxreach <- max(reach)

psi=matrix(0,nrow=length(spp),ncol=dim(u)[1])
richness.p = matrix(0, nrow=length(hpunish1), ncol=dim(u)[1])
for (j in 1:length(hpunish1)) {
  for(i in 1: length(spp)) {
    for (k in 1:dim(u)[1]) {
      psi[i,k] <- plogis(u[k,i] +mu.a1[k]*mean(ele) + mu.a3[k]*pop.ave + mu.a5[k]*hpunish1[j] + 
                    mu.a6[k]*maxreach + mu.a11[k]*income.ave  )  
#      psi[i,k] <- plogis(u[k,i] +a1[k,i]*mean(ele) + a3[k,i]*pop.ave + a5[k,i]*hpunish1[j] + 
#                           a6[k,i]*maxreach + a11[k,i]*income.ave  )  
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


#reach effect when pun is high

####
r.dummy=seq(0,15,length.out=100)
mreach <- mean(sitecov$reach)
sdreach<- sd(sitecov$reach)
hreach1=(r.dummy-mreach)/sdreach
maxpun <- max(punish)

psi=matrix(0,nrow=length(spp),ncol=dim(u)[1])
richness = matrix(0, nrow=length(r.dummy), ncol=dim(u)[1])
for (j in 1:length(r.dummy)) {
  for(i in 1: length(spp)) {
    for (k in 1:dim(u)[1]) {
      psi[i,k] <- plogis(u[k,i] +a1[k,i]*mean(ele) + a3[k,i]*pop.ave + a5[k,i]*maxpun + 
                           a6[k,i]*hreach1[j] + a11[k,i]*income.ave )  
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

##### elevation effcet #######
####
e.dummy=seq(min(sitecov$ele),max(sitecov$ele),length.out=100)
melev <- mean(sitecov$ele)
sdelev <- sd(sitecov$ele)
helev1=(e.dummy-melev)/sdelev

psi=matrix(0,nrow=length(spp),ncol=dim(u)[1])
richness = matrix(0, nrow=length(e.dummy), ncol=dim(u)[1])
for (j in 1:length(e.dummy)) {
  for(i in 1: length(spp)) {
    for (k in 1:dim(u)[1]) {
      psi[i,k] <- plogis(u[k,i] +a1[k,i]*helev1[j] + a3[k,i]*pop.ave + a5[k,i]*punish.ave + 
                           a6[k,i]*mreach + a11[k,i]*income.ave )  
    } 
  }
  richness[j,] <- apply(psi,2,sum)  
}   
richness1<-cbind(apply(richness,1,mean),apply(richness,1,quantile,0.975),apply(richness,1,quantile,0.025))
Reachrichness<-cbind(e.dummy, richness1)
Reachrichness <- as.data.frame(Reachrichness)
names(Reachrichness) <- c("elevation", "richness", "up", "low")

plot(Reachrichness$richness ~ Reachrichness$elevation)

plot(richness ~ elevation, data=Reachrichness, type="l", ylim=c(0,20), 
     ylab="Predicted species richness", xlab= "Elevation(m)", col = "blue", lwd = 3)
lines(low ~ elevation, Reachrichness, type="l", col=gray(0.5), lty=2)
lines(up ~ elevation, Reachrichness,type="l", col=gray(0.5), lty=2)

library(tidyr)
siterichness.est <- matrix(999, nrow = 115, ncol=3)
siterichness.est[,1] <- colMeans(Nsite)
siterichness.sd <- apply(Nsite,2,sd)
siterichness.est[,2] <- siterichness.est[,1] + siterichness.sd
siterichness.est[,3] <- siterichness.est[,1] - siterichness.sd
siterichness.est<- as.data.frame(siterichness.est)
names(siterichness.est) <- c("mean", "higher","lower")

siterichness.est<-cbind(siterichness.est,sitecov)
P <- ggplot(siterichness.est, aes(y=mean, x=ele))
P + geom_point() + theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.line.y=element_line(), axis.line.x=element_line()) +  
  theme(text = element_text(size=10),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +  geom_errorbar(aes(ymax=higher, ymin=lower), alpha=0.3 , width=10) +
        scale_x_continuous(breaks = seq(600, 1850, by = 50)) + scale_y_continuous(breaks = seq(0,15, by = 1)) + 
        geom_smooth(method='lm') + xlab("Elevation(m)") + ylab("Estimated species richness")  

#####

abline(lm(siterichness.est$mean ~ siterichness.est$ele))
fit <- lm(siterichness.est$siterichness.est ~ siterichness.est$dis)
head(siterichness.est)

mean(siterichness.est$mean)
mean(siterichness.sd)


head(siterichness.est)
siterichness.est<- as.data.frame(siterichness.est)
siterichness.est <- gather(siterichness.est, rep, richness, V2:V12501)
ggplot(siterichness.est, aes(V1,richness))  + geom_point(alpha = 1/10)

