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
nzeroes = 25     #### 25 unobserved spp not 20 
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
occ.data.qou <- list('n', 'nzeroes', 'J', 'K', 'X', 'elev', 'pop', 'punish', 'reach',
                     'income','camhour' , 'cam_angle')
#Specify the parameters to be monitored      
occ.params <- c('u', 'v', 'w', 
                'mu.u', 'mu.v', 'tau.u', 'tau.v','omega',
                'mu.a1','mu.a3', 'mu.a5', 'mu.a6', 'mu.a11',
                'mu.b1', 'mu.b2',
                'a1', 'a3', 'a5', 'a6', 'a11', 
                'b1','b2',
                'Z', 'Nsite','N', ### + N
                'p.fit','p.fitnew')         

modelFile='R/model_5_var_ii.txt'  ### dui de 

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

fitppp <- jags.parallel(occ.data.qou, occ.inits, occ.params, modelFile,     ## 10 hours
                        n.chains=5, n.iter=100000, n.burnin=50000, n.thin=20, n.cluster = 5, export_obj_names =c('Z.int') ) # changed

save.image("~/Desktop/data/Multispecies_occupancy/Workplace/VarII_100000_50000_5_20161123.RData") # change path~~~~~# change path~~~~~
