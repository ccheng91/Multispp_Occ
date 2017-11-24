#The model code below is written for program R and uses the R2WinBUGS package 
#to run WinBUGS as well as the reshape package to format the occurrence data.

#It is designed to estimate static species-specific occupancy and detection, 
#constant across sampling locations using the community model. The model also
#estimates the total species richness N, using data augmentation.
#The data are found in the file "occ data.csv". 


##################### GETTING THE DATA SORTED ##############################
rm(list=ls(all=TRUE))
# Read in the occurence data
data1<-read.csv('data/Borneo_data sample.csv',header=TRUE,sep=",")
data1$Occ <- rep(1, dim(data1)[1])

# See the first ten lines of data
# data1[1:10,]

# How many sightings for each species
total.count=tapply(data1$Occ, data1$species, sum)

# Find the number of unique species
uspecies1=sort(as.character(unique(data1$species)))
# length(uspecies1)

# Find the number of unique sampling locations
upoints=sort(as.character(unique(data1$station)))
J=length(upoints)

# Find the number of unique reps
reps=sort(unique(data1$rep))

# Reshape the data using the R package "reshape"
library(reshape)

# The detection/non-detection data is reshaped into a three dimensional 
# array X where the first dimension, j, is the point; the second 
# dimension, k, is the rep; and the last dimension, i, is the species. 
junk.melt=melt(data1,id.var=c("species", "station", "rep"), measure.var="Occ")
X=cast(junk.melt, station ~ rep ~ species)
dimnames(X)=list(upoints,reps,uspecies1)

# Add in the missing lines (dates where sampling did not occur at each site) with NAs
# In this case there are 134 stations and a max of 30 sampling reps
for (i in 1: dim(X)[3]) {
   b = which(X[,,i] > 0) 
   X[,,i][b] = 1  
   X[,,i][-b] = 0  
   X[1:2,20:30,]=NA; X[3:10,19:30,]=NA; X[11,7:30,]=NA; X[12,18:30,]=NA;
   X[13:14,19:30,]=NA; X[15,16:30,]=NA; X[16:18,19:30,]=NA; X[20,27:30,]=NA;
   X[23,2:30,]=NA; X[26,19:30,]=NA; X[29,25:30,]=NA;
   X[31:34,27:30,]=NA; X[43,22:30,]=NA; X[45,3:30,]=NA; X[46,21:30,]=NA;
   X[49,24:30,]=NA; X[50,15:30,]=NA; X[52,23:30,]=NA; X[53,17:30,]=NA;
   X[54,19:30,]=NA; X[56,28:30,]=NA; X[58,25:30,]=NA; X[59,24:30,]=NA; 
   X[61,16:30,]=NA; X[63,21:30,]=NA; X[66,29:30,]=NA; X[79,22:30,]=NA;
   X[80,8:30,]=NA; X[81:84,22:30,]=NA; X[85,20:30,]=NA; X[86,23:30,]=NA;
   X[87,22:30,]=NA; X[88,6:30,]=NA; X[89,22:30,]=NA; X[90,8:30,]=NA; 
   X[91,22:30,]=NA; X[92,17:30,]=NA; X[93,23:30,]=NA; X[94,13:30,]=NA; 
   X[95,18:30,]=NA; X[96,3:30,]=NA; X[97,5:30,]=NA; X[98,13:30,]=NA; 
   X[99,21:30,]=NA; X[100,26:30,]=NA; X[101:102,30,]=NA; X[103,28:30,]=NA; 
   X[104,30,]=NA; X[105,28:30,]=NA; X[106:107,30,]=NA; X[108,4:30,]=NA; 
   X[109,14:30,]=NA; X[110,22:30,]=NA; X[111,29:30,]=NA; X[112,6:30,]=NA; 
   X[113,29:30,]=NA; X[114,16:30,]=NA; X[115,7:30,]=NA; X[116,3:30,]=NA; 
   X[117:118,4:30,]=NA; X[119,21:30,]=NA; X[120,12:30,]=NA; X[121,21:30,]=NA; 
   X[122,8:30,]=NA; X[123:125,21:30,]=NA; X[126,17:30,]=NA; X[127,21:30,]=NA; 
   X[128,9:30,]=NA; X[129,20:30,]=NA; X[130,9:30,]=NA; X[131:132,18:30,]=NA; 
   X[133,2:30,]=NA; X[134,9:30,]=NA   
}

#Remove the NULL species
X=X[,,-24]
uspecies=dimnames(X)[[3]]
n=length(uspecies)


################## SAMPLING COVARIATES #######################################
#Day and camera hours covars
day=X[,,1]
hours=X[,,1]

for (i in 1:nrow(day)){
  for (j in 1:ncol(day)){ 
    a=which(data1$station==rownames(day)[i] & data1$rep==colnames(day)[j])
    if (length(a)>0) {
      day[i,j]=unique(data1$studyday[a])
      hours[i,j]=unique(data1$camhours[a])
    } }
}

#Standardize day 
mday=mean(as.vector(day),na.rm=T)
sdday=sd(as.vector(day),na.rm=T)
day1<-(day-mday)/sdday 

#Standardize camera hours 
mhours=mean(as.vector(hours),na.rm=T)
sdhours=sd(as.vector(hours),na.rm=T)
hours1<-(hours-mhours)/sdhours


########################## OCCUPANCY COVARIATES ########################
covars<-read.csv('covar_data.csv',header=TRUE,sep=",")
a=pmatch(dimnames(X)[[1]],covars$station)
#covars$station[a]

#Elevation
elev=covars$elev[a]
melev=mean(elev,na.rm=T)
sdelev=sd(elev,na.rm=T)
elev1<-(elev-melev)/sdelev 

#Logging
logged=covars$logged[a]
unlogged=covars$unlogged[a]
newlog=covars$newlog[a]
oldlog=covars$oldlog[a]
invage=covars$invage[a]
minvage=mean(invage,na.rm=T)
sdinvage=sd(invage,na.rm=T)
invage1<-(invage-minvage)/sdinvage

#Hunting
hunt=covars$hunt[a]
#hunt=covars$dist[a]
mhunt=mean(hunt,na.rm=T)
sdhunt=sd(hunt,na.rm=T)
hunt1<-(hunt-mhunt)/sdhunt

#Site
site=covars$site1[a]


########################### DATA AUGMENTATION #######################################
# Create all zero encounter histories to add to the detection array X 
# as part of the data augmentation to account for additional 
# species (beyond the n observed species). 

# nzeroes is the number of all-zero encounter histories to be added
nzeroes = 20

# X.zero is a matrix of zeroes, including the NAs for when a point has not been sampled  
X.zero = matrix(0, nrow=J, ncol=max(data1$rep))
X.zero[1:2,20:30]=NA; X.zero[3:10,19:30]=NA; X.zero[11,7:30]=NA; X.zero[12,18:30]=NA;
X.zero[13:14,19:30]=NA; X.zero[15,16:30]=NA; X.zero[16:18,19:30]=NA; X.zero[20,27:30]=NA;
X.zero[23,2:30]=NA; X.zero[26,19:30]=NA; X.zero[29,25:30]=NA;
X.zero[31:34,27:30]=NA; X.zero[43,22:30]=NA; 
X.zero[45,3:30]=NA; X.zero[46,21:30]=NA;
X.zero[49,24:30]=NA; X.zero[50,15:30]=NA; X.zero[52,23:30]=NA; X.zero[53,17:30]=NA;
X.zero[54,19:30]=NA; X.zero[56,28:30]=NA; X.zero[58,25:30]=NA; X.zero[59,24:30]=NA; 
X.zero[61,16:30]=NA; X.zero[63,21:30]=NA; X.zero[66,29:30]=NA; X.zero[79,22:30]=NA;
X.zero[80,8:30]=NA; X.zero[81:84,22:30]=NA; X.zero[85,20:30]=NA; X.zero[86,23:30]=NA;
X.zero[87,22:30]=NA; X.zero[88,6:30]=NA; X.zero[89,22:30]=NA; X.zero[90,8:30]=NA; 
X.zero[91,22:30]=NA; X.zero[92,17:30]=NA; X.zero[93,23:30]=NA; X.zero[94,13:30]=NA; 
X.zero[95,18:30]=NA; X.zero[96,3:30]=NA; X.zero[97,5:30]=NA; X.zero[98,13:30]=NA; 
X.zero[99,21:30]=NA; X.zero[100,26:30]=NA; X.zero[101:102,30]=NA; X.zero[103,28:30]=NA; 
X.zero[104,30]=NA; X.zero[105,28:30]=NA; X.zero[106:107,30]=NA; X.zero[108,4:30]=NA; 
X.zero[109,14:30]=NA; X.zero[110,22:30]=NA; X.zero[111,29:30]=NA; X.zero[112,6:30]=NA; 
X.zero[113,29:30]=NA; X.zero[114,16:30]=NA; X.zero[115,7:30]=NA; X.zero[116,3:30]=NA; 
X.zero[117:118,4:30]=NA; X.zero[119,21:30]=NA; X.zero[120,12:30]=NA; X.zero[121,21:30]=NA; 
X.zero[122,8:30]=NA; X.zero[123:125,21:30]=NA; X.zero[126,17:30]=NA; X.zero[127,21:30]=NA; 
X.zero[128,9:30]=NA; X.zero[129,20:30]=NA; X.zero[130,9:30]=NA; X.zero[131:132,18:30]=NA; 
X.zero[133,2:30]=NA; X.zero[134,9:30]=NA; 

# Xaug is the augmented version of X.  The first n species were actually observed
# and the n+1 through n+nzeroes species are all-zero encounter histories  
Xaug <- array(0, dim=c(dim(X)[1],dim(X)[2],dim(X)[3]+nzeroes))
Xaug[,,(dim(X)[3]+1):dim(Xaug)[3]] = rep(X.zero, nzeroes)
dimnames(X)=NULL
Xaug[,,1:dim(X)[3]] <-  X

# K is a vector of length J indicating the number of reps at each point j  
KK <- X.zero
a=which(KK==0); KK[a] <- 1
K=apply(KK,1,sum, na.rm=TRUE)
K=as.vector(K)


################################# BAYESIAN MODELING #############################
library(R2WinBUGS)

nsite=max(site)

# Create the necessary arguments to run the bugs() command 
sp.data = list(n=n, nzeroes=nzeroes, J=J, K=K, X=Xaug, day=day1, hours=hours1, 
	hunt=hunt1, newlog=newlog, oldlog=oldlog, elev=elev1,site=site, nsite=max(site)) #,
               # ngroup=max(group1), groups=group2) 

# Specify the parameters to be monitored
sp.params = list('u', 'v', 'mu.u', 'mu.v', 'tau.u', 'tau.v', 'omega', 
                 'N', 
                 'mu.b1', 'mu.b2',
                 'mu.a1', 'mu.a2', 'mu.a3', 'mu.a44', 'mu.a55', 
                 #'mu.a6', 'mu.a7', 
		     'sigma.a8',
                 'b1', 'b2', 'rho1', 'rho2',
                 'a1', 'a2', 'a3', 'a4', 'a5', #'a6', 'a7',
                 'Nsite', 'p.fit', 'p.fitnew')

# Specify the initial values
sp.inits = function() {
    omegaGuess = runif(1, n/(n+nzeroes), 1)
    psi.meanGuess = runif(1, .25,1)
    list(omega=omegaGuess,w=c(rep(1, n), rbinom(nzeroes, size=1, prob=omegaGuess)),
               u=runif(n+nzeroes,-1,1), v=runif(n+nzeroes,-1,1),
               Z = matrix(rbinom((n+nzeroes)*J, size=1, prob=psi.meanGuess),
		   nrow=J, ncol=(n+nzeroes)),   
               b1=rnorm(n+nzeroes), b2=rnorm(n+nzeroes),
		   a1=rnorm(n+nzeroes), a2=rnorm(n+nzeroes),
		   a3=rnorm(n+nzeroes), a4=rnorm(n+nzeroes), 
		   a5=rnorm(n+nzeroes), #a6=rnorm(n+nzeroes), a7=rnorm(n+nzeroes),
       	   a8=matrix(rnorm(nsite*(n+nzeroes)), ncol=nsite)
     )           
}

# Run the model and call the results ?fit?
fit = bugs(sp.data,sp.inits,sp.params,"model_130705_correlation.txt",bugs.directory="c:/MyPrograms/WinBUGS14/",       
		debug=TRUE, n.chains=2, n.iter=100, n.burnin=50, n.thin=2) # quick test
		#debug=TRUE, n.chains=3, n.iter=8000, n.burnin=3000, n.thin=8) # overnight run
	 	debug=TRUE, n.chains=3, n.iter=70000, n.burnin=30000, n.thin=20)	# publication quality; ~3 days run time


######################## VIEWING RESULTS ##############################################
# See a summary of the parameter estimates
fit$summary

###Check model fit
p.fit = fit$sims.list$p.fit
p.fitnew =fit$sims.list$p.fitnew
model.fit= length(which(p.fit-p.fitnew>0))/length(p.fit)
#If model.fit between 0.05 and 0.95, then model fit is good
model.fit

# Matrix of prop model runs > or < 0 for each species for the specified parameter (e.g. a3)
psim =fit$sims.list$a5  #matrix with 1col/spp and 1row/run (n.burnin)
propruns <-matrix(999,30,2,byrow=T)
for (i in 1:(nrow(propruns))){
	propruns[i,1]=(length(which(psim[,i]<0)))/nrow(psim)
	propruns[i,2]=(length(which(psim[,i]>0)))/nrow(psim)
}
propruns

# Table with mean, sd, and CI for occupancy for each observed species
u = fit$sims.list$u
occ = plogis(u)
occupancy <- matrix(rep(0), nrow=length(uspecies), ncol=4)
for (i in 1:length(uspecies)) {
  occupancy[i,1] <- mean(occ[,i])
  occupancy[i,2] <- sd(occ[,i])
  occupancy[i,3] <- quantile(occ[,i], prob=0.025)
  occupancy[i,4] <- quantile(occ[,i], prob=0.975)
  }
rownames(occupancy) <- uspecies
colnames(occupancy) <- c("Mean", "SD", "95 lower", "95 higher")

# Detection 
v = fit$sims.list$v
det = plogis(fit$sims.list$v)
detection <- matrix(rep(0), nrow=length(uspecies), ncol=4)
for (i in 1:length(uspecies)) {
  detection[i,1] <- mean(det[,i])
  detection[i,2] <- sd(det[,i])
  detection[i,3] <- quantile(det[,i], prob=0.025)
  detection[i,4] <- quantile(det[,i], prob=0.975)
}
rownames(detection) <- uspecies
colnames(detection) <- c("Mean", "SD", "95 lower", "95 higher")

# Posterior distribution of occupancy for each species in each habitat
a3 = fit$sims.list$a3  #hunting
a4 = fit$sims.list$a4  #newlog
a5 = fit$sims.list$a5  #oldlog
for (i in 1:length(uspecies)) {
  par(mfrow=c(1,3))
  hist(a3[,i], xlab = uspecies[i], main= "Hunting")
  hist(a4[,i], xlab = uspecies[i], main="New log")
  hist(a5[,i], xlab = uspecies[i], main="Old log")
  browser()
  }
Q # quits the browser() function

# Baseline estimates of species-specific occupancy and detection in one of the habitat types
species.occ = fit$sims.list$u
species.det = fit$sims.list$v

# Are effects of new logging bigger than hunting in absolute value?
prob.newlog.hunt = matrix(0,nrow=length(uspecies),ncol=1)
rownames(prob.newlog.hunt)=uspecies
for (i in 1:length(uspecies)){prob.newlog.hunt[i]=(round((length(which((abs(a4[,i])-abs(a3[,i]))>= 0)))/nrow(psim),digits=4))*100}

# Are effects of old logging bigger than hunting in absolute value?
prob.oldlog.hunt = matrix(0,nrow=length(uspecies),ncol=1)
rownames(prob.oldlog.hunt)=uspecies
for (i in 1:length(uspecies)){prob.oldlog.hunt[i]=(round((length(which((abs(a5[,i])-abs(a3[,i]))>= 0)))/nrow(psim),digits=4))*100}

# Does hunting have a greater negative impact than new logging?
prob.neg.hunt1 = matrix(0,nrow=length(uspecies),ncol=1)
rownames(prob.neg.hunt1)=uspecies
for (i in 1:length(uspecies)){prob.neg.hunt1[i]= (round((length(which((a3[,i]<a4[,i]) & a3[,i]<0)))/nrow(psim),digits=4))*100}
prob.neg.hunt1

# Does hunting have a greater negative impact than old logging?
prob.neg.hunt2 = matrix(0,nrow=length(uspecies),ncol=1)
rownames(prob.neg.hunt2)=uspecies
for (i in 1:length(uspecies)){prob.neg.hunt2[i]= (round((length(which((a3[,i]<a5[,i]) & a3[,i]<0)))/nrow(psim),digits=4))*100}
prob.neg.hunt2

######################## GRAPHING RESULTS ##############################################
Nsite= fit$sims.list$Nsite
#hunt is nonscaled param
#hunt1 is scaled param
#Not super informative
#plot(hunt,apply(Nsite,2,mean))
u  = fit$sims.list$u
a1 = fit$sims.list$a1
a2 = fit$sims.list$a2
a3 = fit$sims.list$a3
a4 = fit$sims.list$a4
a5 = fit$sims.list$a5
#a6 = fit$sims.list$a6
#a7 = fit$sims.list$a7

#-------------------------- Line graph of hunting effects in each habitat
h.dummy=seq(0,0.25,0.005)
hhunt1=(h.dummy-mhunt)/sdhunt

#unlogged forest
psi=matrix(0,nrow=length(uspecies1),ncol=dim(u)[1])
richness = matrix(0, nrow=length(hhunt1), ncol=dim(u)[1])
for (j in 1:length(hhunt1)) {
    for(i in 1: length(uspecies1)) {
      for (k in 1:dim(u)[1]) {
      	psi[i,k]=plogis(u[k,i] +a1[k,i]*mean(elev1) +a2[k,i]*(mean(elev1)^2) +a3[k,i]*hhunt1[j] +a4[k,i]*0 
			+a5[k,i]*0 ) #+a6[k,i]*0 +a7[k,i]*0)  
      } }
      richness[j,] = apply(psi,2,sum)  
}   
richness1<-cbind(apply(richness,1,mean),apply(richness,1,quantile,0.975),apply(richness,1,quantile,0.025))

#newlogged forest
psi = matrix(0, nrow=length(uspecies1), ncol=dim(u)[1])
richness = matrix(0, nrow=length(hhunt1), ncol=dim(u)[1])
for (j in 1:length(hhunt1)) {
    for(i in 1: length(uspecies1)) {
      for (k in 1:dim(u)[1]) {
      	psi[i,k]=plogis(u[k,i] +a1[k,i]*mean(elev1) +a2[k,i]*(mean(elev1)^2) +a3[k,i]*hhunt1[j] +a4[k,i]*1 
			+a5[k,i]*0 ) #+a6[k,i]*hhunt1[j]*1 +a7[k,i]*0)  
      } }
      richness[j,] = apply(psi,2,sum)  
}
richness2<-cbind(apply(richness,1,mean),apply(richness,1,quantile,0.975),apply(richness,1,quantile,0.025))

#oldlogged forest
psi = matrix(0, nrow=length(uspecies1), ncol=dim(u)[1])
richness = matrix(0, nrow=length(hhunt1), ncol=dim(u)[1])
for (j in 1:length(hhunt1)) {
    for(i in 1: length(uspecies1)) {
      for (k in 1:dim(u)[1]) {
      	psi[i,k]=plogis(u[k,i] +a1[k,i]*mean(elev1) +a2[k,i]*(mean(elev1)^2) +a3[k,i]*hhunt1[j] +a4[k,i]*0 
			+a5[k,i]*1 ) #+a6[k,i]*0 +a7[k,i]*hhunt1[j]*1)  
      } }
      richness[j,] = apply(psi,2,sum)  
}
richness3<-cbind(apply(richness,1,mean),apply(richness,1,quantile,0.975),apply(richness,1,quantile,0.025))

#everything
richness<-cbind(h.dummy,richness1,richness2,richness3) # unlogged, newlog, oldlog
   

#-------------------------- Line graph of elevation effects in unlogged forest with no hunting
elev.dummy=seq(0,2000,,50)
eelev1=(elev.dummy-melev)/sdelev
psi=matrix(0,nrow=length(uspecies1),ncol=dim(u)[1])
richness = matrix(0, nrow=length(eelev1), ncol=dim(u)[1])
for (j in 1:length(eelev1)) {
    for(i in 1: length(uspecies1)) {
      for (k in 1:dim(u)[1]) {
      	psi[i,k]=plogis(u[k,i] +a1[k,i]*eelev1[j] +a2[k,i]*(eelev1[j]^2) +a3[k,i]*0 +a4[k,i]*0 +a5[k,i]*0)  
      } }
      richness[j,] = apply(psi,2,sum)  
}   
richness4<-cbind(apply(richness,1,mean),apply(richness,1,quantile,0.975),apply(richness,1,quantile,0.025))
erichness<-cbind(elev.dummy,richness4) 


#-------------------------- Bar graph of logging effects with no hunting
# unlogged
psi=matrix(0,nrow=length(uspecies1),ncol=dim(u)[1])
richness = matrix(0, nrow=1, ncol=dim(u)[1])
for (j in 1:1) {
    for(i in 1: length(uspecies1)) {
      for (k in 1:dim(u)[1]) {
      	psi[i,k]=plogis(u[k,i] +a1[k,i]*eelev1[j] +a2[k,i]*(eelev1[j]^2) +a3[k,i]*0 +a4[k,i]*0 +a5[k,i]*0)  
      } }
      richness[j,] = apply(psi,2,sum)  
}   
richness5<-cbind(apply(richness,1,mean),apply(richness,1,quantile,0.975),apply(richness,1,quantile,0.025))

# new logged
psi=matrix(0,nrow=length(uspecies1),ncol=dim(u)[1])
richness = matrix(0, nrow=1, ncol=dim(u)[1])
for (j in 1:1) {
    for(i in 1: length(uspecies1)) {
      for (k in 1:dim(u)[1]) {
      	psi[i,k]=plogis(u[k,i] +a1[k,i]*eelev1[j] +a2[k,i]*(eelev1[j]^2) +a3[k,i]*0 +a4[k,i]*1 +a5[k,i]*0)  
      } }
      richness[j,] = apply(psi,2,sum)  
}   
richness6<-cbind(apply(richness,1,mean),apply(richness,1,quantile,0.975),apply(richness,1,quantile,0.025))

# old logged
psi=matrix(0,nrow=length(uspecies1),ncol=dim(u)[1])
richness = matrix(0, nrow=1, ncol=dim(u)[1])
for (j in 1:1) {
    for(i in 1: length(uspecies1)) {
      for (k in 1:dim(u)[1]) {
      	psi[i,k]=plogis(u[k,i] +a1[k,i]*eelev1[j] +a2[k,i]*(eelev1[j]^2) +a3[k,i]*0 +a4[k,i]*0 +a5[k,i]*1)  
      } }
      richness[j,] = apply(psi,2,sum)  
}   
richness7<-cbind(apply(richness,1,mean),apply(richness,1,quantile,0.975),apply(richness,1,quantile,0.025))

# everything
loggingrichness<-cbind(richness5,richness6,richness7) # unlogged, newlog, oldlog


#-------------------------- Melt plot of posterior distributions of beta coefficients
library(ggplot2)
library(denstrip)
library(lattice)
library(reshape)

mu.a1=fit$sims.list$mu.a1
mu.a2=fit$sims.list$mu.a2
mu.a3=fit$sims.list$mu.a3
mu.a4=fit$sims.list$mu.a4
mu.a5=fit$sims.list$mu.a5
toExport<-cbind(mu.a1,mu.a2,mu.a3,mu.a4,mu.a5)
#write.csv(toExport,file="beta_posteriors.csv")

preds<-as.data.frame(cbind(mu.a1,mu.a2,mu.a3,mu.a4,mu.a5))
# sort effects by median size
idx<-sort(abs(apply(preds,2,median)),index.return=T,decreasing=F)$ix
idx2<-c(1,2,3,4,5)

# apply and sort labels
labs=c('mu.a1','mu.a2','mu.a3','mu.a4','mu.a5')[idx2]

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
























