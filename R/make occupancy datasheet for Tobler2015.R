rm(list=ls(all=TRUE))
data.photo <- read.csv("data/All_mammal_photo.csv", header=TRUE)
sitecovs <- read.csv("data/sitecov_temp_Cheng_reorded.csv")
list <- data.frame(aggregate(data.photo$n ~ data.photo$species + data.photo$camera, FUN =sum)) 
names(list) <- c("spp", "camera", "n")
spp.li <- as.character(unique(list$spp))

occ_data <- matrix(999, nrow = 22*115, ncol = 3 )
occ_data[,1] <-  tolower(sitecovs$Station)
j = 1
k = 115
for (i in 1:22) {      # make spp column, each spp repeat for 115 times
  occ_data[j:k,2] = spp.li[i]  
  j = j + 115
  k = k + 115
} 

for (i in 1:nrow(list)) {  # make occ column, if station & spp matches 
  n <- which(occ_data[,2] == list[i,1] & occ_data[,1] == list[i,2])
  occ_data[n,3] <- list[i,3]
}

occ_data[which(occ_data[,3] == 999),3] <- 0
occ_data <- data.frame(occ_data, stringsAsFactors = FALSE)
names(occ_data) <- c("Station","Species","Count")
occ_data$Count <- as.numeric(occ_data$Count)
sum(list[,3]) # check 
sum(occ_data[,3]) #check


write.csv(occ_data,file="/Users/chencheng/Desktop/data/occupancy/data/occ_data_by_station.csv", row.names = FALSE)
