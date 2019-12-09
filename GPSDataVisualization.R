#This file takes in raw GPS Data, and visualizes the data. Requires Google Maps Static API key to visualize the data
#Install the relevant libraries - do this one time
#install.packages("lubridate")
#install.packages("ggplot2")
#install.packages("data.table")
#install.packages("ggrepel")
#install.packages("dplyr")
#install.packages("data.table")
#install.packages("tidyverse")
#install.packages("ggmap")
#install.packages("raster")
#install.packages("NbClust")
#install.packages("factoextra")
#install.packages("lomb")

# Load the relevant libraries - do this every time
library(lomb)
#library(lubridate)
library(ggplot2)
library(dplyr)
#library(data.table)
library(ggrepel)
library(tidyverse)
library("ggmap")
library(raster)
library(NbClust)
library(factoextra)

#Load Latitude and Logitude data
setwd("C:/Users/VP1050/OneDrive/Documents/1Fall 2019/IE 431 - Senior Design/MyData - App") #Change Diretory
GPSdata <- read.delim("GPSdata2.txt", header = FALSE) #Data input, Change File Name
GPSdata <- as.vector(unlist(GPSdata, use.names=FALSE)) #Each element of vector becomes 1 line of text.
lat <- vector() #Create empty vector for latitude data 
long <- vector() #Create empty vector for longitude data
timeStamp <- vector() #Create empty vector for timestamp data
for(i in 1:length(GPSdata)){ #iterate through each line
  lat[i] <- as.numeric(substr(GPSdata[i], 30,39)) #Make longitude data numeric
  long[i] <- as.numeric(substr(GPSdata[i], 41, 51)) #Make latitdue data numeric
  timeStamp[i] <- substr(GPSdata[i], 0, 19) #Substring cretes timestamp
}
#Add data to df
GPSdf <- data.frame(lat,long, timeStamp) #Dataframe of previously created vectors
GPSdf$timeStamp <- as.POSIXlt(GPSdf$timeStamp, format = "%d-%m-%Y %H:%M:%S") #"POSIXlt" format is a 
  #named list of vectors representing each time unit. 
  #See https://stat.ethz.ch/R-manual/R-devel/library/base/html/DateTimeClasses.html
GPSdf <- na.omit(GPSdf) #Omit NAs
GPSdf <- GPSdf[GPSdf$long <= 0, ] #Ensure that longitude col is negative in case of any parsing issues.
#Convert to military time - Needed if data does not have AM / PM and is not military time. 
#Running this loop assumes data starts at PM
PM <- TRUE
MN <- FALSE
for(i in 1:length(GPSdf$timeStamp)){
  if (PM == TRUE) {
    GPSdf$timeStamp$hour[i] <- GPSdf$timeStamp$hour[i] + 12
    if (GPSdf$timeStamp$hour[i + 1] == 12){
      PM <- FALSE
      MN <- TRUE
    }
  } else if (MN == TRUE) {
    GPSdf$timeStamp$hour[i] <- GPSdf$timeStamp$hour[i] - 12
    if(GPSdf$timeStamp$hour[i + 1] == 1){
      MN <- FALSE
    }
  } else {
    if(GPSdf$timeStamp$hour[i] == 12 & GPSdf$timeStamp$hour[i + 1] == 1){
      PM <- TRUE
    }
  }
}

timeCheckPlot <- plot(c(1:1531),GPSdf$timeStamp$hour, xlab="Data Point", ylab="Hour of Day (Military Time)", axes = FALSE) #Checks that data points are taken at least every hour
axis(2, at = c(0:23)) #Change plot axis
axis(1, at = c(0:1541*100)) #Change plot axis

#Calculate Location Variance
locationVariance <- log(sd(na.omit(lat))^2 + sd(na.omit(long))^2)

#Define Optimal amount of Clusters using NbClust and visualize
GPSdfNoTime <- subset(GPSdf, select = c(lat,long)) #subset only lat and long data
nbc <- NbClust(GPSdfNoTime, method = "kmeans", index = "gap") #Calculates optimal number of clusters
bestNbClust <- as.numeric(nbc$Best.nc[1]) #Get numeric value of optimal number of clusters
silPlot <- fviz_nbclust(GPSdfNoTime, kmeans, method = "silhouette", k.max = 24) + theme_minimal() + ggtitle("The Silhouette Plot") 
  #Visualizes optimal number of clusters

#Execute kmeans with optimal number of clusters
cluster <- kmeans(GPSdfNoTime, centers = bestNbClust)

#Calculate Entropy
clusterLabels <- cluster[1] #Cluster labels of every data point
GPSdf <- cbind(clusterLabels, GPSdf) #add cluster labels to df.
GPSdf <- GPSdf[, !duplicated(colnames(GPSdf))] #remove duplicate column if applicable due to cbind
clusterFreq <- table(clusterLabels) #tally number of data points in each cluster
entropy <- 0
for(i in length(clusterFreq)){ #Saeb et al summation for entropy formula
  p <- clusterFreq[i]/length(clusterLabels[[1]])
  entropy <- entropy + -(p) * log2(p)
}
entropy = entropy[[1]] #get numeric value only 

#Calculate Homestay
GPSHomestay <- GPSdf[GPSdf$timeStamp$hour < 6, ] #Get times between 12AM and 6AM 
clusterFreqHS <- sort(table(GPSHomestay$cluster), decreasing = TRUE) 
  #Tally subset of points in each cluster and sort with highest value first.
homestayCluster <- as.numeric(rownames(as.matrix(clusterFreqHS))[1]) #get home cluster
homeStayPct <- nrow(GPSdf[GPSdf$cluster == homestayCluster, ]) / nrow(GPSdf) #get homestay percentage

#Calculate Circadian Movement
hrsSinceStart <- (as.vector(as.POSIXct(GPSdf$timeStamp)) - as.vector(as.POSIXct(GPSdf$timeStamp))[1])/3600
  #Converts all timesteamps into hourse since start of data collection.
freqLow = 23.5
freqHigh = 24.5
latvTime <- plot(hrsSinceStart, GPSdf$lat) #Visualize sinosoidal pattern across lattitude
longvTime <- plot(hrsSinceStart, GPSdf$long) 
lspPlotLat <- lsp(GPSdf$lat, hrsSinceStart,freqLow, freqHigh) #Lomb Scargle Periodogram lat
lspPlotLong <- lsp(GPSdf$long, hrsSinceStart,freqLow, freqHigh) #Lomb Scargle Periodogram long
Elat <- 1/(freqHigh - freqLow) * sum(lspPlotLat$power) #Saeb et al 
Elong <- 1/(freqHigh - freqLow) * sum(lspPlotLong$power)
circadianMovement <- log10(Elat + Elong) #Saeb et al formula 

#Data.frame of All Calculated Metrics
metrics <- round(c(locationVariance, bestNbClust, entropy, homeStayPct * 100, circadianMovement), 3)
  #Round the vector of all calculated metrics to 3 decimal places
units <- c("N/A", "clusters", "bits", "%", "N/A") #units of each metric 
metrics <- data.frame(metrics, units) #Dataframe of values
row.names(metrics) <-c("Location Variance", "Number of Clusters", "Entropy", "Homestay Percentage", "Circadian Movement")

#Cluster Visualization - No Map
fviz_cluster(cluster, GPSdf[c("long","lat")], ellipse.type = "norm") +
  theme_minimal()
clusterCenters <- as.data.frame(cluster[2])

#Set API Key for ggmaps - will need an API key
ggmap::register_google(key = "AIzaSyBHLEEw0R4syKD6FiygfegB3iPjkRW2IHw")
#Get Map with appropirate parameters 
ggmapObj <- get_googlemap(center = c(lon = mean(clusterCenters$centers.long) , lat = mean(clusterCenters$centers.lat)),
                          zoom = 13, scale = 2,
                          maptype ='terrain',
                          color = 'color')

#Cluster Visualization w/ Map 
ggmap(ggmapObj) + 
  geom_point(aes(x = long, y = lat), data = GPSdf, alpha = 0.1, size = 2.5, colour = 'red') +
  geom_point(aes(x = centers.long, y = centers.lat), data = clusterCenters , alpha = 0.1, size = 20, colour = 'blue')


#Calculate cluster size bounds to be 500m, This method was used to calculate the best number of clusters
  #prior to using sihlouette width. 

#i <- TRUE #True if cluster distance greater than 500m 
# numClusters <- 1
# while (i == TRUE){
#   cluster <- kmeans(GPSdf, centers = numClusters)
#   label <- cluster[1]
#   centers <- cluster[2]
#   numClusters <- numClusters + 1
#   n <- 1
#   k <- TRUE
#   while (k == TRUE && n <= length(label)){
#     v <- label[[1]][n]
#     c <- centers[[1]][v, ]
#     GPSrow <- as.numeric(GPSdf[n, ])
#     if (pointDistance(GPSrow,c, lonlat = TRUE) > 500){
#       k <- FALSE
#     }
#     n <- n + 1 
#     }
#     if (k && n >= length(label)){
#       i <- FALSE
#     }
# }


