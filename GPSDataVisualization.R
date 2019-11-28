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
GPSdata <- read.delim("GPSdata2.txt", header = FALSE) #Change File Name
GPSdata <- as.vector(unlist(GPSdata, use.names=FALSE))
lat <- vector()
long <- vector()
timeStamp <- vector()
for(i in 1:length(GPSdata)){
  lat[i] <- as.numeric(substr(GPSdata[i], 30,39))
  long[i] <- as.numeric(substr(GPSdata[i], 41, 51))
  timeStamp[i] <- substr(GPSdata[i], 0, 19)
}
#Add data to df
GPSdf <- data.frame(lat,long, timeStamp)
GPSdf$timeStamp <- as.POSIXlt(GPSdf$timeStamp, format = "%d-%m-%Y %H:%M:%S") #compare to ct
GPSdf <- na.omit(GPSdf)
GPSdf <- GPSdf[GPSdf$long <= 0, ]
#Convert to military time
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

timeCheckPlot <- plot(c(1:1531),GPSdf$timeStamp$hour, xlab="Data Point", ylab="Hour of Day (Military Time)", axes = FALSE)
axis(2, at = c(0:23))
axis(1, at = c(0:1541*100))
#Calculate Location Variance
locationVariance <- log(sd(na.omit(lat))^2 + sd(na.omit(long))^2)

#Define Optimal amount of Clusters using NbClust and visualize
GPSdfNoTime <- subset(GPSdf, select = c(lat,long))
nbc <- NbClust(GPSdfNoTime, method = "kmeans", index = "gap")
bestNbClust <- as.numeric(nbc$Best.nc[1])
silPlot <- fviz_nbclust(GPSdfNoTime, kmeans, method = "silhouette", k.max = 24) + theme_minimal() + ggtitle("The Silhouette Plot")

#Execute kmeans with optimal number of clusters
cluster <- kmeans(GPSdfNoTime, centers = bestNbClust)

#Calculate Entropy
clusterLabels <- cluster[1]
GPSdf <- cbind(clusterLabels, GPSdf) #add clusters to df.
GPSdf <- GPSdf[, !duplicated(colnames(GPSdf))] #remove duplicat column if applicable
clusterFreq <- table(clusterLabels)
entropy <- 0
for(i in length(clusterFreq)){
  p <- clusterFreq[i]/length(clusterLabels[[1]])
  entropy <- entropy + -(p) * log2(p)
}
entropy = entropy[[1]]

#Calculate Homestay
GPSHomestay <- GPSdf[GPSdf$timeStamp$hour < 6, ]
clusterFreqHS <- sort(table(GPSHomestay$cluster), decreasing = TRUE)
homestayCluster <- as.numeric(rownames(as.matrix(clusterFreqHS))[1])
homeStayPct <- nrow(GPSdf[GPSdf$cluster == homestayCluster, ]) / nrow(GPSdf)

#Calculate Circadian Movement
hrsSinceStart <- (as.vector(as.POSIXct(GPSdf$timeStamp)) - as.vector(as.POSIXct(GPSdf$timeStamp))[1])/3600
freqLow = 23.5
freqHigh = 24.5
latvTime <- plot(hrsSinceStart, GPSdf$lat) #Visualize sinosoidal pattern
longvTime <- plot(hrsSinceStart, GPSdf$long)
lspPlotLat <- lsp(GPSdf$lat, hrsSinceStart,freqLow, freqHigh)
lspPlotLat <- lsp(GPSdf$lat, hrsSinceStart,1, 26)
lspPlotLong <- lsp(GPSdf$long, hrsSinceStart,freqLow, freqHigh)
Elat <- 1/(freqHigh - freqLow) * sum(lspPlotLat$power)
Elong <- 1/(freqHigh - freqLow) * sum(lspPlotLong$power)
circadianMovement <- log10(Elat + Elong)

#Data.frame of All Calculated Metrics
metrics <- round(c(locationVariance, numberOfClusters, entropy, homeStayPct * 100, circadianMovement), 3)
units <- c("N/A", "clusters", "bits", "%", "N/A")
metrics <- data.frame(metrics, units)
row.names(metrics) <-c("Location Variance", "Number of Clusters", "Entropy", "Homestay Percentage", "Circadian Movement")
col.names(metrics) <-c("Metric", "Value", "Unit")
#Cluster Visualization - No Map
fviz_cluster(cluster, GPSdf[c("long","lat")], ellipse.type = "norm") +
  theme_minimal()
clusterCenters <- as.data.frame(cluster[2])

#Set API Key for ggmaps
ggmap::register_google(key = "api key")
#Get Map
ggmapObj <- get_googlemap(center = c(lon = mean(clusterCenters$centers.long) , lat = mean(clusterCenters$centers.lat)),
                          zoom = 13, scale = 2,
                          maptype ='terrain',
                          color = 'color')

#Cluster Visualization w/ Map 
ggmap(ggmapObj) + 
  geom_point(aes(x = long, y = lat), data = GPSdf, alpha = 0.1, size = 2.5, colour = 'red') +
  geom_point(aes(x = centers.long, y = centers.lat), data = clusterMeans , alpha = 0.1, size = 10, colour = 'blue')

#Calculate cluster size bounds to be 500m
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


