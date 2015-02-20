setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files")

library(dismo)
library(rgbif)
library(maptools)
require(maps)
require(mapdata)

speaGBIFraw <- read.table("Spea_hammondii_rawGBIFoutput.txt", header=TRUE, fill = TRUE, sep = "\t", quote="") #raw GBIF data requires these arguments for importing
#contains locality points for Spea hammondii, the western spadefoot toad

#This function, gbifclean(), takes raw output from GBIF or Vertnet and converts
#it into the appropriate format for running in Maxent
#this function should only be used for a single species at a time
gbifclean <- function(gbifraw) {
  latcol <- grep("decimal[Ll]atitude", colnames(gbifraw)) #returns column number for Latitude values
  longcol <- grep("decimal[Ll]ongitude", colnames(gbifraw)) #returns column number for Longitude values
  namecol <- grep("^scientificName$", colnames(gbifraw)) #returns column number for species name(s)
  latlong <- gbifraw[,c(namecol, longcol, latcol)] #subsets the input data to just these 3 columns
  colnames(latlong)[2] <- "Longitude" #changes column name to "Longitude"
  colnames(latlong)[3] <- "Latitude" #changes column name to "Latitude"
  latlong[,1] <- "Spea hammondii" #CHANGE BASED ON INPUT SPECIES! Converts names in column 1 to proper name
  latlong <- latlong[complete.cases(latlong), ] #removes rows with any missing data
  latlong[,2] <- abs(latlong[,2]) * -1 #ensures that longitude values are negative; use as appropriate for region
  latlong[,3] <- abs(latlong[,3]) #ensures that latitude values are positive; use as appropriate for region
  latlong
}
speaGBIFclean <- gbifclean(speaGBIFraw)
head(speaGBIFclean)

#plotgbifclean() then takes the output from gbifclean() (or any dataframe 
#with that speciesname, longitude, latitude format) and plots it on a map
#may have to fix the projection of the "wrld_simpl" map...
plotgbifclean <- function(latlong) {
  data(wrld_simpl) #loads world map
  xlowlim <- min(latlong$Longitude) - 0 #uses values in data to set appropriate limits on plot axes
  xuplim <- max(latlong$Longitude) + 0
  ylowlim <- min(latolng$Latitude) - 0
  yuplim <- max(latlong$Latitude) + 0
  plot(wrld_simpl, xlim = c(xlowlim, xuplim), ylim = c(ylowlim, yuplim), axes=TRUE, ylab="Latitude", xlab="Longitude", cex.axis=0.7)
  points(latlong$Longitude, latlong$Latitude, col='purple', pch=20, cex=0.75) #adds points to map
  box() #puts box around plot
}
plotgbifclean(speaGBIFclean) #plotting the data may reveal aberrant points; these can be fixed later
#points look a bit off; map projection probably different from regular latlong.
#I think WorldClim data isn't projected--just latlong--so it shouldn't affect the Maxent analyses

#library(ggplot2) #I need to play around with this more. Maybe use ggmaps? Will be good to do this to 
#get used to using ggplot and to make more precise and more attractive figures/maps
#mp <- NULL #creates empty thingy
#mapWorld <- borders("world", color="gray50", fill="gray50", xlim=range(spgeo$Longitude), ylim=range(spgeo$Latitude))
#data(spgeo)
#map("usa", xlim=range(spgeo$Longitude), ylim=range(spgeo$Latitude)) + mp
#mp <- ggplot(spgeo, aes(x=Longitude, y=Latitude)) + mapWorld +geom_point(data=spgeo)
#mp
#ggplot(spealoc, aes(x=Longitude, y=Latitude)) + xlim(-125, -115) + ylim(27,39) + geom_point()
#?borders

require(raster)
BClim = #download 19 bioclim variables
  
#generate random subset of points to use as training points
#create buffers around points to generate pseudoabsences
#run maxent via dismo package
#plot output
#repeat
