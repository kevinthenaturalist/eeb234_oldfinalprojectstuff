#submitted EEB234 project idea: Using primarily/mostly R, I will be implementing species distribution modeling 
#on Spea hammondii (western spadefoot toad) localities retrieved from GBIF. Using Maxent via the "dismo" package 
#in R I will use 19 bioclimate variables in the present and in past climate scenarios to model past distributions 
#of S. hammondii, as well as its congener Spea multiplicata, to see if their ranges once overlapped in the past to 
#allow for ancient hybridization (currently they do not overlap). Final product will be modeled species 
#distributions of both species for the Last Glacial Maximum (LGM) and Last Interglacial (LIG). with these 
#projections as models for extreme climatic scenarios of the past million years since S. hammondii and 
#S. multiplicata have diverged.
#GBIF downloads come with dozens of columns of information; it is easy to select the desired 
#columns (latitude and longitude only) using R, but I will also present a way to remove unwanted 
#columns using regular expressions.

# SDM with dismo/maxent etc.
# but first manipulate GBIF data with python to get appropriate columns and stuff. Regex etc.
# as guides: http://thebiobucket.blogspot.com/2011/11/retrieve-gbif-species-distribution-data.html
# http://www.molecularecologist.com/2013/04/species-distribution-models-in-r/
# http://cran.r-project.org/web/packages/dismo/dismo.pdf

#download localities using rgbif, then use python+regex to edit it to species, longitude, and latitude
#then follow this SDM R procedure: http://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf

#do maxent projections with JUST the southern pop of S. hammondii, and then compare with S. multiplicata projection

setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files")
spea <- read.table("Spea_hammondii_herpnetlocalities.csv", header=TRUE, sep=",")
head(spea)
spealoc <- spea[,2:3] #puts only 2nd and 3rd column (long and lat in this case) into new dataframe
spealoc

library(dismo)
library(rgbif)
library(maptools)
require(maps)
require(mapdata)

dim(spea) #dimensions of spea dataframe
colnames(spea) #names of columns in spea dataframe

#next steps: get just lat-long of points, plot on map
spgeo <- subset(spea, !is.na(Longitude) & !is.na(Latitude)) #leaves only those records with long&lat data
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-120,-110), ylim=c(25,50), axes=TRUE, col="light yellow")
points(spgeo$Longitude, spgeo$Latitude, col='purple', pch=20, cex=0.75)
points(spgeo$Longitude, spgeo$Latitude, col='blue', cex=0.75) #manual suggests plotting points again to add a border for better visibility
box()

require(raster)
BClim = #download 19 bioclim variables
  
#generate random subset of points to use as training points
#create buffers around points to generate pseudoabsences
#run maxent via dismo package
#plot output
#repeat
