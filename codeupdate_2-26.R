#Using R to run a Species Distribution Model with climate projections for the western spadefoot toad, _Spea hammondii_
setwd("C:/Users/Kevin/Google Drive/UCLA Courses or Lab meetings etc/EEB 234/Final Project files")


library(ENMeval)
library(dismo)
library(maptools)
library(maps)
library(mapdata)
library(ggplot2)
library(rJava)
library(rgdal)


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
write.table(speaGBIFclean, file = "Speahammondii_GBIF_longlat.txt", row.names=FALSE)

head(speaGBIFclean)

#plotgbifclean() then takes the output from gbifclean() (or any dataframe 
#with that speciesname, longitude, latitude format) and plots it on a map
#may have to fix the projection of the "wrld_simpl" map...
plotgbifclean <- function(latlong) {
  data(wrld_simpl) #loads world map
  xlowlim <- min(latlong$Longitude) - 0 #uses values in data to set appropriate limits on plot axes
  xuplim <- max(latlong$Longitude) + 0
  ylowlim <- min(latlong$Latitude) - 0
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

#Using dismo! Referring here to the dismo vignette
spealoc <- speaGBIFclean[,2:3] #only need lat and long columns for maxent
head(spealoc)
plot(spealoc)
plot(wrld_simpl, add=T, border='blue', lwd=2)

#delete aberrant points
speaLocfix <- subset(spealoc, subset=((Latitude > 28.5) & (Longitude < -120 | Latitude < 39.3) & (Longitude < -118.5 | Latitude < 36.7)))
#removes points outside Spea hammondii's known range; likely other Spea species that look very similar
dim(speaLocfix)
plot(speaLocfix)
plot(wrld_simpl, add=T, border='blue', lwd=2)

#subsampling data - not required; other methods of dealing with bias
#could make this a function...
speaspdf <- SpatialPointsDataFrame( speaLocfix[ c("Longitude", "Latitude") ], data = data.frame(speaLocfix), proj4string = CRS("+proj=longlat +datum=WGS84")) #makes a "spatial points dataframe"...
r <- raster(speaspdf) #converting to raster only works if I project the localities first using line above???
r
?SpatialPointsDataFrame

res(r) <- 0.25
?res
r <- extend(r, extent(r)+1)
speasel <- gridSample(speaspdf, r, n=1) #produces subsample of n=1 point from each raster grid; this will be occurrence data input
p <- rasterToPolygons(r)
plot(p, border='gray')
points(speaspdf)
points(speasel, cex=1, col='red', pch='x')
#write.table(speasel, file="") #creates a file of this subsampled data
?gridSample
?brick
#speasel, the subsampled set, can be used as the ENMeval input...

###generating pseudoabsence data
circ <- circles(speaspdf, d=60000, lonlat=TRUE) #creates circles of 50km radius around presence points; will draw pseudoabsences from these buffers
discirc <- gUnaryUnion(circ@polygons) #dissolve circle features into one feature
discirc2 <- spTransform(discirc, CRS("+proj=longlat +datum=WGS84")) #changes CRS to match that of speaspdf
#create smaller circles to erase from bigger circles to limit sampling near presence points
smallcirc <- circles(speaspdf, d=10000, lonlat=TRUE)
dissmallcirc <- gUnaryUnion(smallcirc@polygons)
dissmallcirc2 <- spTransform(dissmallcirc, CRS("+proj=longlat +datum=WGS84"))
discirc3 <- erase(discirc2, dissmallcirc2) #removes circles of radius disssmallcirc from discirc2 polygons
#sample pseudoabsence points from discirc3
discirc3rast <- rasterize(discirc3, bio1Sham) #must convert to raster before you can use another raster as a mask
#the point of this is to trim the layer at the coast so that there are no pseudoabsences in the ocean
maskedcirc <- mask(discirc3rast, bio1Sham)
plot(maskedcirc)

plot(bio1Sham)
plot(discirc3rast)
discirc3cut <- rasterToPolygons(maskedcirc)
plot(discirc3cut)
?rasterToPolygons
pseudoabs <- spsample(discirc3cut, 500, type='random', iter=25) #random sample of points from within the circles; used as background points for maxent
#tooclose <- zerodist2(circsamp, speaspdf, zero=5.0) #can use this or can use the erase() method below
#pseudoabs <- circsamp[-c(tooclose),]
points(pseudoabs)
#convert pseudoabs to a matrix
pseudoabspts <- as.matrix(pseudoabs)


bio1 <- raster('bioclim2.5/wc2-5/bio1.bil') #will use this as the mask to get rid of points in the ocean
bio1Sham <- crop(bio1, ShamRangenarrow)

#then clip/mask to remove points in ocean...
par(mfrow=c(1,1)); plot(discirc3); points(pseudoabs, pch=19, cex=0.005); points(speaspdf, col='red', pch=19, cex=0.3)

require(raster)
BClim2_5 = getData("worldclim", var="bio", res=2.5, path="bioclim2.5/") #download 19 bioclim variables, 2.5arcmin resolution
#may be easier/faster to use my pre-cropped files on other computer
ShamRange = extent(-125.5, -90.5, 15.5, 48.5) #crop Bioclim layers to this extent
#do I crop to an extent to compare with other Spea species?
ShamRangenarrow = extent(-125.5, -110, 25, 45)
bclim2.5Shamnarrow = crop(BClim2_5, ShamRangenarrow)
writeRaster(bclim2.5Shamnarrow, filename="bioclim2.5/ShamnarrowBC_2.5.grd", overwrite=T)
bclim2.5Shamnarrow = brick("bioclim2.5/ShamnarrowBC_2.5.grd")
library(ENMeval)

pseudoabscoords <- coordinates(pseudoabs) #pulls coordinates as matrix from a SpatialPoints object

#shamoccENM <- cbind(speaLocfix[,1], speaLocfix[,2]) #have to convert long/lat to matrix in this way before running ENMevaluate
shamnarrowbcENMeval <- ENMevaluate(speasel, bclim2.5Shamnarrow, bg.coords=pseudoabscoords, method="randomkfold", kfolds = 2)
#can use n.bg to set random background points; may be worth pursuing method in molecularecologist.com of 
#getting background points only from areas withon xx km of a presence point
#bg.coords is user-inputted pseudoabsences, which I painstakingly generated...

data(shamnarrowbcENMeval)
shamnarrowbcENMeval@results
plot(shamnarrowbcENMeval@predictions[[which (shamnarrowbcENMeval@results$delta.AICc == 0) ]])
plot(shambcENMeval@predictions[[which (shambcENMeval@results$delta.AICc == 0) ]], xlim=c(-125.5, -110), ylim=c(25, 40))
points(shambcENMeval@occ.pts)
shamnarrowbcENMeval@overlap #? see manual about this
#***I should definitely subset the data before running, such that only 1 or 2 points are in any given cell...
sham_p = kfold(speasel, 5) #vector of group assignments, splitting speasel into 5 eval groups
sham_a = kfold(pseudoabscoords, 5) #same for the background points
test = 2
colnames(pseudoabscoords)[1] <- "Longitude"
colnames(pseudoabscoords)[2] <- "Latitude"
?stack

#making a stack of Future layers for climate projection
#CMIP2070CN30s <- getData('CMIP5', var='bio', res=2.5, rcp=85, model='CN', year=70, path='bioclim2070_CN_rcp85') #future; function goes to wrong URL...
#CMIP2070 <- stack(bc2070files) #to easily do stack this way, have to move the files to the working directory...
#cmip2070sham <- crop(CMIP2070, ShamRangenarrow)

#cmipbio1 <- crop(raster("cn85bi701.tif"), ShamRangenarrow, filename="bio1.tif")

#function to crop and rename 19 tif layers for use as projectionlayers
bc2070files <- list.files(path='bioclim2070_CN_rcp85/cmip5/2-5m/cn85bi70') #download and bring files into R manually
#have to then move these files into the working directory for the following loop to work
for (i in 1:19) {
  crop(raster(bc2070files[i]), ShamRangenarrow, filename=paste("bio",i,".tif", sep=""), overwrite=TRUE)
}

#change raster format from TIFF to .bil (ugh this is annoying)
#again have to move files into working directory for loop to work. move to own folder after
bc2070tifs <- list.files(path='bioclim2070_CN_rcp85')
for (i in 1:19) {
  writeRaster(raster(bc2070tifs[i]), filename=paste("bio",i,".bil",sep=""), format="EHdr")
}

#turn these into a stack...
bc2070bil <- list.files(path="bioclim2070_CN_rcp85", ".bil$") #makes list of only filenames that end in .bil
bc2070bilstack <- stack(bc2070bil) #makes raster stack
writeRaster(bc2070bilstack, filename="bioclim2070_CN_rcp85.grd", overwrite=T) #saves stack to a file

#names(cmipbio1)
#writeRaster(cmip2070sham, filename='CMIP2070_CN2-5Sham.grd', overwrite=T) #re-export the raster stack as a .grd. This takes a long time; crop everything first
#library(rgdal)
#cmip2070shampr <- projectRaster(cmip2070sham, crs=projection(bclim2.5Shamnarrow)) #change projection for consistency
#cmip2070 <- brick("CMIP2070_CN2-5.grd") #load if you need to recall this
#ShamCmip2070 <- cmip2070shampr #making copy for testing...
#names(ShamCmip2070) <- names(bclim2.5Shamnarrow) #maxent() requires projectionlayers to be titled "bio1, bio2,..."
#names(ShamCmip2070)

train_p = speasel[sham_p!=test, c("Longitude", "Latitude")]
train_a = pseudoabscoords[sham_a!=test, c("Longitude", "Latitude")]
test_p = speasel[sham_p==test, c("Longitude", "Latitude")]
test_a = pseudoabscoords[sham_a!=test, c("Longitude", "Latitude")]
maxentrun <- maxent(bclim2.5Shamnarrow, speasel, a=pseudoabscoords, path="maxentrun", args=c("-J", "-P"))
#for projecting onto other layers/scenarios: add arg "projectionlayers=layersdirectory"
maxentfuture <- maxent(bclim2.5Shamnarrow, speasel, a=pseudoabscoords, path="maxentrunfuture", args=c("-J", "-P", "projectionlayers=bioclim2070_CN_rcp85"))
maxentfuturenopseudo <- maxent(bclim2.5Shamnarrow, speasel, path="maxentrunfuturenopseudo", args=c("-J", "-P", "projectionlayers=bioclim2070_CN_rcp85"))
#try a different future model... this one seems unusual...


e = evaluate(test_p, test_a, maxentrun, bclim2.5Shamnarrow)
e
par(mfrow=c(1,2))
pred_me = predict(maxentrun, bclim2.5Shamnarrow) #could I use this to do a climate projection?
pred_metest = predict(maxentrun, bclim2.5Shamnarrow)
plot(pred_me, 1, cex=0.5, legend=T, mar=par("mar"), main="Predicted presence of western spadefoots")
plot(pred_metest, 1, cex=0.5, legend=T, mar=par("mar"), main="uncertain future prediction thing, don't trust")
pred_fut2 <- predict(maxentfuture, bclim2.5Shamnarrow)
plot(pred_fut2, 1, cex=0.5, legend=T, mar=par("mar"), main="futuuuure")
#the run variable matters; the raster stack in there doesn't... I don't think...

?predict
plot(maxentrun)
response(maxentrun)
