zz <- file("HLC_output.Rout", open="wb")
sink(zz)
sink(zz, type = "message")


# Loading libraries
library(sf)
library(spdep)
library(tmap)
library(units)
library(Matrix)

# Importing Devon HLC shapefile
zipfile <- "https://archaeologydataservice.ac.uk/catalogue/adsdata/arch-2090-1/dissemination/zip/rawhlc.zip"
## subject to https://archaeologydataservice.ac.uk/advice/termsOfUseAndAccess.xhtml
## It can be cited by https://doi.org/10.5284/1032952

td <- tempdir()
download.file(zipfile, destfile=file.path(td, "rawhlc.zip"))
fls <- unzip(file.path(td, "rawhlc.zip"), exdir=td, overwrite=TRUE)
devon_hlc <- st_read(file.path(td, "rawhlc.shp"), crs=27700)

## clip to Torridge District boundary

torridge_bys <- "https://raw.githubusercontent.com/digital-land/boundary-collection/master/collection/local-authority/E07000046/index.geojson"
## subject to https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/
bdy <- st_read(torridge_bys)


## initially created with sf linked to GEOS 3.9.0dev, with OverlayNG 
## without OverlayNG, intersection fails because of topology errors
## requiring st_make_valid to repair some polygons

hlc <- st_intersection(st_make_valid(devon_hlc), st_transform(st_geometry(bdy), 27700))

## Create ordered factor


hlc$class <- ifelse(hlc$PERIOD1 %in% c("Modern", "Post-medieval"), hlc$PERIOD1, "Medieval")
hlc$class <- ordered(hlc$class, levels=c("Modern", "Post-medieval", "Medieval"))

## Map classes

HLC_map <- tm_shape(hlc) + 
  tm_fill("class", palette="viridis", title="Character\nTypes") + 
  tm_compass(position=c("left", "bottom")) + 
  tm_scale_bar(breaks=c(0, 5, 10), position=c(0.375, 0.00001))

# tiff("Torridge_HLC.tiff", width=15, height=15, units="cm", res=300)
# HLC_map
# dev.off()
jpeg("Torridge_HLC.jpeg", width=15, height=15, units="cm", res=300)
HLC_map
dev.off()

## Create neighbours
nb1 <- poly2nb(hlc, snap=4, row.names=as.character(hlc$ID))
hlc.nb <- nblag(nb1, 3) ## higher orders
hlc.mat <- as(nb2listw(nblag_cumul(hlc.nb), style="B"), "CsparseMatrix")

## Join-Count Statistics

jc.hlc <- vector(mode="list", length=length(hlc.nb))
jc.hlc.p <- vector(mode="list", length=length(hlc.nb))

for (i in 1:length(hlc.nb)) {
  jc.hlc[[i]] <- joincount.multi(hlc$class, nb2listw(hlc.nb[[i]]))
  jc.hlc.p[[i]] <- pnorm(jc.hlc[[i]][,4], lower.tail=FALSE)
}

## Exporting output

jcs <- do.call("rbind", jc.hlc)[-c(7, 14, 21),]
jcps <- do.call("c", jc.hlc.p)[-c(7, 14, 21)]

(jc_out <- data.frame(order=rep(c("First", "Second", "Third"), each=6), JCS=rownames(jcs), as.data.frame(cbind(jcs, pvalue=jcps)), row.names=NULL))

write.csv(jc_out, "jc_out.csv", row.names=FALSE)

###########################################
## Boots' LICD (from Bivand et al. 2017) ##
###########################################

# For higher-order neighbour weight

#### STEP 1: local composition
(p <- as.matrix(summary(hlc$class))/nrow(hlc)) #probabilities of each "type"
areas <- aggregate(st_area(hlc), list(hlc$class), sum)
areas$x <- set_units(areas$x, "km2")
areas$props <- drop_units(areas$x/sum(areas$x))
areas

adata <- as.numeric(hlc$class) #factor no longer necessary, now numeric

source("local_JC0.R")

res <- local_JC0(obj=hlc, lagsmat=hlc.mat, varname="class", numvar=adata, p=p)
local_comp <- res[[1]]
JC.pvalue_seq <- res[[2]]


#### STEP 2: local configuration

local_config <- matrix(0,length(adata),1)
colnames(local_config) <- c("cluster-dispersion")

for (j in 1:length(adata)){#for cluster is 1, for dispersion -1, otherwise 0
  if (min(JC.pvalue_seq[j,])<1-(1-0.05)^(1/3)){ ###CHANGE
    ifelse(which(JC.pvalue_seq[j,]==min(JC.pvalue_seq[j,]), arr.ind = T)==1,local_config[j]<- 1, ifelse(which(JC.pvalue_seq[j,]==min(JC.pvalue_seq[j,]), arr.ind = T)==3,local_config[j]<- -1, local_config[j] <- 0))
  } # clump 1,dispersion -1, other 0
}

# Combination of local composition and local configuration
Type <- character(length=length(adata))
C <- cbind(local_comp, local_config)
for (i in 1:length(adata)){
  ifelse(C[i,1] == 1 && C[i,2] == 1, Type[i] <- "Cluster",
         ifelse(C[i,1] == 1 && C[i,2] == 0, Type[i] <- "Clump",
                ifelse(C[i,1] == -1 && C[i,2] == 0, Type[i] <- "Outlier",
                       ifelse(C[i,1] == 0 && C[i,2] == -1, Type[i] <- "Dispersed",
                              ifelse(C[i,1] == -1 && C[i,2] == -1, Type[i] <- "Outlier in dispersion area",
                                     Type[i] <- "No cluster")))))
}


# Plot LICD - TIFF + JPEG
Type1 <-  Type
hlc$Type <- Type
is.na(Type1) <- Type1 == "No cluster"
hlc$Type1 <- factor(Type1)
LICD_map <- tm_shape(hlc) + 
  tm_fill("Type1", palette="viridis", title="LICD", textNA="No cluster") + 
  tm_compass(position=c("left", "bottom")) + 
  tm_scale_bar(breaks=c(0, 5, 10), position=c(0.375, 0.00001))
# tiff("Torridge_LICD.tiff", width=15, height=15, units="cm", res=300)
# LICD_map
# dev.off()
jpeg("Torridge_LICD.jpeg", width=15, height=15, units="cm", res=300)
LICD_map
dev.off()


# HLC & LICD

both <- LICD_map + tm_facets("class", nrow=2)

jpeg("Torridge_HLC_LICD.jpeg",width=30,height=25,units="cm",res=300)
both
dev.off()

# mapview installed from "r-spatial/mapview" after #336 #327 #323
library(mapview)
packageVersion("mapview")
if (unname(sf_extSoftVersion()["GDAL"]) >= "3.1.0") mapviewOptions(fgb = FALSE)
file.remove("HLC_map.html")
file.remove("HLC_map.zip")
cl <- mapview(hlc, zcol="class")
ty <- mapview(hlc, zcol="Type")
mapshot(cl + ty, url = paste0(getwd(), "/HLC_map.html"))
zip("HLC_map.zip", "HLC_map.html")
file.remove("HLC_map.html")



sessionInfo()
sf_extSoftVersion()
sink(type = "message")
sink()

