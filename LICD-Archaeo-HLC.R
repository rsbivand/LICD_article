# Set workspace
library(rstudioapi)
cat("\f")
rm(list=ls())
current_path<-getActiveDocumentContext()$path
setwd(dirname(current_path));print(getwd())

# Loading libraries
library(maptools)
library(rgdal)
library(spdep)
library(ggplot2)
library(ggspatial)
library(gridExtra)
library(foreach) #necessary for parallel computing - many spatial units
library(doParallel)

# Importing shapefile
## Character types aggregated into 3 cats: 
## medieval, post-medieval, modern
hlc<-readOGR("rawhlc/Torridge_HLC_3b.shp")
hlc$class<-ifelse(hlc$PERIOD1=="Modern",1,ifelse(hlc$PERIOD1=="Post-medieval",2,3))
hlc$class[is.na(hlc$class)]<-3

# Create neighbours
hlc.nb <- nblag(poly2nb(hlc),3) ##1 order
hlc.mat <- as(listw2mat(nb2listw(hlc.nb[[1]],style="B")),"sparseMatrix") + 
  as(listw2mat(nb2listw(hlc.nb[[2]],style="B")),"sparseMatrix")+ 
  as(listw2mat(nb2listw(hlc.nb[[3]],style="B")),"sparseMatrix")

# Join-Count Statistics

## Set column with analysed data in datafile
clm<-length(names(hlc))
adata<-factor(as.data.frame(hlc[,clm])[,1]) ##object with "levels" (factor), first column is not $cat

jc.hlc<-list()
jc.hlc.p<-list()

for(i in 1:3){
  jc.hlc[[i]]<-joincount.multi(adata,nb2listw(hlc.nb[[i]]))
}

for(i in 1:3){
  jc.hlc.p[[i]]<-pnorm(jc.hlc[[i]][,4],lower.tail=F)
}

## Exporting output

capture.output(cbind(jc.hlc[[1]],pvalue=as.vector(jc.hlc.p[[1]])),file="Output/JCS_HLC_lag1.txt")
capture.output(cbind(jc.hlc[[2]],pvalue=as.vector(jc.hlc.p[[2]])),file="Output/JCS_HLC_lag2.txt")
capture.output(cbind(jc.hlc[[3]],pvalue=as.vector(jc.hlc.p[[3]])),file="Output/JCS_HLC_lag3.txt")


###########################################
## Boots' LICD (from Bivand et al. 2017) ##
###########################################

# For higher-order neighbour weight

#### STEP 1: local composition
p <- (as.matrix(summary(adata)))/length(adata) #probabilities of each "type"
adata <- as.numeric(adata) #factor no longer necessary, now numeric

### Parallel computing block is starting ###

# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster of parallel computing
cl <- makeCluster(no_cores)
registerDoParallel(cl)

system.time({
  ## Routine 1 in parallel computing ##
  # cluster has 5 columns: 1 - number of units of the "type" as the unit j in the "window
  # 2 - probability of the "type", 3 - "window" size, 4 - probability of "absence" of "type"
  
  cluster <- foreach (i=1:length(adata), .combine=rbind) %dopar% {
    require(Matrix)
    c1 <- (hlc.mat[i,] %*% ifelse(adata==adata[i],1,0))+1
    c2 <- p[adata[i]]
    c3 <- sum(hlc.mat[i,])+1
    c4 <- sum(dbinom(0:c1, c3, c2)) ###CHANGE!
    c5 <- sum(dbinom(c1:c3, c3, c2)) ###CHANGE!
    cbind(c1, c2, c3, c4, c5)
  }
  
  cluster[is.nan(cluster)]<- 1
  ## End of routine 1 ##
  
  ### Cluster-outlier analysis -> result of local composition ###
  local_comp <- ifelse(cluster[,5]<1-(1-0.05)^(1/2), 1, ifelse(cluster[,4]<1-(1-0.05)^(1/2), -1, 0))
  # 1 for cluster, -1 for outlier, 0 - neither cluster, nor outlier
  
  ## Routine 2 in parallel computing ##
  # We calculate of bivariate local JC's: 1 for units of the same "type" as j, 0 otherwise
  JC.pvalue <- foreach (j=1:length(adata), .combine=rbind) %dopar% {
    require(Matrix)
    require(spdep)
    hlc.mat.1 <- hlc.mat[j,] #extracting a row from weights matrix
    ktore <- which(hlc.mat.1!=0, arr.ind = T) #looking for neighbours of j
    hlc.1 <- hlc[c(j,ktore),] #adding j to the list of its neighbours
    #datafile.1.B <- nb2listw(poly2nb(datafile.1), style="B") #weights matrix for j and its neighbours
    hlc.1.B <- mat2listw(hlc.mat[c(j,ktore), c(j,ktore)], row.names = NULL, style="B")
    datawork <- as.data.frame(hlc.1) #working data
    adata01 <- matrix(0,length(ktore)+1,1) #"types" of units and j "type"
    A <- matrix()
    
    for (i in 1:length(ktore)+1){
      adata01[i,1] <- ifelse(datawork[i,clm]==datawork[1,clm],1,0)
      #"Type" of unit is 1 if it is as the "type" of j, 0 otherwise
    } 
    adata01[1,1]<-1 
    
    if (any(adata01 == 0)){
      #if any units is different "type" from j, then I calculate JC
      A <- joincount.multi(as.factor(adata01), hlc.1.B)
      JC.pvalue <- cbind(1-pnorm(A[2,4]), 1-pnorm(A[1,4]), pnorm(A[3,4])) #p-value JC
    } else {
      #if all are "type" as j, then JC is out of sense
      JC.pvalue <- cbind(0, 1, 1)
    }
    
  }
  
  JC.pvalue[is.nan(JC.pvalue)]<- 1
  colnames(JC.pvalue) <- c("1:1", "0:0", "1:0")
  ## End of routine 2 ##
  
  stopCluster(cl)
  ### The end of parallel computing block ###
})

#### STEP 2: local configuration

local_config <- matrix(0,length(adata),1)
colnames(local_config) <- c("cluster-dispersion")

for (j in 1:length(adata)){#for cluster is 1, for dispersion -1, otherwise 0
  if (min(JC.pvalue[j,])<1-(1-0.05)^(1/3)){ ###CHANGE
    ifelse(which(JC.pvalue[j,]==min(JC.pvalue[j,]), arr.ind = T)==1,local_config[j]<- 1, ifelse(which(JC.pvalue[j,]==min(JC.pvalue[j,]), arr.ind = T)==3,local_config[j]<- -1, local_config[j] <- 0))
  } # clump 1,dispersion -1, other 0
}

# Combination of local composition and local configuration
Type <- character(length=0)
C <- cbind(local_comp, local_config)
for (i in 1:length(adata)){
  ifelse(C[i,1] == 1 && C[i,2] == 1, Type[i] <- "Cluster",
         ifelse(C[i,1] == 1 && C[i,2] == 0, Type[i] <- "Clump",
                ifelse(C[i,1] == -1 && C[i,2] == 0, Type[i] <- "Outlier",
                       ifelse(C[i,1] == 0 && C[i,2] == -1, Type[i] <- "Dispersed",
                              ifelse(C[i,1] == -1 && C[i,2] == -1, Type[i] <- "Outlier in dispersion area",
                                     Type[i] <- "No cluster")))))
}

## Preparing for shapefile
ID.data <- sapply(slot(hlc, "polygons"), function(x) slot(x, "ID"))
hlc.id<-cbind(hlc,ID.data)
names(hlc.id)[length(names(hlc.id))] <- c("ID_object")
X<-cbind(as.data.frame(hlc.id@data["ID_object"]),as.data.frame(cbind(C,Type)))
hlc.id@data = data.frame(hlc.id@data, X[match(hlc.id@data[,length(names(hlc.id))], X[,1]),])

#hlc.id <- cbind(hlc,ID.data)
#names(hlc.id)[length(names(hlc.id))] <- c("ID_object")

## Merging data and geometry
#X <- cbind(as.data.frame(hlc.id@data[,7]), as.data.frame(cbind(C, Type)))
#hlc.id@data = data.frame(hlc.id@data, X[match(hlc.id@data[,7], X[,1]),])
#names(hlc.id)[length(names(hlc.id))-1] <- c("local_config")

## Writing data
#writePolyShape(hlc.id,"Torridge_LICD")

# Plot HLC - TIFF + JPEG
tiff("Output/Torridge_HLC.tiff",width=12,height=10,units="cm",res=300)
ggplot()+
  layer_spatial(hlc,aes(fill=as.factor(class)),colour=NA)+
  scale_fill_viridis_d(name="Character\nTypes",label=c("Modern","Post-\nMedieval","Medieval"))+
  theme_light()+
  annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                         pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                         style=north_arrow_fancy_orienteering())
dev.off()

jpeg("Output/Torridge_HLC.jpeg",width=12,height=10,units="cm",res=300)
ggplot()+
  layer_spatial(hlc,aes(fill=as.factor(class)),colour=NA)+
  scale_fill_viridis_d(name="Character\nTypes",label=c("Modern","Post-\nMedieval","Medieval"))+
  theme_light()+
  annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                         pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                         style=north_arrow_fancy_orienteering())
dev.off()


# Plot LICD - TIFF + JPEG
hlc.id$Type[hlc.id$Type=="No cluster"]<- NA
tiff("Output/Torridge_LICD.tiff",width=14,height=10,units="cm",res=300)
ggplot()+
  layer_spatial(hlc.id,aes(fill=Type),colour=NA)+
  scale_fill_viridis_d(name="LICD",na.value="#DCDCDC",
                       label=c("Clump","Cluster",
                               "Dispersed","Outlier","Outlier \nin dispersion area","No cluster"))+
  annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                         pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                         style=north_arrow_fancy_orienteering())+
  theme_light()
dev.off()

jpeg("Output/Torridge_LICD.jpeg",width=14,height=10,units="cm",res=300)
ggplot()+
  layer_spatial(hlc.id,aes(fill=Type),colour=NA)+
  scale_fill_viridis_d(name="LICD",na.value="#DCDCDC",
                       label=c("Clump","Cluster",
                               "Dispersed","Outlier","Outlier \nin dispersion area","No cluster"))+
  annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                         pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                         style=north_arrow_fancy_orienteering())+
  theme_light()
dev.off()


# HLC & LICD

HLCLICD<-ifelse(hlc.id$class==1 & hlc.id$Type=="Clump","Clump Modern",
                ifelse(hlc.id$class==1 & hlc.id$Type=="Cluster", "Cluster Modern",
                       ifelse(hlc.id$class==1 & hlc.id$Type=="Dispersed","Dispersed Modern",
                              ifelse(hlc.id$class==1 & hlc.id$Type=="Outlier","Outlier Modern",
                                     ifelse(hlc.id$class==1 & hlc.id$Type=="Outlier in dispersion area","Outlier Dispersed Modern",
                                            ifelse(hlc.id$class==2 & hlc.id$Type=="Clump","Clump Post-Med.",
                                                   ifelse(hlc.id$class==2 & hlc.id$Type=="Cluster","Cluster Post-Med.",
                                                          ifelse(hlc.id$class==2 & hlc.id$Type=="Dispersed","Dispersed Post-Med.",
                                                                 ifelse(hlc.id$class==2 & hlc.id$Type=="Outlier","Outlier Post-Med.",
                                                                        ifelse(hlc.id$class==2 & hlc.id$Type=="Outlier in dispersion area","Outlier Dispersed Post-Med.",
                                                                               ifelse(hlc.id$class==3 & hlc.id$Type=="Clump","Clump Medieval",
                                                                                      ifelse(hlc.id$class==3 & hlc.id$Type=="Cluster","Cluster Medieval",
                                                                                             ifelse(hlc.id$class==3 & hlc.id$Type=="Dispersed","Dispersed Medieval",
                                                                                                    ifelse(hlc.id$class==3 & hlc.id$Type=="Outlier","Outlier Medieval",
                                                                                                           ifelse(hlc.id$class==3 & hlc.id$Type=="Outlier in dispersion area","Outlier Dispersed Medieval",NA)))))))))))))))

tiff("Output/Torridge_HLC_LICD.tiff",width=14,height=10,units="cm",res=300)
ggplot()+
  layer_spatial(hlc.id,aes(fill=HLCLICD),colour=NA)+
  scale_fill_viridis_d(direction=-1,na.value="#DCDCDC",name="Classes + LICD")+
  annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                         pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                         style=north_arrow_fancy_orienteering())+
  theme_light()
dev.off()

jpeg("Output/Torridge_HLC_LICD.jpeg",width=14,height=10,units="cm",res=300)
ggplot()+
  layer_spatial(hlc.id,aes(fill=HLCLICD),colour=NA)+
  scale_fill_viridis_d(direction=-1,na.value="#DCDCDC",name="Classes + LICD")+
  annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                         pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                         style=north_arrow_fancy_orienteering())+
  theme_light()
dev.off()

# LICD divided by classes

Mod<-ifelse(hlc.id$class==1 & hlc.id$Type=="Clump","Clump",
            ifelse(hlc.id$class==1 & hlc.id$Type=="Cluster","Cluster",
                   ifelse(hlc.id$class==1 & hlc.id$Type=="Dispersed","Dispersed",
                          ifelse(hlc.id$class==1 & hlc.id$Type=="Outlier","Outlier",
                                 ifelse(hlc.id$class==1 & hlc.id$Type=="Outlier in dispersion area","Outlier Dispersed",NA)))))
PostMed<-ifelse(hlc.id$class==2 & hlc.id$Type=="Clump","Clump",
                ifelse(hlc.id$class==2 & hlc.id$Type=="Cluster","Cluster",
                       ifelse(hlc.id$class==2 & hlc.id$Type=="Dispersed","Dispersed",
                              ifelse(hlc.id$class==2 & hlc.id$Type=="Outlier","Outlier",
                                     ifelse(hlc.id$class==2 & hlc.id$Type=="Outlier in dispersion area","Outlier Dispersed",NA)))))
Med<-ifelse(hlc.id$class==3 & hlc.id$Type=="Clump","Clump",
                ifelse(hlc.id$class==3 & hlc.id$Type=="Cluster","Cluster",
                       ifelse(hlc.id$class==3 & hlc.id$Type=="Dispersed","Dispersed",NA)))


tiff("Output/HLC_Class_LICD.tiff",width=28,height=20,units="cm",res=300)
A<-ggplotGrob(ggplot()+
                layer_spatial(hlc.id,aes(fill=Mod),colour=NA)+
                scale_fill_viridis_d(direction=-1,na.value="#DCDCDC",name="Modern")+
                annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                                       pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                                       style=north_arrow_fancy_orienteering())+
                theme_light())
B<-ggplotGrob(ggplot()+
                layer_spatial(hlc.id,aes(fill=PostMed),colour=NA)+
                scale_fill_viridis_d(direction=-1,na.value="#DCDCDC",name="Post-Medieval")+
                annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                                       pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                                       style=north_arrow_fancy_orienteering())+
                theme_light())
C<-ggplotGrob(ggplot()+
                layer_spatial(hlc.id,aes(fill=Med),colour=NA)+
                scale_fill_viridis_d(direction=-1,na.value="#DCDCDC",name="Medieval")+
                annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                                       pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                                       style=north_arrow_fancy_orienteering())+
                theme_light())
lay<-rbind(c(1,1,2,2),c(0,3,3,0))
grid.arrange(A,B,C, layout_matrix=lay)
dev.off()

jpeg("Output/HLC_Class_LICD.jpeg",width=28,height=20,units="cm",res=300)
A<-ggplotGrob(ggplot()+
                layer_spatial(hlc.id,aes(fill=Mod),colour=NA)+
                scale_fill_viridis_d(direction=-1,na.value="#DCDCDC",name="Modern")+
                annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                                       pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                                       style=north_arrow_fancy_orienteering())+
                theme_light())
B<-ggplotGrob(ggplot()+
                layer_spatial(hlc.id,aes(fill=PostMed),colour=NA)+
                scale_fill_viridis_d(direction=-1,na.value="#DCDCDC",name="Post-Medieval")+
                annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                                       pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                                       style=north_arrow_fancy_orienteering())+
                theme_light())
C<-ggplotGrob(ggplot()+
                layer_spatial(hlc.id,aes(fill=Med),colour=NA)+
                scale_fill_viridis_d(direction=-1,na.value="#DCDCDC",name="Medieval")+
                annotation_north_arrow(location="br",height=unit(0.9,"cm"),width=unit(0.6,"cm"),
                                       pad_x=unit(0.5,"cm"),pad_y=unit(0.3,"cm"),
                                       style=north_arrow_fancy_orienteering())+
                theme_light())
lay<-rbind(c(1,1,2,2),c(0,3,3,0))
grid.arrange(A,B,C, layout_matrix=lay)
dev.off()
