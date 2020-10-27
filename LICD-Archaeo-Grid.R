# Set workspace
library(rstudioapi)
cat("\f")
rm(list=ls())
current_path<-getActiveDocumentContext()$path
setwd(dirname(current_path));print(getwd())

# Import data fro CRAN
library(GISTools)
library(archdata)
data("BarmoseI.grid")
data("BarmoseI.pp")

# Convert to spatial polygon df
coordinates(BarmoseI.grid)<-c("East","North")
gridded(BarmoseI.grid)<-T
barmose<-as(BarmoseI.grid,"SpatialPolygonsDataFrame")

# Convert to spatial point df
cl<-BarmoseI.pp[,3]=="10"
cores<-SpatialPoints(data.frame(East=BarmoseI.pp$East[cl],North=BarmoseI.pp$North[cl]))

# count cores in grid, create categorised value
barmose$cores<-as.vector(poly.counts(cores,barmose))
barmose$class<-(barmose$cores>0)+0

# Create neighbours

## Contiguity neighbours l-order, l = 3 - USE ONLY WHEN WINDOW IS LARGER THAN 1ST ORDER!
barmose.nb <- nblag(poly2nb(barmose),2) ##1 order
#barmose.mat <-as(listw2mat(nb2listw(barmose.nb[[1]])), "sparseMatrix")
#2 order
barmose.mat <- as(listw2mat(nb2listw(barmose.nb[[1]],style="B")),"sparseMatrix")+
  as(listw2mat(nb2listw(barmose.nb[[2]],style="B")),"sparseMatrix")
#3 order
#barmose.mat <- as(listw2mat(nb2listw(barmose.nb[[1]],style="B")),"sparseMatrix") + as(listw2mat(nb2listw(barmose.nb[[2]],style="B")),"sparseMatrix") + as(listw2mat(nb2listw(barmose.nb[[3]],style="B")),"sparseMatrix")

## Set column with analysed data in datafile
clm <- 3
adata <- factor(as.data.frame(barmose[,clm])[,1]) #object with "levels" (factor), first column is not $cat

# Join-Count Statistics
## JC for contiguity

jc.barmose<-list()
jc.barmose.p<-list()

for(i in 1:2){
  jc.barmose[[i]]<-joincount.multi(adata,nb2listw(barmose.nb[[i]]))
}

for(i in 1:2){
  jc.barmose.p[[i]]<-pnorm(jc.barmose[[i]][,4],lower.tail=F)
}

## Exporting output

capture.output(cbind(jc.barmose[[1]],pvalue=as.vector(jc.barmose.p[[1]])),
               file="Output/JCS_Barmose_lag1.txt")
capture.output(cbind(jc.barmose[[2]],pvalue=as.vector(jc.barmose.p[[2]])),
               file="Output/JCS_Barmose_lag2.txt")


###########################################
## Boots' LICD (from Boots 2003) ##
###########################################

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
  # 2 - probability of the "type", 3 - "window" size, 4 - P(X>=x), 5 - P(X<=x)
  
  cluster <- foreach (i=1:length(adata), .combine=rbind) %dopar% {
    require(Matrix)
    c1 <- (barmose.mat[i,] %*% ifelse(adata==2,1,0))+ifelse(adata[i]==2,1,0)
    c2 <- p[2]
    c3 <- sum(barmose.mat[i,])+1
    c4 <- sum(dbinom(c1:c3, size=c3, prob=c2))
    c5 <- sum(dbinom(0:c1, size=c3, prob=c2))
    cbind(c1, c2, c3, c4, c5)
  }
  
  cluster[is.nan(cluster)]<- 1
  ## End of routine 1 ##
  
  ### Custer-outlier analysis -> result of local composition ###
  sc <- 1-(1-0.05)^(1/2) #Sidak correction, 0.05 level of significance
  local_comp <- ifelse(cluster[,4]< sc, 1, ifelse(cluster[,5]< sc, 0, -1))
  # 1 for black, 0 for white, -1 - black-white
  
  #### STEP 2: local configuration
  
  ## Routine 2 in parallel computing ##
  ## We built empirical distribution of joincounts on regular lattices##
  ## We don't have Boots data about distributions, nor Tinkler (1977)
  ## We use permutational approach
  JC.pvalue <- foreach (j=1:length(adata), .combine=rbind) %dopar% {
    require(Matrix)
    library(Rlab)
    require(spdep)
    barmose.mat.1 <- barmose.mat[j,] #extracting a row from weights matrix
    ktore <- which(barmose.mat.1!=0, arr.ind = T) #looking for neighbours of j
    barmose.1 <- barmose[c(j,ktore),] #adding j to the list of its neighbours
    barmose.1.B <- nb2listw(poly2nb(barmose.1, queen=FALSE), style="B") #weights matrix for j and its neighbours
    indices <- summary(barmose.1.B)
    adata01 <- barmose.1@data #substraction data related to window
    if (any(adata01[,clm] != adata01[1,clm])){
      #if any unit is different "type" from j, then I calculate JC
      A <- joincount.multi(as.factor(adata01[,clm]), barmose.1.B)
      A[is.nan(A)]<- 1
    } else {
      A <- matrix(0,3,1)
      if (adata01[1,clm] == 1){
        A[,1] <- c(0,indices$S0/2,0)
      } else {
        A[,1] <- c(indices$S0/2,0,0)
      }
    }
    JC.distrib <- cbind(A[2,1],A[1,1],A[3,1])
    
    # Building empirical distribution - using permutations
    for (s in 1:length(c(j,ktore))^2){
      vector_01 <- as.factor(sample(adata01[,clm])) # permutation of original 0-1 vector
      if (any(vector_01 != vector_01[1])){
        B <- joincount.multi(vector_01, barmose.1.B) # JC for permutation
        JC.distrib <- rbind(JC.distrib, cbind(B[2,1],B[1,1],B[3,1]))
      } else { #JC in the case where all are "0" or "1": The permutations are out of sense
        if (vector_01[1] == 1){
          JC.distrib <- rbind(JC.distrib, cbind(indices$S0/2,0,0)) # all connections are BB or WW
        } else {
          JC.distrib <- rbind(JC.distrib, cbind(0,indices$S0/2,0))  
        }
      }
    }
    # empirical distribution
    JC.pvalue <- (1/(1+s))*cbind(length(which(JC.distrib[,1]>=A[2,1])),length(which(JC.distrib[,2]>=A[1,1])),length(which(JC.distrib[,3]>=A[3,1])))
  }
  ## End of routine 2 ##
  
  colnames(JC.pvalue) <- c("1:1X>=x","0:0X>=x", "1:0X>=x")
  
  local_config <- matrix(3,length(adata),1)
  scJC <- 1-(1-0.05)^(1/3) # Sidak correction JC - 3 tests!
  #scJC <- 0.05 #standard
  
  ### Routine 3 Local configuration
  
  local_config <- foreach (j=1:length(adata), .combine=rbind) %dopar% {#for black is 1, for white is 0, for black-white is -1, otherwise -2
    if (min(JC.pvalue[j,])<scJC){
      ifelse(which(JC.pvalue[j,]==min(JC.pvalue[j,]), arr.ind = T)==1,local_config[j]<- 1,
             ifelse(which(JC.pvalue[j,]==min(JC.pvalue[j,]), arr.ind = T)==3,local_config[j]<- -1,
                    local_config[j]<- 0 ))
    } else {
      local_config[j]<- -2}
  }
  ## End of routine 3 ##
  colnames(local_config) <- c("cluster-dispersion")
  stopCluster(cl)
})
### The end of parallel computing block ###

# Combination of local composition and local configuration
Type <- character(length=0)
C <- cbind(local_comp, local_config)
for (i in 1:length(adata)){
  ifelse(C[i,1] == 1 && C[i,2] == 1, Type[i] <- "Hot Clump",
         ifelse(C[i,1] == 1 && (C[i,2] == -2 || C[i,2] == 0), Type[i] <- "Hot only",
                ifelse(C[i,1] == 0 && (C[i,2] == -2 || C[i,2] == 1), Type[i] <- "Cold only",
                       ifelse(C[i,2] == -1, Type[i] <- "Dispersed only",
                              ifelse(C[i,1] == 0 && C[i,2] == 0, Type[i] <- "Cold clump",
                                     ifelse(C[i,1] == -1 && C[i,2] == 1, Type[i] <- "Clump only (H)",
                                            ifelse(C[i,1] == -1 && C[i,2] == 0, Type[i] <- "Clump only (C)",
                                                   Type[i] <- "No cluster")))))))
}

#### STEP 3: Exporting data

## Preparing for shapefile
ID.data <- sapply(slot(barmose, "polygons"), function(x) slot(x, "ID"))
barmose.id <- cbind(barmose,ID.data)
names(barmose.id)[length(names(barmose.id))] <- c("ID_object")

## Merging data and geometry
X <- cbind(as.data.frame(barmose.id@data[["ID_object"]]), as.data.frame(cbind(C, Type)))
names(X)[1] <- c("ID_object")
barmose.LICD <- merge(barmose.id, X, by="ID_object")
names(barmose.LICD)[length(names(barmose.LICD))-1] <- c("local_config")

## Writing data
#writePolyShape(barmose.LICD, "Barmose_LICD")

## Plot Grid - TIFF + JPEG

tiff("Output/Barmose_Grid_Cores.tiff",width=12,height=10,units="cm",res=300)
ggplot()+
  layer_spatial(barmose,aes(fill=log10(Debitage)),colour="white")+
  layer_spatial(cores)+
  scale_fill_viridis_c(name="Debitage\n(log count)")+
  theme_light()
dev.off()

jpeg("Output/Barmose_Grid_Cores.jpeg",width=12,height=10,units="cm",res=300)
ggplot()+
  layer_spatial(barmose,aes(fill=log10(Debitage)),colour="white")+
  layer_spatial(cores)+
  scale_fill_viridis_c(name="Debitage\n(log count)")+
  theme_light()
dev.off()

# Plot results - TIFF + JPEG

tiff("Output/Barmose_Class_LICD.tiff",width=24,height=10,units="cm",res=300)
A<-ggplotGrob(ggplot()+
                layer_spatial(barmose,aes(fill=factor(class)),colour="white")+
                scale_fill_viridis_d(begin=0.2,end=0.7,name="Classes",
                                     labels=c("No cores", "Core(s)"))+
                theme_light())
B<-ggplotGrob(ggplot()+
                layer_spatial(barmose.LICD,aes(fill=Type),colour="white")+
                scale_fill_viridis_d(begin=0.2,end=1,direction=-1,name="LICD")+
                theme_light())
grid.arrange(A,B,layout_matrix=rbind(c(1,2)))
dev.off()

jpeg("Output/Barmose_Class_LICD.jpeg",width=24,height=10,units="cm",res=300)
A<-ggplotGrob(ggplot()+
                layer_spatial(barmose,aes(fill=factor(class)),colour="white")+
                scale_fill_viridis_d(begin=0.2,end=0.7,name="Classes",
                                     labels=c("No cores", "Core(s)"))+
                theme_light())
B<-ggplotGrob(ggplot()+
                layer_spatial(barmose.LICD,aes(fill=Type),colour="white")+
                scale_fill_viridis_d(begin=0.2,end=1,direction=-1,name="LICD")+
                theme_light())
grid.arrange(A,B,layout_matrix=rbind(c(1,2)))
dev.off()

# Plot Cores + LICD - TIFF + JPEG

LICDCore<-ifelse(barmose.LICD$class==0 & barmose.LICD$Type=="Clump only (C)","No Core Clump",
                 ifelse(barmose.LICD$class==1 & barmose.LICD$Type=="Clump only (C)", "Core Clump",
                        ifelse(barmose.LICD$class==0 & barmose.LICD$Type=="Cold only", "No Core Cold",
                               ifelse(barmose.LICD$class==1 & barmose.LICD$Type=="Cold only","Core Cold ",
                                      ifelse(barmose.LICD$class==0 & barmose.LICD$Type=="Hot only","No Core Hot ",
                                             ifelse(barmose.LICD$class==1 & barmose.LICD$Type=="Hot only","Core Hot",
                                                    ifelse(barmose.LICD$class==0 & barmose.LICD$Type=="No cluster","No Core No Cluster ","Core No Cluster")))))))


tiff("Output/Barmose_LICD_Cores.tiff",width=14,height=10,units="cm",res=300)
ggplot()+
  layer_spatial(barmose,aes(fill=LICDCore),colour="white")+
  scale_fill_viridis_d(direction=-1,name="Classes + LICD")+
  theme_light()
dev.off()

jpeg("Output/Barmose_LICD_Cores.jpeg",width=14,height=10,units="cm",res=300)
ggplot()+
  layer_spatial(barmose,aes(fill=LICDCore),colour="white")+
  scale_fill_viridis_d(direction=-1,name="Classes + LICD")+
  theme_light()
dev.off()
