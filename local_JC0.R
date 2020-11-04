## early prototype function, both internals and argument names will change

local_JC0 <- function(obj, lagsmat, varname, numvar, p) {

## cluster has 5 columns: 1 - number of units of the "type" as the unit j in the "window" 2 - probability of the "type", 3 - "window" size, 4 - probability of "absence" of "type"

c1 <- c2 <- c3 <- c4 <- c5 <- numeric(length(adata))
for (i in 1:length(adata)) {
    c1[i] <- (lagsmat[i,] %*% ifelse(adata==adata[i], 1, 0)) + 1
    c2[i] <- p[adata[i]]
    c3[i] <- sum(lagsmat[i,])+1
    c4[i] <- sum(dbinom(0:c1[i], c3[i], c2[i])) ###CHANGE!
    c5[i] <- sum(dbinom(c1[i]:c3[i], c3[i], c2[i])) ###CHANGE!
}
cluster_seq <- cbind(c1, c2, c3, c4, c5)
cluster_seq[is.nan(cluster_seq)]<- 1

## Cluster-outlier analysis -> result of local composition ###
## Sidak correction, 0.05 level of significance
sc <- 1-(1-0.05)^(1/2)
local_comp_seq <- ifelse(cluster_seq[,5] < sc, 1, ifelse(cluster_seq[,4] < sc, -1, 0))

## We calculate of bivariate local JC's: 1 for units of the same "type" as j, 0 otherwise

clm <- varname
JC.pvalue_seq <- matrix(0, ncol=3, nrow=length(adata))

for (j in 1:length(adata)) {
    lagsmat.1 <- lagsmat[j,] #extracting a row from weights matrix
    ktore <- which(lagsmat.1 !=0, arr.ind = TRUE) #looking for neighbours of j
    hlc.1 <- obj[c(j,ktore),] #adding j to the list of its neighbours
    hlc.1.B <- mat2listw(lagsmat[c(j,ktore), c(j,ktore)], row.names = NULL, style="B")
    datawork <- as.data.frame(hlc.1) #working data
    adata01 <- matrix(0,length(ktore)+1,1) #"types" of units and j "type"
    
    for (i in 1:length(ktore)+1){
      adata01[i,1] <- ifelse(datawork[i,clm]==datawork[1,clm],1,0)
      #"Type" of unit is 1 if it is as the "type" of j, 0 otherwise
    } 
    adata01[1,1]<-1 
    
    if (any(adata01 == 0)){
      #if any units is different "type" from j, then I calculate JC
      A <- joincount.multi(as.factor(adata01), hlc.1.B)
      JC.pvalue_seq[j,] <- c(1-pnorm(A[2,4]), 1-pnorm(A[1,4]), pnorm(A[3,4])) #p-value JC
    } else {
      #if all are "type" as j, then JC is out of sense
      JC.pvalue_seq[j,] <- c(0, 1, 1)
    }
}
JC.pvalue_seq[is.nan(JC.pvalue_seq)]<- 1
colnames(JC.pvalue_seq) <- c("1:1", "0:0", "1:0")
list(local_comp_seq, JC.pvalue_seq)
}

