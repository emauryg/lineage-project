## R-code to perform preliminary hierchical clustering analysis
library(ape)
library(phangorn)
library(seqinr)
library(dendextend)
library(adephylo)

distance_matrix <- matrix(0, length(samples), length(samples))
rownames(distance_matrix) <- gsub("1465-cortex_1-","", samples)
colnames(distance_matrix) <- gsub("1465-cortex_1-","", samples)

for ( cell1 in rownames(binary_matrix)){
     for (cell2 in rownames(binary_matrix)){
          d = sum(binary_matrix[cell1,]!=binary_matrix[cell2,] )
          distance_matrix[cell1, cell2]  <- d
     }
}

     d_cells <- dist(distance_matrix)
     hc_cells <- hclust(d_cells, method="ward.D2")
     dend <- as.dendrogram(hc_cells)
     dend <- set(dend, "labels_cex", 0.5)
     #dend <- color_branches(dend, k=4)
     plot(dend, horiz=TRUE, edgePar=list(cex=4))

## trying to do bootstrapping
temp <- t(as.matrix(d_cells))
temp <- temp[,ncol(temp):1]
par(mar=c(1,5,5,1))
image(x=1:16, y=1:16, temp, col=rev(heat.colors(100)), xaxt="n", yaxt="n", xlab="", ylab="")
axis(side=2, at=1:16, lab=rownames(binary_matrix), las=2, cex.axis=0.5)
axis(side=3, at=1:16, lab=rownames(binary_matrix), las=3, cex.axis=0.5)


snvs <- as.phyDat(t(binary_matrix), type="USER", levels=c(0,1))
## Checking NJ
tre <- nj(d_cells)
plot(tre, cex=0.6)
axisPhylo()
title("NJ tree Bulk Phasing")

x <- as.vector(d_cells)
y <- as.vector(as.dist(cophenetic(tre)))
plot(x,y, xlab="original pairwise distances", ylab="pairwise distance on the tree", main="Is NJ appropriate?", pch=20,cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2

## Checking UPGMA
tre2 <- as.phylo(hclust(d_cells, method="average"))
y<- as.vector(as.dist(cophenetic(tre2)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is UPGMA appropriate?", pch=20, cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2

plot(tre2, cex=0.5)
title("UPGMA tree Bulk Phasing")

binary.dist <- function(binary_matrix){
     for ( cell1 in rownames(binary_matrix)){
     for (cell2 in rownames(binary_matrix)){
          d = sum(binary_matrix[cell1,]!=binary_matrix[cell2,] )
          distance_matrix[cell1, cell2]  <- d
     }
     return(dist(distance_matrix))
}

}

## Bootstrapping using the NJ method
myBoots <- boot.phylo( tre, binary_matrix,  function(x) nj(binary.dist(x)), B=10000 )

plot(tre, edge.width=2)
title("NJ tree + bootstrap values")
axisPhylo()
nodelabels(myBoots, cex=0.6)

## Collapse not well supported nodes into multifurcations
temp <- tre
N <- length(tre2$tip.label)
toCollapse <- match(which(myBoots <7000) +N, temp$edge[,2])
temp$edge.length[toCollapse] <- 0
tre3 <- di2multi(temp, tol=0.00001)
plot(tre3, edge.width=2)
title("NJ tree after collapsing weak nodes")
axisPhylo()
