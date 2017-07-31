## Using t-SNE algorithm to classify cluster
install.packages("Rtsne")
install.packages("vegan")
library(vegan)
library(Rtsne)
library(ape)
library(phangorn)

train <- read.table("binary_matrix_PhaseCalling_BrainB.txt", sep="\t")

j_train <- vegdist(train, binary=TRUE, method="jaccard")
labels <-as.factor(gsub("neuron_MDA_","",rownames(train)))
colors=rainbow(length(unique(train$label)))
names(colors) = unique(train$label)

tsne <- Rtsne(j_train, is.distance=TRUE,dims=2, perplexity=5, verbose=TRUE, max_iter=5000)
d_tsne_1 = as.data.frame(tsne$Y)

# k-means
d_tsne_1_original = d_tsne_1
fit_cluster_kmeans=kmeans(scale(d_tsne_1),4)
d_tsne_1_original$cl_kmeans=factor(fit_cluster_kmeans$cluster)

# Hierarchical cluster model
tre <- nj(j_train)
myBoots <- boot.phylo(tre, train, function(x) nj(vegdist(x, binary=TRUE, method="jaccard")), B=1000)
temp<-tre
N <- length(tre$tip.label)
toCollapse <- match(which(myBoots< 700)+N, temp$edge[,2])
temp$edge.length[toCollapse] <- 0
tre2 <- di2multi(temp, tol=0.00001)
fit_cluster_herarchical=as.hclust.phylo(tre2)
d_tsne_1_original$cl_hierarchical=factor(cutree(fit_cluster_herarchical, k=4))

## plotting
plot_cluster=function(data, var_cluster, palette){
     ggplot(data, aes_string(x="V1", y="V2", color=var_cluster))+
     geom_point(size=2) +
     geom_text(aes(label=train$label), hjust=0, vjust=0)+
     guides(color=guide_legend(override.aes=list(size=5)))+
     xlab("") + ylab("") +
     ggtitle("") +
     theme_light(base_size=12) +
     theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction="horizontal",
          legend.position="bottom",
          legend.box="horizontal")+
     scale_colour_brewer(palette=palette)
}

#plot_k=plot_cluster(d_tsne_1_original,"cl_kmeans", "Accent")
plot_h=plot_cluster(d_tsne_1_original, "cl_hierarchical", "Set1")
#gridExtra lib...
library(gridExtra)
pdf("t-sne_clustering_kmeans_hclust.pdf")
grid.arrange( plot_h, ncol=1)
dev.off()




# pdf("t-sne_clustering.pdf")
# plot(tsne$Y, t='n', main='tsne', xlab="tsne-1", ylab="tsne-2")
# text(tsne$Y, labels=train$label, col=colors[train$label])
# dev.off()
