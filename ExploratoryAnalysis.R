### This R script is an exploratory analysis of single cell data. We mainly do Quality Control and PCA t-SNE. Here we use four different datasets come from four different tissues but the same techniques.
library(ggplot2)
### Four tissues: Pancreas, Liver, Embryo, Brain.
liver=read.csv("./Liver/GSE81252/GSE81252_norm_counts.csv",header=T)
embryo=read.csv("./Embryo/GSE65481/GSE65481_discrete_counts.csv",row.names=1)
brain = read.csv("./Brain/Fetal Brain Tissue/GSE76381/GSE76381_norm_counts.csv",row.names=1)
pancreas = read.csv("./Pancreas/GSE83139/GSE83139_norm_counts.csv",row.names=1)
### Basic quality control
nonzeros_genepancreas = apply(pancreas,1, function(x){length(which(x >0))})
nonzeros_cellpancreas = apply(pancreas,2, function(x){length(which(x >0))})
nonzeros_genebrain = apply(brain,1, function(x){length(which(x >0))})
nonzeros_cellbrain = apply(brain,2, function(x){length(which(x >0))})
nonzeros_geneembryo = apply(embryo_RPM,1, function(x){length(which(x >0))})
nonzeros_cellembryo = apply(embryo_RPM,2, function(x){length(which(x >0))})
nonzeros_geneliver = apply(liver,1, function(x){length(which(x >0))})
nonzeros_cellliver = apply(liver,2, function(x){length(which(x >0))})

### filtering: left the genes which ahave at least ecpressed in three cells, and
### left the cells which have at least 200 genes expressd
pancreas_f=pancreas[,which(nonzeros_cellpancreas>=100)]
pancreas_f=pancreas_f[which(nonzeros_genepancreas>=3),]
brain_f=brain[,which(nonzeros_cellbrain>=100)]
brain_f=brain_f[which(nonzeros_genebrain>=3),]
embryo_f=embryo_RPM[,which(nonzeros_cellembryo>=100)]
embryo_f=embryo_f[which(nonzeros_geneembryo>=3),]
liver_f=liver[,which(nonzeros_cellliver>=100)]
liver_f=liver_f[which(nonzeros_geneliver>=3),]

### merge the dataset and create the label
### commongenes 12716 genes
merge1= merge(pancreas_f, embryo_f,all = TRUE) ### when merge, shoud have one column named genename
merge2= merge(brain_f, liver_f,all = TRUE)
merge= merge(merge1, merge2,all = TRUE)
common=intersect(rownames(merge1),rownames(merge2))
rownames(merge)=merge[,1]
merge=merge[,-1]
merge[is.na(merge)] <- 0
merge_t=t(merge)
### create labels for the dataset
merge_t=as.data.frame(merge_t)
merge_t$label=c(rep("pancreas",times=635),rep("embryo",times=22),rep("brain",times=4029),rep("liver",times=777))

## quantile normalization
d=merge1
df=d
for(i in 1:ncol(df)){
      s=sample(35220,35220,replace = F)
      d[s,i]=qqnorm(df[s,i],plot.it = F)$x
      }
mergeQN = d
## mean expression per cell density plot
meanQN=apply(mergeQN,2,mean)
meanQN_data=cbind.data.frame(meanQN,label)
ggplot(meanQN_data, aes(meanQN, color=label)) +geom_density(alpha=.5)+ggtitle("Density plot of mean expression level across groups")


### PCA analysis
### PCA analysis
### ### standardlization check this for gene for cell
merge.pca=prcomp(t(merge))
label = c(rep("pancreas",times=635),rep("embryo",times=22),rep("brain",times=4029),rep("liver",times=777))
# visualize in 2D the first two principal components and color by tissue type
sgCol <- rainbow(length(levels(label)))[label]
data=cbind(merge.pca$x[,1],merge.pca$x[,2],merge.pca$x[,3],merge.pca$x[,4],merge.pca$x[,5],
merge.pca$x[,6],merge.pca$x[,7],merge.pca$x[,8],merge.pca$x[,9],merge.pca$x[,10])

colnames(data)=c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")
###plot(merge.pca$x[,1], merge.pca$x[,3], col=sgCol, pch=16, main='PCA')
data=as.data.frame(data)
data$label=label
ggplot(data,aes(x=PCA1,y=PCA2,color=label))+geom_point(size=1)


### tsne plot
library(Rtsne)
d <- stats::dist(t(merge))
set.seed(0) # tsne has some stochastic steps (gradient descent) so need to set random
tsne_out <- Rtsne(d, is_distance=TRUE, perplexity=10, verbose = TRUE)
label = c(rep("pancreas",times=635),rep("embryo",times=22),rep("brain",times=4029),rep("liver",times=777))
# visualize in 2D the first two principal components and color by tissue type
sgCol <- rainbow(length(levels(label)))[label]
y=tsne_out$Y
colnames(y)=c("y1","y2")
y=as.data.frame(y)
library(ggplot2)
ggplot(y, aes(x=y1,y=y2, color=label)) +geom_point(size=1,alpha=0.3)+ggtitle("tSNE plot")


### hirachical clustering
v <- apply(df, 1, var)
vi <- names(sort(v)[1:100])
hc <- hclust(dist(t(merge[vi,])))
# visualize as heatmap
data=apply(df[vi,], as.numeric)
heatmap(data, Rowv=NA, Colv=as.dendrogram(hc), ColSideColors = sgCol, col=colorRampPalette(c('blue', 'white', 'red'))(100))


### iCA plot
y_ica=ics$A
y_ica=as.matrix(y_ica)
colnames(y_ica)=c("y1","y2")
library(ggplot2)
ggplot(y_ica, aes(x=y1,y=y2, color=label)) +geom_point(size=1,alpha=0.3)+ggtitle("ICA plot")



