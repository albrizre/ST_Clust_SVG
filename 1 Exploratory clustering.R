library(sf)
library(spdep)
library(tmap)
library(sparklyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(mgcv)
library(cluster)
library(factoextra)
library(RColorBrewer)
library(scuttle)
library(BiocFileCache) 
library(data.table)
library(reshape)
library(zCompositions)
library(gridExtra)
library(grid)
options(scipen=10000)

# Set your working directory
setwd("G:/Mi unidad/Investigacion/Milano/SVG density/Github/")
setwd("...")

# Brain data loading
# brain <- LoadData("stxBrain", type = "anterior1")
load("Data/brain.rda")

# Extract counts and coordinates
matrix <- as.matrix(brain@assays$Spatial@counts)
coordinates <- cbind(brain@images$anterior1@coordinates$row,brain@images$anterior1@coordinates$col)

# We select the 1000 genes with greater variability
matrix<-matrix[names(sort(apply(matrix,1,sd),decreasing=TRUE)[1:1000]),]
dim(matrix)
matrix_sample<-data.frame(coordinates,total_counts = colSums(matrix))
rownames(matrix_sample)<-colnames(matrix)
sce <- SingleCellExperiment(as.matrix(matrix))
sce@assays@data@listData$counts <- as.matrix(unlist(sce@assays@data@listData),
                                             nrow = nrow(matrix))
sce@assays@data@listData[["counts"]] <- sce@assays@data@listData[[1]]
# Data normalization (logNormCounts)
sce <- logNormCounts(sce)

dim(sce@assays@data$counts)
dim(sce@assays@data$logcounts)

# Imputation of zeros
imput_0 <- zCompositions::cmultRepl(sce@assays@data$logcounts,z.warning = 1)
# save(imput_0,file = "imput_0.rda") # (cmultRepl can be computationally intensive, better to save)

# Obtain density values and transform to clr
density_values <- imput_0
density_values_clr <- t(apply(density_values,1,function(r) log(r)-mean(log(r))))
dim(density_values_clr)
save(density_values_clr, file = "density_values_clr.rda")

# k-means clustering of genes according to logcounts and density values
logcounts <- sce@assays@data$logcounts
km.logcounts<-list()
for (i in 2:10){
  print(i)
  km.logcounts[[i]]<-kmeans(logcounts,i,nstart = 30)
}
withinss.logcounts<-c()
for (i in 2:10){
  withinss.logcounts[i]<-km.logcounts[[i]]$tot.withinss
}
plot(2:10,withinss.logcounts[2:10],type = "b")

km.dens<-list()
for (i in 2:10){
  print(i)
  km.dens[[i]]<-kmeans(density_values,i,nstart = 30)
}
withinss.dens<-c()
for (i in 2:10){
  withinss.dens[i]<-km.dens[[i]]$tot.withinss
}
plot(2:10,withinss.dens[2:10],type="b")

# We choose k = 5 for comparison

# Clusters based on logcounts
table(km.logcounts[[5]]$cluster)
logcounts_clusters=data.frame(x=coordinates[,1],y=coordinates[,2],
                        c1=apply(logcounts[km.logcounts[[5]]$cluster==1,],2,mean),
                        c2=apply(logcounts[km.logcounts[[5]]$cluster==2,],2,mean),
                        c3=apply(logcounts[km.logcounts[[5]]$cluster==3,],2,mean),
                        c4=apply(logcounts[km.logcounts[[5]]$cluster==4,],2,mean),
                        c5=apply(logcounts[km.logcounts[[5]]$cluster==5,],2,mean)) 
logcounts_clusters_long <- reshape::melt(logcounts_clusters, id = c("x","y"), variable_name = "Cluster")
logcounts_clusters_long$Cluster=as.character(logcounts_clusters_long$Cluster)
logcounts_clusters_long$Cluster[logcounts_clusters_long$Cluster=="c1"]="Cluster 1"
logcounts_clusters_long$Cluster[logcounts_clusters_long$Cluster=="c2"]="Cluster 2"
logcounts_clusters_long$Cluster[logcounts_clusters_long$Cluster=="c3"]="Cluster 3"
logcounts_clusters_long$Cluster[logcounts_clusters_long$Cluster=="c4"]="Cluster 4"
logcounts_clusters_long$Cluster[logcounts_clusters_long$Cluster=="c5"]="Cluster 5"

plots=list()
for (c in 1:5){
  max_aux=max(logcounts_clusters_long[logcounts_clusters_long$Cluster==paste0("Cluster ",c),"value"])
  max_round=round(max_aux,0)
  vector_aux=c(0.5,ifelse(max_aux<max_round,max_round-0.5,max_round))
  plots[[c]]=ggplot(data=logcounts_clusters_long[logcounts_clusters_long$Cluster==paste0("Cluster ",c),],aes(x=y,y=-x,col=value))+
    geom_point(size=0.5)+
    scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")),name="Value")+
    theme_bw()+
    xlab(NULL)+ylab(NULL)+
    theme(text=element_text(size=14),
          plot.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    ggtitle(paste0("Cluster ",c))
}
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
             plots[[5]],
             nrow = 3,top=textGrob('Gene clusters for the log-normalized count data', gp = gpar(fontsize = 26, fontface = 'bold', align = "left")))

# Clusters based on density values
table(km.dens[[5]]$cluster)
density_values_clusters=data.frame(x=coordinates[,1],y=coordinates[,2],
                                c1=apply(density_values[km.dens[[5]]$cluster==1,],2,mean),
                                c2=as.numeric(density_values[km.dens[[5]]$cluster==2,]), # apply cannot be used because there is one element only
                                c3=apply(density_values[km.dens[[5]]$cluster==3,],2,mean),
                                c4=apply(density_values[km.dens[[5]]$cluster==4,],2,mean),
                                c5=apply(density_values[km.dens[[5]]$cluster==5,],2,mean)) 
density_values_clusters_long <- reshape::melt(density_values_clusters, id = c("x","y"), variable_name = "Cluster")
density_values_clusters_long$Cluster=as.character(density_values_clusters_long$Cluster)
density_values_clusters_long$Cluster[density_values_clusters_long$Cluster=="c1"]="Cluster 1"
density_values_clusters_long$Cluster[density_values_clusters_long$Cluster=="c2"]="Cluster 2"
density_values_clusters_long$Cluster[density_values_clusters_long$Cluster=="c3"]="Cluster 3"
density_values_clusters_long$Cluster[density_values_clusters_long$Cluster=="c4"]="Cluster 4"
density_values_clusters_long$Cluster[density_values_clusters_long$Cluster=="c5"]="Cluster 5"

plots=list()
for (c in 1:5){
  max_aux=max(density_values_clusters_long[density_values_clusters_long$Cluster==paste0("Cluster ",c),"value"])
  max_round=round(max_aux,0)
  vector_aux=c(0.5,ifelse(max_aux<max_round,max_round-0.5,max_round))
  plots[[c]]=ggplot(data=density_values_clusters_long[density_values_clusters_long$Cluster==paste0("Cluster ",c),],aes(x=y,y=-x,col=value))+
    geom_point(size=0.5)+
    scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")),name="Value")+
    theme_bw()+
    xlab(NULL)+ylab(NULL)+
    theme(text=element_text(size=14),
          plot.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    ggtitle(paste0("Cluster ",c))
}
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
             plots[[5]],
             nrow = 3,top=textGrob('Gene clusters for the density data', gp = gpar(fontsize = 26, fontface = 'bold', align = "left")))
