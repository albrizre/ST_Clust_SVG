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

# Set your working directory
# setwd("...")

# Brain data loading
# InstallData("stxBrain") 
brain <- LoadData("stxBrain", type = "anterior1")

# Load matrix and coordinates
load("matrix_genes.rda")
load("coordinates_spots.rda")

# Load clr values
load("density_values_clr.rda")

### Auxiliary functions
inner_product <- function (f,g){
  sum(f*g)
}

clr2density <- function(clr){ # back-transform a clr to a density
  exp(clr)/sum(exp(clr))
}

######################### SIMPLICIAL FPCA ###############################

clr_values <- t(density_values_clr)

# Compute the mean
m1 <- apply(clr_values, 1, mean)
n1 <- dim(clr_values)[2]

# Numerical way to get SFPCs
lmod <- clr_values
lmod <- lmod
SS <- (n1-1)/n1*cov(t(lmod)) # we can divide by n or n-1
S <- SS
pc <- eigen(S)

# Set value of K
K <- 100

# Check norm equal to 1
norms <- c()
for (k in 1:K){
  norms <- c(norms,inner_product(pc$vectors[,k],pc$vectors[,k]))
}

# Compute the centred observations
clr_values.c <- NULL
for(i in 1:n1)
  clr_values.c <- cbind(clr_values.c, clr_values[,i]-m1)

# Compute the scores
sc <- matrix(NA, ncol = K, nrow = n1)
for(i in 1:n1)
{
  for(j in 1:K)
    sc[i,j] <- inner_product(clr_values.c[,i],pc$vec[,j])
}

# Choose the truncation order K: look at the boxplots of the scores and the scree plot 
boxplot(sc, las = 1, col = 'gold', main = 'Principal Components')
plot(1:K, cumsum(pc$val)[1:K]/sum(pc$val), type = 'b', pch = 20, ylim = c(0,1),
     xlab = 'K', ylab = '% Variance')
abline(h = .95, lty = 2, col = 'grey')
abline(h = .98, lty = 2, col = 'grey')

plot(1:K, pc$values[1:K], type = 'b', pch = 20, 
     xlab = 'K', ylab = 'Eigenvalue')
abline(h = 1, lty = 2, col = 'red')
which.min(pc$values>1)

# Checking clr and SFPCA properties
dim(density_values_clr)
j=1;inner_product(density_values_clr[j,],1) # clr integrate = 0
j=1;inner_product(pc$vectors[,j],pc$vectors[,j]) # norm = 1
j1=1;j2=2;inner_product(pc$vectors[,j1],pc$vectors[,j2]) # orthogonality
j=1;inner_product(pc$vectors[,j],1) # zero-integral

# Projected densities  
clr_values.p <- NULL
for(i in 1:n1)
{
  tmp <- m1
  for(j in 1:K)
    tmp <- tmp+sc[i,j]*(pc$vec[,j])
  clr_values.p <- cbind(clr_values.p, (tmp))
}
d.p <- NULL
for(i in 1:n1)
{
  d.p <- cbind(d.p,clr2density(clr_values.p[,i]))
}

# Plot check
data_gen_aux <- data.frame(x=coordinates[,1],y=coordinates[,2],
                          clr_values.p=as.numeric(clr_values.p[,2]),
                          d.p=as.numeric(d.p[,2]))
ggplot(data=data_gen_aux,aes(x=y,y=-x,col=d.p))+
  geom_point(size=1.5)+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")),name = "d.p")
ggplot(data=data_gen_aux,aes(x=y,y=-x,col=clr_values.p))+
  geom_point(size=1.5)+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")),name = "clr_values.p")

# SFPCA object
sfpca <- list()
sfpca$values <- pc$values
sfpca$harmonics <- pc$vectors
sfpca$scores <- sc

# Plot SFPCs
data_gen_aux <- c()
for (k in 1:12){
  data_gen_aux <- rbind(data_gen_aux,data.frame(x=coordinates[,1],
                          y=coordinates[,2],
                          SFPCA=pc$vec[,k],
                          Comp=paste0("SFPC ",k)))
}
data_gen_aux$Comp <- factor(data_gen_aux$Comp,levels=paste0("SFPC ",1:12))
ggplot(data=data_gen_aux,aes(x=y,y=-x,col=SFPCA))+
  facet_wrap(~Comp)+
  geom_point(size=0.2)+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")),name="Value")+
  theme_bw()+
  xlab(NULL)+ylab(NULL)+
  ggtitle("Maps of the first 12 SFPCs")+
  theme(text=element_text(size=26),
        plot.title = element_text(face="bold"),
        legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Plot percentage of variability and distribution of the scores
df_perc=data.frame(K=1:100,CumVar=100*cumsum(pc$val)[1:K]/sum(pc$val))
ggplot(data=df_perc,aes(x=K,y=CumVar))+
  geom_line(size=1,col="gray30")+
  theme_bw()+
  xlab("Number of SFPC")+ylab("Percentage of variance")+
  ggtitle("Variability explained by the first 100 SFPCs")+
  scale_x_continuous(limits=c(1,100),breaks=c(1,seq(11,100,10)-1))+
  scale_y_continuous(limits=c(0,100))+
  theme(text=element_text(size=26),
        plot.title = element_text(face="bold"),
        legend.position = "right")

df_scores=c()
for (i in 1:K){
  df_scores=rbind(df_scores,data.frame(K=i,Score=sc[,i]))
}
ggplot(df_scores, aes(x=as.factor(K), y=Score)) + 
  geom_boxplot(fill="gray80")+
  theme_bw()+
  xlab("Number of SFPC")+ylab("SFPCA scores")+
  scale_x_discrete(breaks=c(1,seq(11,100,10)-1))+
  ggtitle("Distribution of SFPCA scores")+
  theme(text=element_text(size=26),
        plot.title = element_text(face="bold"),
        legend.position = "right")

# Export for model-based clustering
save(sfpca,file="sfpca.rda")

# Scores for the uniform density and export
d0 <- rep(1/ncol(matrix),ncol(matrix))
clr0 <- log(d0)-mean(log(d0))
clr0_c <- clr0-m1
K <- 1000
sc_d0 <- rep(NA,K)
for (j in 1:K){
  sc_d0[j] <- inner_product(clr0_c,pc$vec[,j])
}
plot(sc_d0)
save(sc_d0,file = "sc_d0.rda")
