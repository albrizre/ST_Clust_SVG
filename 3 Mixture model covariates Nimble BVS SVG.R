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
library(nimble)

# Set your working directory
# setwd("...")

# Load nimble code for the model
source("cluster_model_basic_covariates_BVS_SVG.R")

# Load matrix and coordinates
load("matrix_genes.rda")
load("coordinates_spots.rda")

# Load gene information and create binary covariates

load("Data/GO_tbl.rda")
GOs=unique(GO_tbl$goslim_goa_accession)
GOs=gsub(":","",GOs)
GOs=paste0("X_",GOs)
GO_covariates=matrix(0,nrow=nrow(matrix),ncol=length(GOs))
colnames(GO_covariates)=GOs
rownames(GO_covariates)=toupper(rownames(matrix))
aux=rep(1,nrow(matrix))
names(aux)=rownames(GO_covariates)
GO_tbl$goslim_goa_accession=gsub(":","",GO_tbl$goslim_goa_accession)
for (j in 1:ncol(GO_covariates)){
  print(j)
  GO_name=gsub("X_","",colnames(GO_covariates)[j])
  find=which(GO_tbl$goslim_goa_accession==GO_name)
  genes_found=GO_tbl$mgi_symbol[find]
  GO_covariates[which(rownames(GO_covariates)%in%toupper(genes_found)),j]=1
}
GO_covariates=GO_covariates[which(rownames(GO_covariates)%in%toupper(rownames(matrix))),]
print(sum(GO_covariates))

# SFPCA

load("sfpca.rda")
load("sc_d0.rda")

# Nimble mixture model

N_comp = 12
N_genes = nrow(sfpca$scores)
N_covariates = ncol(GO_covariates)

for (N_clusters in c(3)){ # Vary the number of clusters
  
  constants <- list(
    
    N_genes = N_genes,
    N_comp = N_comp,
    N_clusters = N_clusters,
    N_covariates = N_covariates,
    mu_0 = rep(0, N_comp),
    wish_V = diag(0.01, N_comp),
    tau_0 = diag(0.01, N_comp),
    GO_covariates = GO_covariates,
    sc_d0 = sc_d0
    
  )
  
  y <- sfpca$scores[,1:N_comp]
  data <- list(y = y)
  inits <- function() list(tau = diag(1,N_comp),
                           z = sample(1:(N_clusters+1), N_genes, replace = T),
                           mus = matrix(rep(0, N_comp*(N_clusters+1)), ncol = (N_clusters+1)),
                           omega = t(matrix(rep(rep(1/(N_clusters+1),(N_clusters+1)), N_genes), 
                                            ncol = N_genes)),
                           w = rep(1,N_genes),
                           alpha = matrix(0, nrow = N_genes, ncol = (N_clusters+1)),
                           beta = matrix(0, nrow = N_covariates, ncol = (N_clusters+1)),
                           I = matrix(1, nrow = N_covariates, ncol = (N_clusters+1)),
                           sd_alpha = rep(1,N_clusters+1))
  
  code <- cluster_model_basic_covariates_BVS_SVG
  print(Sys.time())
  mcmc.output <- nimbleMCMC(code, data = data, inits = inits, constants = constants,
                            monitors = c("mus", 
                                         "z", 
                                         "tau",
                                         "alpha",
                                         "I",
                                         "beta",
                                         "sd_alpha"), 
                            niter = 30000, nburnin = 15000, nchains = 1, thin = 10,
                            summary = TRUE, WAIC = TRUE)
  print(Sys.time())
  save(mcmc.output,file=paste0("model_",N_clusters,"_BVS_SVG.rda"))
}


