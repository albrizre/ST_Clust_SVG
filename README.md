R codes and files for fitting the models proposed in the paper *Biologically-informed clustering of gene expression density patterns and identification of spatially variable genes*, by Álvaro Briz-Redón and Alessandra Menafoglio, published in Journal of Agricultural, Biological, and Environmental Statistics.

File description:

- 1 Exploratory clustering.R allows performing an exploratory clustering of the genes based on their spatial expression across the brain tissue considered in the study, considering the data in the form of log-normalized counts and density values. This corresponds to Section 2.3 of the paper.
- 2 SFPCA gene densities.R allows computing the simplicial functional principal components (SFPCs) of the analyzed data. This mainly corresponds to the results displayed in Section 4.1 of the paper.
- 3 Mixture model covariates Nimble BVS SVG.R allows fitting the main clustering model proposed in the paper, which takes as input the values of the main SFPCs obtained with the previous code. It contains the code for fitting the model through the <tt>nimble<tt> R package.
- cluster_model_basic_covariates_BVS_SVG.R is the code of the main model in <tt>nimble<tt> .
