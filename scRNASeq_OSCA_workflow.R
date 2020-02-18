# this code borrows from "Orchestrating Single-Cell Analysis with Bioconductor"
# by Rob Amezquita et al. (https://osca.bioconductor.org/
# formally published in https://www.nature.com/articles/s41592-019-0654-x)

library(BiocManager)
install(c("SingleCellExperiment","scater","scran","uwot","Rtnse"))

library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)
        
###########################################
# GENERATE DUMMY DATA FOR REMAINING EXAMPLES
###########################################
counts_matrix <- data.frame(cell_1 = rpois(10, 10), 
                            cell_2 = rpois(10, 10), 
                            cell_3 = rpois(10, 30))
rownames(counts_matrix) <- paste0("gene_", 1:10)
counts_matrix <- as.matrix(counts_matrix)

sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
counts(sce)