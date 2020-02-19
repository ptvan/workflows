# this code borrows from "Orchestrating Single-Cell Analysis with Bioconductor"
# by Rob Amezquita et al. (https://osca.bioconductor.org/
# formally published in https://www.nature.com/articles/s41592-019-0654-x)

library(BiocManager)
install(c("SingleCellExperiment","scater","scran","uwot","Rtnse", "scRNASeq"))

library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)
library(scRNASeq)
        
#################################################
# GENERATE DUMMY DATA TO WORK WITH THE sce CLASS
#################################################
counts_matrix <- data.frame(cell_1 = rpois(10, 10), 
                            cell_2 = rpois(10, 10), 
                            cell_3 = rpois(10, 30))
rownames(counts_matrix) <- paste0("gene_", 1:10)
counts_matrix <- as.matrix(counts_matrix)

sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
counts(sce)

sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)
logcounts(sce)
assays(sce)

cell_metadata <- data.frame(batch = c(1, 1, 2))
rownames(cell_metadata) <- paste0("cell_", 1:3)

rowRanges(sce)

sce <- addPerFeatureQC(sce)
sce[c("gene_1", "gene_4"), ]

sizeFactors(sce)

#########################
# WORKING WITH REAL DATA
#########################

sce <- MacoskoRetinaData()

library(scater)
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
sce <- sce[, !filtered$discard]

# Normalization.
sce <- logNormCounts(sce)

# Feature selection.
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)

# Dimensionality reduction.
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row=hvg)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

# Clustering.
g <- buildSNNGraph(sce, use.dimred = 'PCA')
sce$clusters <- factor(igraph::cluster_louvain(g)$membership)

# Visualization.
plotUMAP(sce, colour_by="clusters")