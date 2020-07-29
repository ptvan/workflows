# this code borrows from "Orchestrating Single-Cell Analysis with Bioconductor"
# by Rob Amezquita et al. (https://osca.bioconductor.org/
# formally published in https://www.nature.com/articles/s41592-019-0654-x)

library(BiocManager)
# install(c("SingleCellExperiment","scater","scran","uwot","Rtnse", "scRNASeq","DropletUtils", "EnsDb.Hsapiens.v86"))

library(BiocFileCache)
library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)
library(scRNAseq)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(pheatmap)

        
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

######################################
# MOUSE DATA FROM Macosko et al 2016
######################################

sce <- MacoskoRetinaData()

library(scater)
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
sce <- sce[, !filtered$discard]

### normalization.
sce <- logNormCounts(sce)

### feature selection.
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)

### dimensionality reduction.
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row=hvg)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

### clustering.
g <- buildSNNGraph(sce, use.dimred = 'PCA')
sce$clusters <- factor(igraph::cluster_louvain(g)$membership)

### visualization.
plotUMAP(sce, colour_by="clusters")


########################
# 4K PMBC DATA FROM 10X 
########################
### load data
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir=file.path(tempdir(), "pbmc4k"))
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)

### gene annotation
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID")

### cell detection
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

### QC
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.pbmc <- sce.pbmc[,!high.mito]

### normalization
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

### model variance
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)

### dimensional reduction
sce.pbmc <- denoisePCA(sce.pbmc, subset.row=top.pbmc, technical=dec.pbmc)
sce.pbmc <- runTSNE(sce.pbmc, dimred="PCA")
sce.pbmc <- runUMAP(sce.pbmc, dimred="PCA")

### clustering
g <- buildSNNGraph(sce.pbmc, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.pbmc) <- factor(clust)

### find markers
## using t-test, matching the cluster labels stored in the SCE object (k = 18 in this case)
markers.pbmc <- findMarkers(sce.pbmc, test.type="t", groups=colLabels(sce.pbmc))

# explore one of the clusters
rownames(markers.pbmc[[1]])
cluster9 <- markers.pbmc[["9"]]
cluster9.best <- cluster9[cluster9$Top <=5, ]
cluster9.best.logFCs <- getMarkerEffects(cluster9.best)
pheatmap(cluster9.best.logFCs, breaks=seq(-5, 5, length.out=101))
