# this code borrows from "Orchestrating Single-Cell Analysis with Bioconductor"
# by Rob Amezquita et al. (https://osca.bioconductor.org/
# formally published in https://www.nature.com/articles/s41592-019-0654-x)

library(BiocManager)
# install(c("SingleCellExperiment","scater","scran","uwot","Rtnse", "scRNASeq","DropletUtils", "EnsDb.Hsapiens.v86", "org.Mm.eg.db", "SingleR"))

library(BiocFileCache)
library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)
library(scRNAseq)
library(SingleR)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(org.Mm.eg.db)
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
set.seed(100)
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

### doublet detection
## by clustering
dbl.cls <- doubletCluster(sce.pbmc)

## by simulation
dbl.dens <- doubletCells(sce.pbmc, d=ncol(reducedDim(sce.pbmc)))
sce.pbmc$DoubletScore <- log10(dbl.dens+1)
plotTSNE(sce.pbmc, colour_by="DoubletScore")
plotColData(sce.pbmc, x="label", y="DoubletScore", colour_by="label")

### cell detection
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

### find markers / DEGs
## using t-test, matching the cluster labels stored in the SCE object (k = 18 in this case)
markers.pbmc <- findMarkers(sce.pbmc, test.type="t", groups=colLabels(sce.pbmc))

# explore one of the clusters
rownames(markers.pbmc[[1]])
cluster9 <- markers.pbmc[["9"]]
colnames(cluster9.best)

# using Top to filter
cluster9.best <- cluster9[cluster9$Top <=5, ]

# create a heatmap of logFCs
cluster9.best.logFCs <- getMarkerEffects(cluster9.best)
pheatmap(cluster9.best.logFCs, breaks=seq(-5, 5, length.out=101))

# finding instead cluster-specific markers
markers.pbmc.up3 <- findMarkers(sce.pbmc, pval.type="all", direction="up")
cluster9.specific <- markers.pbmc.up3[["9"]]
colnames(cluster9.specific) # no "Top" column

## using Wilcoxon-Mann-Whitney (WMW) test
markers.pbmc.wmw <- findMarkers(sce.pbmc, test="wilcox", direction="up")

# this testing regime returns FDR and AUCs, so we need to change getMarkerEffects() accordingly:
colnames(markers.pbmc.wmw[[9]])
cluster9.wmw <- markers.pbmc.wmw[[9]]
cluster9.wmw.best <- cluster9.wmw[cluster9.wmw$Top <= 5,]
AUCs <- getMarkerEffects(best.set, prefix="AUC")

# similarly plotting a heatmap of AUCs...
pheatmap(AUCs, breaks=seq(0, 1, length.out=21),
         color=viridis::viridis(21))

### cell type annotation
## assign cell labels using a reference
ref <- BlueprintEncodeData()
pred <- SingleR(test=sce.pbmc, ref=ref, labels=ref$label.main)

# seeing annotated cell types
table(pred$labels)
plotScoreHeatmap(pred)

# checking for low quality assignments
sum(is.na(pred$pruned.labels))
plotScoreDistribution(pred)

# comparing assignments with cluster numbers, adding 10 to avoid strong color jumps with just 1 cell
tab <- table(Assigned=pred$pruned.labels, Cluster=colLabels(sce.pbmc))
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))

## assign cell labels using genesets
# load, QC and normalize the Zeisel brain data
sce.zeisel <- ZeiselBrainData()
sce.zeisel <- aggregateAcrossFeatures(sce.zeisel, 
                                      id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))
rowData(sce.zeisel)$Ensembl <- mapIds(org.Mm.eg.db, 
                                      keys=rownames(sce.zeisel), keytype="SYMBOL", column="ENSEMBL")
stats <- perCellQCMetrics(sce.zeisel
                          , subsets=list(Mt=rowData(sce.zeisel)$featureType=="mito"))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", 
                                              "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[,!qc$discard]
clusters <- quickCluster(sce.zeisel)
sce.zeisel <- computeSumFactors(sce.zeisel, cluster=clusters) 
sce.zeisel <- logNormCounts(sce.zeisel)
wilcox.z <- pairwiseWilcox(sce.zeisel, sce.zeisel$level1class, 
                           lfc=1, direction="up")
markers.z <- getTopMarkers(wilcox.z$statistics, wilcox.z$pairs,
                           pairwise=FALSE, n=50)

### Analyzing T-Cell Receptor (TCR) repertoire
## read in TCR data
bfc <- BiocFileCache(ask=FALSE)
tcr.data <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples/cell-vdj/3.1.0",
  "vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_t_filtered_contig_annotations.csv"))
tcr <- read.csv(tcr.data, stringsAsFactors=FALSE)

# merge with PBMC SCE
tra <- tcr[tcr$chain=="TRA",]
sce.pbmc$TRA <- split(data.frame(tra), factor(tra$barcode, sce.pbmc$Barcode))
