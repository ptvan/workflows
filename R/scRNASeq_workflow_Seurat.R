library(Seurat)
library(SeuratData)
library(glmGamPoi)
library(presto)

### Load the 3k PBMC dataset
# data is from https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
pbmc.data <- Read10X(data.dir = "~/working/datasets/scRNASeq/pbmc3k/filtered_gene_bc_matrices/hg19/")

# Seurat counts are stored in Sparse Matrices
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

head(pbmc@meta.data, 5)

# Alternatively, use the processed data from the SeuratData package
InstallData("pbmc3k")
pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")

# Low-quality / dying cells often exhibit extensive mitochondrial contamination
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

### Normalization
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

# find most variable features
top10 <- head(VariableFeatures(pbmc), 10)

# calculate module score
my_module <- c("ASXL1","KIT","NPM1","FLT3")
pbmc <- AddModuleScore(
  pbmc,
  features = list(c(my_module))
)

### General plotting 
# plot variable features with and without labels
p1 <- VariableFeaturePlot(pbmc)
p2 <- LabelPoints(plot = p1, points = top10_manual, repel = TRUE)
p1 + p2

# expressions of genes across different cell subsets
RidgePlot(pbmc, features=c("CCL5","IL32"))

### Dimensional Reduction
## center and scale data 
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

## run PCA, glance quickly at results
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), verbose = FALSE)
print(pbmc[["pca"]])
ElbowPlot(pbmc)

# visualize features
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:2, cells = 500, balanced = TRUE)

## cluster the cells, Satija lab recommends 0.4-1.2 for 3k cells
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5)

## run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

## plot existing DimReduc slot
p3 <- DimPlot(pbmc, reduction="pca")
p4 <- DimPlot(pbmc, reduction="umap")
p3 + p4

### Find differentially expressed features
## list markers for cluster1
cluster1.markers <- FindMarkers(pbmc, ident.1 =1)
head(cluster2.markers, n = 5)

# markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

## DEG conditional to a gene marker (comparing GeneX+ cells vs. GeneX- cells)
## IMPORTANT: the next 2 examples below depend on made-up metadata:
# a `treatment` column with 2 possible values: `control` or `treated`
# CALCR+ vs. CALCR-
Calcr_cells <- WhichCells(pbmc, expression = Calcr > 0, slot = "counts")
pbmc$Calcr_expr <- ifelse(colnames(pbmc) %in% Gfral_cells, "pos", "neg")
DEG_Calcr <- FindMarkers(so
                         , ident.1="pos", ident.2="neg"
                         , group.by = "Calcr_expr"
                         , test.use="wilcox"
                         , assay = "RNA"
                         , slot = "counts"
)

# Gipr+ cells, control vs. treated
pbmc_Gipr <- subset(pbmc, subset = Gipr_expr == "pos", slot = "counts")
DEG_Gipr_fasting <- FindMarkers(so_Gipr
                                , ident.1="control", ident.2="treated"
                                , group.by = "treatment"
                                , test.use="wilcox"
                                , assay = "RNA"
                                , slot = "counts"
)