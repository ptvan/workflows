library(Biobase)
library(ggplot2)
library(limma)
library(edgeR)
library(DESeq2)
library(GSEABase)
library(dplyr)     
library(stringr)  
library(biomaRt)
library(org.Hs.eg.db)
library(topGO)

######################
# GENERATE DUMMY DATA
######################

# get a bunch of human gene names
genes <- unlist(lookUp(as.character(1:50000), 'org.Hs.eg', 'SYMBOL')) 
genes <- sample(unique(genes[!is.na(genes)]))
names(genes) <- NULL

# generate some dummy subjects, p001 to p050
ptids <- paste0("p",str_pad(as.character(1:50), 3, "left", "0"))

# generate dummy RNASeq counts
counts <- matrix(sample(1:1e6, length(genes)*length(ptids)), nrow=length(genes), ncol=length(ptids)*2)

# each subject has "MEDIA" and "STIM" samples
samples <- c(paste0(ptids,"_MEDIA"), paste0(ptids,"_STIM"))
colnames(counts) <- samples
rownames(counts) <- genes


# create ExpressionSet
dat <- ExpressionSet(assayData=as.matrix(counts))

# make up metadata/covariates
# NOTE: this is a great way to show confounders/batch effects !!!
anno <- as.data.frame(cbind(samples, rep(c("M","F"), length(samples)/2)))
colnames(anno) <- c("sample", "gender")
anno$ptid <- gsub("_MEDIA|_STIM", "", anno$sample)
anno$stim <- gsub("p[0-9]{3}_", "", anno$sample)

# at least *some* of our covariates aren't completely confounded
anno$age <- sample(c(5:60), 100, replace=TRUE)
anno$vaccStatus <- sample(rep(c("Y","N","N","Y","UNKNOWN"), nrow(anno)/5))
anno$school <- factor(rep(c("Springfield","Hogwarts","Xavier","Xavier","Springfield"),2))

# !!! IMPORTANT: make sure the samples' metadata are in the same order
# as the columns of the count matrix !!! This is especially critical
# if there has been merging of the annotations, which reorders the samples silently
anno <- anno[order(anno$sample),] 
counts <- counts[,order(colnames(counts))] 

# for coloring MDS plots later
colvac <- anno$vaccStatus

# cut() will discretize and create factors for us
colage <- cut(anno$age, breaks=c(0,10,30,60), laels=c("children","teens","middleage", "elderly"))

##############################
# MODEL MATRIX & NORMALIZATION
##############################

# create model matrices for different models and corresponding labels for contrasts
mmatrix_stim <- model.matrix(~stim, data = anno)
mmatrix_full <- model.matrix(~stim+school, data=anno)
cons <- c("stim", "school")

# normalize using voom
normy <- calcNormFactors(dat)
libNorm <- colSums(exprs(dat))*normy
data_voomed <- voom(exprs(dat), design=mmatrix_full, plot=FALSE, lib.size=libNorm)
ranCor <- duplicateCorrelation(v, design=mmatrix_full, block=anno$subject)$consensus.correlation

###############################
# QUERYING GENE INFORMATION
###############################

# look up gene symbols' chromosomal location on ENSEMBL
# explicitly set which server we use, since mirrors can go down for maintenance
ensemblMart <- useMart("ensembl"
              ,host = "www.ensembl.org"
              ,ensemblRedirect = FALSE)
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensemblMart)

chrom <- c(1:22,"X","Y")
attrs <- c("hgnc_symbol", "external_gene_name", "chromosome_name", "gene_biotype")


# NOTE: this query can be a bit slow, so likely we'll want to save results
# rather than running it everytime
chrMap <- getBM(attributes=attrs
                ,filters=c("chromosome_name")
                ,values=list(chrom)
                ,mart=ensembl)

# get only protein coding genes
proteinCodingGenes <- subset(chrMap, gene_biotype=="protein_coding")$external_gene_name

## alternatively, we can fetch the gene names from a GTF file, eg. from UCSC hg38
library(refGenome)
hg38gtf <- ensemblGenome()
read.gtf(hg38gtf, "UCSC_hg38_genes.gtf"
         , useBasedir = FALSE)

allGenes <- unique(hg38gtf@ev$gtf$gene_name)

exon <- extractFeature(hg38gtf, "exon")
exonicGenes <- unique(exon@ev$gtf$gene_name)

# mapping mouse <-> human genes via ENSEMBL gene_id
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 

# extra orthology columns so we can filter the mouse genes if necessary
attrs <- c("ensembl_gene_id"
           ,"mmusculus_homolog_ensembl_gene"
           ,"mmusculus_homolog_perc_id_r1"
           ,"mmusculus_homolog_orthology_type"
           ,"mmusculus_homolog_subtype"
           ,"mmusculus_homolog_perc_id"
           ,"mmusculus_homolog_associated_gene_name"
) 

mouse <- getBM( attrs
                ,filters="with_mmusculus_homolog"
                ,values =TRUE
                , mart = human
                , bmHeader=FALSE)

# two separate getBM() calls since BioMart doesn't allow
# gene- and transcript-level queries in the same call
attrs2 <- c("ensembl_gene_id"
            ,"hgnc_symbol"
            ,"external_gene_name")

E2GN <- getBM( attrs2
               ,filters="with_mmusculus_homolog"
               ,values =TRUE
               , mart = human
               , bmHeader=FALSE)

map <- merge(mouse, E2GN, by="ensembl_gene_id")
map$mmusculus_homolog_associated_gene_name <- toupper(map$mmusculus_homolog_associated_gene_name)
map[c("hgnc_symbol","mmusculus_homolog_associated_gene_name")]

##############################
# FILTERING GENES
##############################

# these are the genes actually in our data matrix
expressedGenes <- rownames(v$E)

# get genes on Y chromosome, remove them from data matrix ...
chrYgenes <- intersect(expressedGenes, subset(chrMap, chromosome_name=="Y")$external_gene_name)
eDatnoY <- ExpressionSet(assayData=as.matrix(counts[!rownames(counts) %in% chrYgenes,]))

###############################
#  M.D.S. PLOTTING
###############################

# ... and see how the data separates without sex genes
# plotMDS can either take a voom-transformed structure like `vDatnoY` below...
plotMDS(vDatnoY, col=rainbow(length(unique(colvac)))[as.numeric(colvac)],  main="Vaccinaition status, chrY genes removed", pch=1)
legend("bottomleft", c("yes","no", "unknown"), pch=1, col=c("#FF0000FF", "#00FFFFFF", "#0000FFFF"))

# ... or just a bare expression set, `eDatnoY`
plotMDS(eDatnoY, col=rainbow(length(unique(colage)))[as.numeric(colage)])


################################################## 
# DIFFERENTIALLY-EXPRESSED GENES (D.E.G.) ANALYSIS
##################################################

# fit different models
fit1 <- lmFit(data_voomed, mmatrix_stim, block=anno$ptid, correlation=ranCor)
fit2 <- lmFit(data_voomed, mmatrix_full, block=anno$ptid, correlation=ranCor)
fit2 <- eBayes(fit2, trend=FALSE)

# generate the table of top DEGs from the linear model fit
# for contrasts with 2 factors, sort by t-values
# for contrasts with >2 factors, sort by F-values

allOut <- list()
allOut[[1]] <- topTable(fit2, number=nrow(data_voomed), coef="stimSTIM", sort="P")
allOut[[2]] <- topTable(fit2, number=nrow(data_voomed), coef=grep("school",colnames(fit2)), sort="F")

names(allOut) <- cons


################################################## 
# GENE SET ENRICHMENT ANALYSIS (G.S.E.A)
##################################################

# load GMT files (http://software.broadinstitute.org/gsea/msigdb/collections.jsp)
# GMT files can be concatenated together
FDRCut <- 0.2

gmtFile <- "c2c7.concatenated.v6.0.symbols.gmt"
minGeneSetCut <- 5

geneSet <- getGmt(gmtFile)
geneIds <- geneIds(geneSet)
setsIndices <- ids2indices(geneIds, rownames(vDat$E))
setsIndices <- setsIndices[sapply(setsIndices, length) > minGeneSetCut]

combo <- comboTab <- list()
comboGSEA <- list()
for(i in 1:length(cons)) {
  res <- camera(vDat, setsIndices, design=designMat, contrast=cons[i], sort=TRUE)
  indo <- res$FDR <= FDRCut
  if(sum(indo) > 0) {
    combo[[labs[i]]] <- res[indo, c("Direction","PValue","FDR")]
    comboTab[[labs[[i]]]] <- datatable(res[indo, c("Direction","PValue","FDR")], caption=labs[i]
                                       , extensions = 'Buttons'
                                       , options = list(dom = 'Bfrtip',
                                                        buttons = c('csv', 'excel'))
    )
  }
  comboGSEA[[i]] <- res
}

###################################
# GENE ONTOLOGY (G.O.) ANALYSIS
###################################

# using TopGO
# which in turn requires annotation from org.Hs.eg.db
# read in a table of genes, logFC, AveExpr, t, P.Value, adj.P.Val and B 
# from limma's topTable()
favGenes <- readRDS("DEG_FDR020.Rds")

# topGO expects an named vector of p.values, and a function to select top genes
# to create a GOdata object
GOinput <- favGenes$P.Value
names(GOinput) <- rownames(favGenes)
selection <- function(allScore){ return(allScore < 0.05)}
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")

GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=GOinput,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection,
              nodeSize=10)

# run test on on the GOdata object, tests are listed by topGO:whichTests()
# which includes Fisher, Kolmogoov-Smirnov (KS) and others
# these tests are for the GO enrichment, and *not* for differential expression 
GOoutput <- runTest(GOdata
                    , algorithm = "classic"
                    , statistic = "ks")

GOtable <- GenTable(GOdata, KS=GOoutput, orderBy="KS", topNodes=20)
