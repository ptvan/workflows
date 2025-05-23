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
library(magrittr)
library(reshape2)
library(sva)
library(scBatch)
library(splines)

##############################################
# ALIGN READS TO GENOME USING RNASeqPipelineR
##############################################
devtools::install_github("RGLab/RNASeqPipelineR")
library(RNASeqPipelineR)
library(parallel)

PREFIX <- "/home/pvan/data/"

createProject(project_name = "RNASeq_demo" ,path = PREFIX, load_from_immport = FALSE)
loadProject(project_dir=PREFIX, name="RNASeq_demo")

utils_dir <- "/fh/fast/gottardo_r/10_ref_files/Reference_Genome/Homo_sapiens/UCSC/hg38/"
setGenomeReference(utils_dir, "hg38.fa")

buildReference(path=utils_dir
               , gtf_file="UCSCDec2016.gtf"
               , fasta_file="hg38.fa"
               , isoformsFile="UCSCKnownIsoformsDec2016.txt"
               , doSTAR=TRUE
               , threads=6
               , name="hg38")

# note: RNASeqPipelineR also supports Kallisto for alignment
# which is faster and more accurate than STAR
# below are options for paired reads

AlignmentSTAR(parallel_threads=4
              , star_threads=2
              , paired=TRUE
              , paired_pattern=c("_R1.fq", "_R2.fq"))

setGenomeReference(utils_dir, "hg38")

# nchunks should be multiple of parallel_threads
RSEMCalculateExpression(parallel_threads=8,bowtie_threads=1,paired=TRUE, nchunks=16,
                        slurm=FALSE, fromBAM=TRUE, fromSTAR=TRUE)
RSEMAssembleExpressionMatrix(force=TRUE)

annotateUCSC(genome="hg38", force=TRUE)

anno <- read.csv("/home/pvan/data/metadata.csv")

runFastQC(ncores=2)
qc_matrix <- QualityControl(paired=TRUE)
mData <- mergeData(mergeAnnotation=FALSE)

counts <- sumDuplicates(mData$counts, mData$featureData, NULL)



#########################################################################
# ALTERNATIVELY, USE kallisto` or `salmon` TO QUANTIFY WITHOUT ALIGNING
#########################################################################
## !!! IMPORTANT: PSEUDOALIGNERS WON'T DETECT NOVEL TRANSCRIPTS OUTSIDE INDICES !!!
## download kallisto transriptome indices from Pachter lab GitHub, eg.
## wget https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz

# cd ~/working/Databases/homo_sapiens
# kallisto index -i transcriptome.idx transcriptome.fa
# for fn in data/sample{01..40};
#do
# samp=`basename ${fn}`
# echo "Processing sample ${samp}"
# kallisto quant -i transcriptome.idx -o output ${fn}/${samp}_1.fastq.gz  ${fn}/${samp}_2.fastq.gz
# done 

counts <- read.csv("abundance.tsv")

###########################################
# ALTERNATELY, GENERATE DUMMY COUNTS MATRIX
###########################################

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

# for coloring MDS plots later, cut() will discretize and create factors for us

colage <- cut(anno$age
              , breaks=c(0,10,30,60), laels=c("children","teens","middleage", "elderly"))


##############################
# MODEL MATRIX & NORMALIZATION
##############################

# create model matrices for different models and corresponding labels for contrasts
# alternatively, we can use `makeContrasts()`
mmatrix_stim <- model.matrix(~stim, data = anno)
mmatrix_full <- model.matrix(~stim+school, data=anno)
  
cons <- c("stim", "school")

# for a continuous variable like age, we can put the variable directly in
# which assumes linear trend
mmatrix_age <- model.matrix(~age, data=anno)

# assuming non-linear age effect, we can fit a spline
# NOTE: the spline uses up degrees of freedom for the fit
agespline <- ns(anno$age, df=5)
mmatrix_age <- model.matrix(~agespline, data=anno)


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
# plotMDS can either take a voom-transformed structure like `vDatnoY`
plotMDS(vDatnoY, col=rainbow(length(unique(colvac)))[as.numeric(colvac)],  main="Vaccinaition status, chrY genes removed", pch=1)
legend("bottomleft", c("yes","no", "unknown"), pch=1, col=c("#FF0000FF", "#00FFFFFF", "#0000FFFF"))

# ... or just a bare expression set, `eDatnoY`
# either way, we have to supply a vector of colors for each covariate we're facetting  
plotMDS(eDatnoY, col=rainbow(length(unique(colage)))[as.numeric(colage)])

# alternatively, we can just get the coordinates, merge metadata and use ggplot
# which also avoids rerunning plotMDS multiple times
tmp <- plotMDS(eDatnoY)
mds <- data.frame(tmp$x, tmp$y) 
mds$age <- anno$age

ggplot(mds, aes(col=age)) +
  geom_point() +
  labs(title="MDS, colored by age")


###############################
#  BATCH EFFECT ADJUSTMENT
###############################
# using ComBat-Seq
batch <-  c(rep(1, 50), rep(2,50))
counts_adjusted <- ComBat_seq(counts, batch=batch, group=NULL)

# alternatively, using scBatch
counts_adjusted <- QuantNorm(counts
                             , batches
                             , method='row/column'
                             , cor_method='pearson'
                             , logdat=F
                             , standardize = T
                             ,tol=1e-4)

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

######################
# D.E.G. WITH KINSHIP 
######################
# we can fit a familyID covariate using limma, which may not show any DEGs
# alternatively, can explicitly model kinship with a pairwise matrix if available

library(coxme)
library(kinship2)
library(parallel)

data(sample.ped)

# using an example kinship matrix from `kinship2`
pedAll <- pedigree(id = sample.ped$id, dadid = sample.ped$father, momid = sample.ped$mother, 
                   sex = sample.ped$sex, famid = sample.ped$ped)
ped2basic <- pedAll["2"]
kin <- kinship(ped2basic)

# make up corresponding counts data for the subjects in the pedigree
# modeling voomed data would make more sense
counts_k <- matrix(sample(1:1e6, length(genes)*ncol(kin)), nrow=length(genes), ncol=ncol(kin))
rownames(counts_k) <- genes
colnames(counts_k) <- colnames(kin)

# make up corresponding metadata
anno_k <- data.frame(cbind(colnames(kin), c(rep("Y",7),rep("N",7))))
colnames(anno_k) <- c("ptid", "vaccStatus")

out.mat <- data.frame(counts_k[,c(1,2)])
rownames(out.mat) <- genes
colnames(out.mat) <- c("pval","sigma")
out.mat[,"gene"] <- genes

# thanks to Kim Dill-McFarland for original code on methylation data

fit.kin.results <- data.frame()

fit.kin.results <- rbindlist(foreach(i=1:nrow(RSTR.M.kin)) %dopar%{
  print(i)
  
    gene <- rownames(counts_k)[i]
    
    fit <- lmekin(counts_k[i,] ~ vaccStatus + (1|ptid), 
       data=anno_k, varlist=as.matrix(kin))
    beta <- fit$coefficients$fixed
    nvar <- length(beta)
    nfrail <- nrow(fit$var) - nvar
    se <- sqrt(diag(fit1$var)[nfrail + 1:nvar])
    p.kin <- signif(1 - pchisq((beta/se)^2, 1), 2)[2]
    sigma.kin <- fit$sigma
    
    out.mat[gene,]$pval <- p.kin
    out.mat[gene,]$sigma <- sigma.kin
})

# add Benjamini-Hochberg FDR to make results comparable to limma
fit.kin.results$FDR <- p.adjust(fit.kin.results$pval, method="BH")


######################################
#  DIMENSIONAL REDUCTION / CLUSTERING
######################################

# we could reduce the dimensions of the data down to help with modeling 
# using WGCNA
library(WGCNA)
# filter permissively, FDR < 0.3
iGenes03 <- topTable(fit2, number=nrow(data_voomed), coef="stimSTIM", sort="P") %>% dplyr::filter(adj.P.Val < 0.3)

exprFDR03 <- vDat$E[iGenes03$gene,]
# sft <- pickSoftThreshold(t(exprFDR03), powerVector=c(1:20), verbose=5)
# pwr <- sft$powerEstimate
modules <- blockwiseModules(t(exprFDR03), power=8, 
                            networkType="signed", TOMType="signed",
                            minModuleSize=20, maxBlockSize=600, deepSplit=4,
                            numericLabels=TRUE,
                            saveTOMFileBase="TOM-blockwise")
WGCNAmodules <- modules$colors
# skip module 0 (genes that don't fit in any module) 
WGCNAmodules <- WGCNAmodules[WGCNAmodules > 0]


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
