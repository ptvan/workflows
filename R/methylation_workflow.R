library(ChAMP)
library(ChAMPdata)

testDir <- system.file("extdata", package = "ChAMPdata")
data <- champ.load(testDir, 
                   arraytype = "450K")

# perform QC and save the associated plots
champ.QC(beta = data$beta,
         pheno = data$pd$Sample_Group,
         mdsPlot = TRUE,
         densityPlot = TRUE,
         dendrogram = TRUE,
         PDFplot = TRUE,
         Rplot = TRUE,
         Feature.sel = "None",
         resultsDir="./CHAMP_QCimages/")

# perform normalization
norms <- champ.norm(beta = data$beta,
                    arraytype = "450K",
                    cores = 4)

# spot-check data after normalization
QC.GUI(beta = norms)


# plot SVD to look at covariates
champ.SVD(beta = norms,
          pd = data$pd)

# perform batch correction using ComBat
datCombat <- champ.runCombat(beta = norms,
                            pd = data$pd,
                            batchname = c("Slide"))

# perform Differentially Methylated Probes (DMP) analysis
# ChAMP uses limma to model differences
diffProbes <- champ.DMP(beta = norms,
                    pheno = data$pd$Sample_Group
                    )

# Differentially Methylated Regions (DMR) analysis
diffRegions <- champ.DMR(beta = norms,
                         pheno = data$pd$Sample_Group,
                         method = "Bumphunter")

# Differentially Methylated Blocks (DMB)
diffBlocks <- champ.Block(beta = norms,
                          pheno = data$pd$Sample_Group,
                          arraytype = "450k")
# GSEA on DMPs and DMRs
probeGSEA <- champ.GSEA(beta = norms,
                        DMP = diffProbes[[1]],
                        DMR = diffRegions,
                        arraytype = "450K",
                        adjPval = 0
                        )
