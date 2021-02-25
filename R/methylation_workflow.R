library(ChAMP)
library(ChAMPdata)

testDir=system.file("extdata",package="ChAMPdata")
data <- champ.load(testDir,arraytype="450K")

champ.QC(beta = data$beta,
         pheno = data$pd$Sample_Group,
         mdsPlot = TRUE,
         densityPlot = TRUE,
         dendrogram = TRUE,
         PDFplot = TRUE,
         Rplot = TRUE,
         Feature.sel = "None",
         resultsDir="./CHAMP_QCimages/")

norms <- champ.norm()