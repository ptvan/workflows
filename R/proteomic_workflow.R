library(MSstats)

setwd("~/working/_scratch")

### normalize data and filter low-quality values
mat.filtered <- dataProcess(raw = DDARawData
                    , normalization = "equalizeMedians"
                    , summaryMethod = "linear"
                    , featureSubset = "highQuality"
                    , censoredInt = NULL
                    , MBimpute = FALSE
                    )

head(mat.filtered$FeatureLevelData)
head(mat.filtered$ProteinLevelData)

### exploratory plots (PDFs saved to working dir)
dataProcessPlots(data = mat.filtered
                 , type = "Profileplot"
                 , ylimUp = 35
                 , featureName = "NA"
                 , width = 5
                 , height = 5
                 , address = "DDA2009_linear_")

### Differentially Abundant Proteins analysis
## create model matrix
comparison1 <- matrix(c(-1,1,0,0,0,0),nrow=1) 
comparison2 <- matrix(c(0,-1,1,0,0,0),nrow=1) 
comparison3 <- matrix(c(0,0,-1,1,0,0),nrow=1) 
comparison4 <- matrix(c(0,0,0,-1,1,0),nrow=1) 
comparison5 <- matrix(c(0,0,0,0,-1,1),nrow=1) 
comparison6 <- matrix(c(1,0,0,0,0,-1),nrow=1)

comparison <- rbind(comparison1,comparison2,comparison3,comparison4,comparison5,comparison6) 
rownames(comparison) <- c("C2-C1","C3-C2","C4-C3","C5-C4","C6-C5","C1-C6")
colnames(comparison) <- c("C1", "C2", "C3", "C4", "C5", "C6")

## perform across-group comparison
mat.filtered.comparisons <- groupComparison(contrast.matrix = comparison, 
                                       data = mat.filtered)
head(mat.filtered.comparisons$ComparisonResult)

## extract significant proteins by p-val
significant.proteins <- with(mat.filtered.comparisons
                             ,ComparisonResult[ComparisonResult$adj.pvalue < 0.05,]
                             )

## verify model assumptions
# QQ-plots
modelBasedQCPlots(mat.filtered.comparisons
                  ,type = "QQPlots"
                  ,width = 5
                  ,height = 5
                  ,address = "DDA2009_proposed_")
# residual plots
modelBasedQCPlots(mat.filtered.comparisons
                  ,type = "ResidualPlots"
                  ,width = 5
                  ,height = 5
                  ,address = "DDA2009_proposed_")

### Visualizing differentially abundant proteins
## Volcano plots
groupComparisonPlots(data = mat.filtered.comparisons$ComparisonResult
                     , type = 'VolcanoPlot'
                     , width=5
                     , height=5
                     , address="DDA2009_proposed_")

## Group p-value heatmap
groupComparisonPlots(data =  mat.filtered.comparisons$ComparisonResult
                     , type = 'Heatmap'
                     , width = 5
                     , height =5
                     )

## Comparison plot
groupComparisonPlots(data =  mat.filtered.comparisons$ComparisonResult
                     , type = 'ComparisonPlot'
                     , width = 5
                     , height =5
)

### Power calculation
sampleSize.table <- designSampleSize(data = mat.filtered.comparisons$FittedModel
                                        , numSample=TRUE, desiredFC=c(1.25, 3), FDR=0.05, power=0.8)