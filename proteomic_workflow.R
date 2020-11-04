library(MSstats)

# load raw data
data(SRMRawData)

# process raw data
quant <- dataProcess(SRMRawData)
head(quant$ProcessedData)

# exploratory plots (saved to current dir by default)
dataProcessPlots(quant, type="ProfilePlot")