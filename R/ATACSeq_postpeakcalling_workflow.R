library(ChIPseeker)


DATA_DIR <- "~/working/raw_data/output"
BROADPEAK_FILES <- list.files(DATA_DIR, pattern = "broadPeak")
GAPPEDPEAK_FILES <- list.files(DATA_DIR, pattern = "gappedPeak")

# sample2 and sample3 have zero-length peak files, so we exclude them from analysis
BROADPEAK_FILES <- BROADPEAK_FILES[-c(2,3)]
GAPPEDPEAK_FILES <- GAPPEDPEAK_FILES[-c(2,3)]

SAMPLE_NAMES <- substr(BROADPEAK_FILES, 1,7)
setwd(DATA_DIR)

GR_broadpeaks <- lapply(as.list(BROADPEAK_FILES), readPeakFile)
names(GR_broadpeaks) <- SAMPLE_NAMES

GR_gappedpeaks <- lapply(as.list(GAPPEDPEAK_FILES), readPeakFile)
names(GR_gappedpeaks) <- SAMPLE_NAMES

covplot(GR_broadpeaks, weightCol = "V7")
