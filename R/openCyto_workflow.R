# this file illustrates some parts of a typical flow-cytometry workflow using 
# openCyto and associated tools. It also shows some lesser-known tips and tricks
# but is by no means comprehensive

# Phu T. Van, Gottardo Lab, FHCRC, first created September 2015
# edited and expanded sporadically since, usually when new features are added

library(openCyto)
library(data.table)
library(ggplot2)
library(gdata)
library(gridExtra)
library(latticeExtra)

path <- '/home/pvan/flowProject'
dataPath <- '/home/pvan/flowProject/data'
setwd(path)

# this project has each patient's data under a separate dir named for the patientID
fcs_files <- list.files(path, recursive=TRUE, pattern="fcs")
ptids <- list.files(dataPath)

# read in Excel sheet containing metadata
# here I used gdata:::read.xls() but many alternative functions exist
# this example has clinical variables "controller_status", among others
meta <- read.xls("clinical_data.xlsx")
meta$controller_status <- gsub(" ", "", meta$controller_status)
meta$has_data <- gsub(" ", "", meta$has_data)
meta <- as.data.table(meta)
meta$ptid <- as.character(meta$ptid)
setkey(meta, ptid)

# read all fcs_files into an ncdfFlowset (uses NetCDF library to store flowSets on disk)
# regular flowSets are memory-resident and can be huge
# DO NOT MOVE OR DELETE THIS .cdf FILE OR THINGS WILL BREAK
fs  <- read.ncdfFlowSet(fcs_files, ncdfFile = file.path(path, "output", "flowProject.cdf"))

# get marker/channel information from the flowSet
markers <- pData(parameters(fs[[1]]))  
markers <- markers[,c(1:2)]
markers <- data.table(markers)

# make list of channels that need to be transformed
chnls <- as.vector(markers$name)
names(chnls) <- markers$desc

# for non-cyTOF data you would also need to compensate here...
comp <- spillover(fs[[1]])[["SPILL"]]
chnls <- colnames(comp)
comp <- compensation(comp)
gs <- compensate(gs, comp)


# estimate parameters of logicle transformation (originated by flowJo) from the data,
# then transform flowSet 
fs_trans <- fs[seq(along=fs)]
fr1 <- fs[["sample1.fcs"]]
tlist <- estimateLogicle(fr1, channels = chnls, type = "data")
fs_trans <- transform(fs_trans, tlist)

# ... alternatively, since CyTOF data shouldn't have zeros, 
# you can use flowCore:::arcsinhTransform() ie.
# x <- asinh(a+b*x) + c, a = shift, b = cofactor
chnls <- as.vector(markers$name)
names(chnls) <- markers$desc
asinTrans <-flowWorkspace::asinhtGml2_trans( T = 25000, M = 4, A = 0)
tlist_asinh <- transformerList(chnls, asinTrans)


# compare before-vs-after transformation, "after" should be more spread out 
p0 <- densityplot(~Rh103Di, fs[1], main="CD3 marker, raw data", margin=T)
p1 <- densityplot(~Rh103Di, fs_trans[1], main="transformed", margin=T)
grid.arrange(p0,p1, nrows=2)

# transformation looks good, make an empty (no gates) gatingSet from our transformed data
gs <- GatingSet(fs_trans)

# the above illustrates how to transform a flowSet, but most of the time, you want to
# perform the transformation directly on the gatingSet so the transformation gets saved
gs <- transform(gs, asinTrans)

# you should index explicitly by sample name when referring to flowSets and gatingSets
# since numeric indices could change
sampleNames(fs)
sampleNames(gs) # should be the same, both should match pData(fs) or pData(gs)
fs[['first_fcs_file.fcs']] # good 
fs[[1]] # BAD !
fs[c(1:10)] # WORSE !!! These are likely NOT the samples you think they are

# if you need to refer to lists of samples, do this:
good_samples <- c("first_fcs_file.fcs", "second_fcs_file.fcs")
fs[good_samples]

# gatingSets uses phenoData structure from bioC to store metadata
pd <- pData(gs)

# clean metadata and flag samples with their stimulation, so we can facet later
# can also conceivably get this from .fcs file headers using flowCore:::read.FCSheader()
# but users rarely fill this out
pd$ptid <- pd$name
pd$ptid <- substr(pd$ptid, start=0, stop=6)
pd$antigen[grep("ENV|Env|env", pd$name )] <- "ENV"
pd$antigen[grep("GAG|Gag|gag", pd$name )] <- "GAG"
pd$antigen[grep("SEB|Seb|seb", pd$name )] <- "SEB"
pd$antigen[grep("neg|unstim|NA|na", pd$name )] <- "unstim"

# `pd` here is just a data.frame, so all applicable operations are available
pd <- merge(pd, meta[,.(ptid,controller_status,neut_status)], by="ptid")

# write the metadata to the gatingSet
pData(gs) <- pd

# handle missing markers: most tersely, 
# and complimentary to the `markers <- pData(parameters(fs[[1]])) call above:
markernames(gs)

# if this is empty or doesn't contain expected markers, you have a problem
# double-check against the pData() result above, all the <desc> fields should be non-empty
# eg. <CD4>, ,<Perforin>, <IFNg>, and *NOT* <NA>
# to fix this, create a named vector of fluorophores and respective markers
# note the <> in channel names as per convention
chnls <- c("<APC-Cy7-A>"
           ,"<PerCP-Cy5-5-A>"
           ,"<PE-Cy7-A>"
           ,"<PE-A>"
           ,"<APC-Cy5-5-A>"
           ,"<FITC-A>"
)

mkrs <- c("CD3"       
          ,"CD8"       
          ,"IFNg"      
          ,"IL2"       
          ,"GranzymeB" 
          ,"Perforin" 
)

names(mkrs) <- chnls

# ... then assign it to the gatingSet.
markernames(gs) <- mkrs

# save the gatingSet, breathe sigh of relief that you have survived this far
# note that at this point you have created a gatingSet, but since 
# you have not gated your data, this gatingSet is empty
save_gs(gs, "output/gs_auto")

# read gatingSet back in
gs <- load_gs("output/gs_auto/")

# load gatingTemplate from a .CSV
gtFile <- "cyTOF_gt.csv"
gt <- gatingTemplate(gtFile)

# OPTIONALLY, and a bit more advanced
# you can define custom preprocessing or gating functions
# which can then be registered with openCyto

# custom preprocessing function, this creates a static (not data-driven)
# 1D gate over the `Tb159Di` dimension, gating out the negative population
# its children will have only Tb159Di+ cells
staticGate <- function(fs, gs, gm, channels, groupBy, isCollapse, ...){
  g <- rectangleGate(filterId="cd164exclude", Tb159Di=c(0.15,Inf))  
  sapply(sampleNames(fs), function(sn)g)
}

# this function simply returns the results in an openCyto-compatible way
ppGate <- function(fr, pp_res, channels, ...){
  pp_res  
}

registerPlugins(ppGate, "ppGate")
registerPlugins(staticGate, "staticGate", type = "preprocessing")

# the gating template has to be changed accordingly to refer to these functions
# eg. this line:
# cd164exclude,cd164exclude,cd8,Tb159Di,ppGate,,,,staticGate,


# perform gating using the template. Could take a while if you have a lot of data
# so definitely do this on a cluster or big multi-core machine if you can
# Notice we can also optionally gate from a particular gate in the hierarchy
# to save time if nothing upstream has changed
gating(gt, gs
       , mc.cores = 4
       , parallel_type = "multicore"
       #, start ="dna"
)


#### exploring a difficult-to-gate sample using flowClust
#### 
library(flowClust)

# extract the flowSet from the gatingSet, then get one flowFrame
fs <- as.flowSet(getData(gs, "cd8"))
fr <- fs[["TC022_UNS.fcs"]]
chnl <- c("Tb159Di","Sm149Di")

# run flowClust on the problematic sample, try up to 5 clusters
cl <- flowClust(fr, varNames = chnl, K=c(1:5), usePrior = "no", lambda=1,trans=0,nu=Inf)

# plot the results, 2 of the clusters should give good data separation
plot(cl, data=fr)

# extract the prior, print out the data structure
# which we can then put into the gating template
# kappa < 1 means the clustering will be more data-driven
prior <- flowClust2Prior(g[[5]],kappa = 0.5)
dump("prior",file="")


## ALTERNATIVELY, you can read in a gating hierarchy from GatingML (eg. from Cytobank)
library(CytoML)
xmlfile <- system.file("my_Cytobank_GatingML_file.xml")
gs <- cytobank2GatingSet(xmlfile, fcs_files) # using the same fcs_files above

## MORE ALTERNATIVELY, you can read in a flowJo XML workspace
# create the workspace object
ws <- openWorkspace("my_flowJo_workspace.xml")

# parse the workspace
# default `execute=T` applies transformations,compensations and calculates gates
# `name` specifies which FlowJo samples to process
gs <- parseWorkspace(ws , execute=F, name=3)

# more options for parseWorkspace()
# `extend_val` captures negative gate coordinates that sometimes occur when
# users gated samples by hand
# `leaf.bool` specifies whether Boolean leaf nodes should be imported
# often for ICS data with many markers these nodes are quite extensive and 
# suck up time and computing resources, and will be reconstructed on-the-fly 
# later anyway
gs <- parseWorkspace(ws, extend_val=-100)
gs <- parseWorkspace(ws, leaf.bool=FALSE)

# plots the gating hierarchy, which now shouldn't be empty
plot(gs)
plot(gs, fontsize=16) #useful with too many nodes and names are too small to read

# many analyses operate on cell proportions (eg. # CD8+IFNg+ cells / total CD8 cells)
# get some population statistics
getPopStats(gs)
getPopStats(gs, sub="CD4+")
getProp(gs[["first_fcs_file.fcs"]], "CD3")

# gatingSets only contain the gates, the flow event data are stored in flowSets
# but thankfully, flowSets are attached to gatingSets, so you can get the event data:
# note that this data is in the state you used to create the gatingSet (in our case, this is transformed data)

# this gets you a flowFrame of flow events for the first .fcs file, *DOWNSTREAM* of the `CD4+` gate
fr <- getData(gs[["first_fcs_file.fcs"]], "CD4+")

# openCyto has various functions for manipulating flowFrames, some safer than others
fr <- openCyto:::.truncate_flowframe(fr, channels = "Rh103Di", min = 1)

# this plots a histogram of the data using base R
x <- as.vector(exprs(fr)[, "Rh103Di"])
hist(x)

# you can set various parameters used by flowWorkspace
flowWorkspace.par.set("plotGate", list("default.y" = "Rh103Di"))
flowWorkspace.par.set("plotGate", list("type" = "histogram"))

# removing a gate, which also removes all downstream gates
# Afterwards when we run gating() again, first gate to be gated will be CD3
Rm("CD4+", gs)

# some basic plotting,  plotGate() uses R lattice graphics, has lots of options 
# more details and examples at:
# http://www.bioconductor.org/packages/release/bioc/vignettes/flowWorkspace/inst/doc/plotGate.html
plotGate(gs, "CD3", type="densityplot")
plotGate(gs[["first_fcs_file.fcs"]], "CD3", main="dotplot of one sample's CD3 2D gate")
useOuterStrips(plotGate(gs, "live"
                      , type="densityplot"
                      , cond="stim+ptid"
                      , main="density plot of `live` gate faceted by stimulation and patientID from pData"))

# in conjunction with plotting above, you can also look at the cutpoints 
# for a 1D gate to identify potential outliers, eg.:

get1DGateCutpoints <- function(gs, gate){
  cuts <- lapply(gs, function(x)getGate(x,gate))
  clean <- unlist(cuts)
  for (i in 1:length(cuts)){
    if(cuts[[i]]@min == Inf || cuts[[i]]@min ==-Inf){
      clean[i] <- cuts[[i]]@max
    } else {
      clean[i] <- cuts[[i]]@min
    }
  }
  return(clean)
}
# get the cutpoints
cuts <- get1DGateCutpoints(gs, "CD8+")
# here most other samples' cutpoint==0.1
hist(as.numeric(cuts))
bad_samples <- names(which(cuts < 0.1 )) 

# you can also exclude samples from the gatingSet by subsetting
# below removes all unstimulated (control) samples specified in gatingSet's phenoData
# one more reason to fill out your pData(gs) correctly !!!
gs <- subset(gs, !stim %in% c("unstim") )

# digging a bit deeper, you can select a subset of cells from your gatingSet, 
# plot their expression of various markers on 2D plots to see if there are any
# biases/codependencies
fs <- getData(gs)
xyplot(Ir191Di ~ Pt195Di, fs[[1]] , smooth=F, margin=F)

# manually create a gate, add it to the gating Hierarchy (just 1 sample to save time)
g <- rectangleGate(filterId="rect", Pt195Di=c(0,0.4), Ir191Di=c(0,6))
add(gs[[1]], g, parent="root", name="backgate")
recompute(gs[['first_fcs_file.fcs']], "backgate")

# plot the newly-added manual gate
plotGate(gs[[1]], "backgate", xbin=64)

# get the cells inside this gate, plot their expression (of CD4 vs. CD8 in this case)
dat <- getData(gs[['first_fcs_file.fcs']], "backgate")
xyplot(In115Di ~ Nd145Di, dat , smooth=F, margin=F)

# alternatively, overlay cells from the new gate on existing gates to see possible trends
plotGate(gs[['first_fcs_file.fcs']], "CD4+CD8-", overlay="backgate",xbin=64)

# use ggcyto to plot using ggplot2 syntax
# NOTE: ggcyto is meant to make quick plots for a few samples, mostly for
# publication, plotting an entire gatingSet will be EXTREMELY SLOW !!!
library(ggcyto)
autoplot(gs, c("CD4", "CD8"), bins=64) # basic 2D plot
ggplot(fs, aes(x = `<Nd145Di>`, y = `<In115Di>`)) + facet_wrap(~name) + geom_hex(bins = 64)

# changes to the gatingSet are only in memory, so you need to explicitly save changes
# NOTE: you will not get a warning when overwrite=TRUE, so be careful !!!
save_gs(gs, "output/gs_auto", overwrite=TRUE)

# once we have the flow data properly processed and gated in a gatingSet
# sometimes we need the single-cell expression data for various purposes
# eg: we would like to run COMPASS (Lin et al, Nature Biotech 2015)
# to identify polyfunctional cell subsets in ICS experiments
# alternatively we can also perform dimension-reduction using Rtsne

# flowWorkspace provides getSingleCellExpression():

# nodes in the gatingSet from which the single-cell data will be extracted
# needs to specify full path if gating hierarchy is complicated 
# with duplicate node names
nodes <- c("CD8+IL2+"
           ,"CD8+Perforin+"
           ,"CD8+GranzymeB+"
           ,"CD8+IFNg+"
)

# mapping of the nodes to the markers
map <- list("CD8+IL2+" = "IL2"  
            , "CD8+Perforin+" = "Perforin"  
            , "CD8+GranzymeB+" = "GranzymeB" 
            , "CD8+IFNg+" = "IFNg" 
)

# getting the actual single-cell expression
# `swap` specifies whether nodes and marker names should be swapped
expr <-  getSingleCellExpression(gs
                                 ,nodes = nodes
                                 ,map = map                                    
                                 ,swap = F
)

# save the data so you can create a COMPASSContainer with it later
# (a bit outside the scope of openCyto, consult COMPASS documentation)
saveRDS(expr, file="output/my_single_cell_expression.Rds")