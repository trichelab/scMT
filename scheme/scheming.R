library(chromVAR)
# ttriche/sesamizeGEO
library(sesamizeGEO)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)

schemePath <- "scheme"
BAMpath <- "BAMs"
oldwd <- getwd() 

# load peaks 
setwd(schemePath) 
peakfile <- "scheme.hg19.bed"
peaks <- getPeaks(peakfile)
setwd(oldwd)

# load BAMs 
setwd(BAMpath)
BAMs <- list.files(patt=".*hg19.*bam$")
colDat <- DataFrame(name=elts(BAMs, z=1:3), subect=elts(BAMs, z=1),
                    celltype=elts(BAMs, z=2), cell=elts(BAMs, z=3), BAM=BAMs)
setwd(oldwd)

# process everything 
frags <- getCounts(BAMs, peaks, paired=TRUE, colData=colDat)
frags <- addGCBias(frags, genome=BSgenome.Hsapiens.UCSC.hg19)
frags <- filterSamples(frags, min_depth=1500, min_in_peaks=0.1, shiny=FALSE)
colnames(frags) <- frags$name
saveRDS(frags, file="frags.rds")

# project everything 
library(pcaMethods)

# plot everything 
library(plotly)

