# ttriche/sesamizeGEO
library(sesamizeGEO)
library(ATACseeker) 
library(Homo.sapiens) # this needs to be a Depends or Imports

BAMpath <- "BAMs"
oldwd <- getwd()

# load BAMs 
setwd(BAMpath)
frags <- readRDS("../frags.rds") # prefiltered
BAMs <- list.files(patt=".*hg19.*bam$")
names(BAMs) <- elts(BAMs, ".", 1) 
BAMs <- BAMs[intersect(colnames(frags), names(BAMs))] 
colDat <- DataFrame(name=elts(BAMs, z=1:3), subect=elts(BAMs, z=1),
                    celltype=elts(BAMs, z=2), cell=elts(BAMs, z=3), BAM=BAMs)

# load QC (actually don't) 
bams <- List(lapply(BAMs, atacPairedEnd)) 
saveRDS(bams, file="bams.rds") 

# load rtracklayer
library(rtracklayer)
blacklist <- 
  import(system.file("extdata/blacklists/hg19.blacklist.ENCFF001TDO.bed.gz", 
                     package="ATACseeker"), 
         genome = "hg19")

# [re]load
library(csaw) # ATACseeker should load this?
stdChrs <- paste0("chr", c(1:22, "X", "Y"))
load.bam.params <- readParam(minq=20, 
                             dedup=TRUE, 
                             BPPARAM=MulticoreParam(16), 
                             discard=blacklist, 
                             restrict=stdChrs, 
                             pe="both")

#Load bams from bulk ATAC-seq phenotype data from above
bams.load <- BAMs
#Can estimate insert sizes
pesize <- getPESizes(bams.load[1], param=load.bam.params)
#Add the ext = frag.len option to the windowCounts function if using SE data
data <- windowCounts(bams.load, width = 150, param = load.bam.params)
#Bin 1kb windows for background filtering
binned <- windowCounts(bams.load, bin=TRUE, width=1000, param=load.bam.params)
filter.stat <- filterWindows(data, background=binned, type="global")

# simple stupid and effective?
library(ggplot2) 
library(reshape2)
library(ggthemes) 
bg <- filter.stat$back.abundances
fg <- filter.stat$abundances
abundances <- data.frame(value=c(bg,fg), 
                         abundance=c(rep("background", length(bg)), 
                                     rep("foreground", length(fg))))
cutoff <- sum(sapply(split(abundances$value, abundances$abundance), median))/2

# QC check:
if (FALSE) { 
  # plot the cutoff
  ggplot(melt(abundances), 
       aes(x=value, fill=abundance)) + 
  geom_density(alpha=0.5) + 
  scale_x_log10(limits=c(7.2, 7.7)) + 
  scale_y_continuous(trans="log1p") + 
  geom_vline(xintercept=cutoff) +
  theme_tufte() + 
  NULL 
}

keep <- filter.stat$filter > (cutoff - median(bg))
table(keep)

#Filter the data if the above cutoff looks reasonable
library(sesamizeGEO) # for elts() 
filtered.data <- data[keep,] 
filtered.data$subject <- elts(colnames(filtered.data), "_", 1)
filtered.data$fraction <- elts(colnames(filtered.data), "_", 2)
filtered.data$subfraction <- elts(colnames(filtered.data), "_", 1:2)
saveRDS(filtered.data, file="filtered.data.rds") 

# load counts
library(compartmap) 
filtered.data <- readRDS("filtered.data.rds") 
selfNamed <- function(x) {
  names(x) <- x
  return(x) 
}
splitSE <- function(x, y) {
  yvals <- selfNamed(unique(colData(x)[,y]))
  lapply(yvals, function(i) x[, which(colData(x)[,y] == i)])
}

# let's look at HOXA and HOXB because duh
chrs <- c("chr7", "chr17")

# let's try various resolutions as well
ress <- c(1e6, 5e5, 1e5, 5e4)

# also shrink against the rest of the cell fraction
multiChrMultiRes <- function(x, chrs, ress, type="atac", ...) {
  chrs <- selfNamed(chrs)
  ress <- selfNamed(ress)
  lapply(chrs, function(chr)
         lapply(ress, function(res)
                getCompartments(x, chrs=chr, res=res, type=type, ...)))
}

# do it
library(parallel)
options("mc.cores"=8) # no idea why, but less is more
resultsByFraction <- lapply(splitSE(filtered.data, "fraction"),
                            multiChrMultiRes, chrs=chrs, ress=ress, parallel=T)
saveRDS(resultsByFraction, 
        file="../compartmap/scATAC.compartMap.resultsByFraction.rds")

