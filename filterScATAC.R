# install.packages("BiocManager") 
library(BiocManager)

# BiocManager::install("trichelab/MTseeker")
library(MTseeker)

# preloaded & called using MTseeker::getMT and MTseeker::callMT
# note: reloading & recalling with renamed BAMs for easier labeling
scATAC <- readRDS("scATACvariants.rds")

# mask common FP sites per Triska methods
ffpStarts <- c(302, 513, 3105)
ffpEnds <- c(315, 525, 3110)

# added to MTseeker data as fpFilter_Triska
fpFilter_Triska <- GRanges(seqnames=rep("chrM", 3),
                           ranges=IRanges(start=ffpStarts, end=ffpEnds),
                           strand="*", 
                           seqinfo=seqinfo(scATAC[[1]]))
fpFilter <- subset(gaps(fpFilter_Triska), strand=="*")

# limit to high-quality calls in regions not commonly subject to false positives
filteredScATAC <- filt(MVRangesList(lapply(scATAC, subsetByOverlaps, fpFilter)))
LSCATAC <- filteredScATAC[grep("singles-SU", names(filteredScATAC))]

# find variants with low population heteroplasmy and high cellular heteroplasmy
variants <- sort(table(Reduce(c, lapply(LSCATAC, names))))
recurrent <- rev(names(subset(variants, variants > 2)))
LSCVAFs <- lapply(LSCATAC, function(x) {
                    present <- intersect(recurrent, names(x))
                    res <- x[present]$VAF
                    names(res) <- present
                    return(res)
                  })
overallVAFs <- lapply(filteredScATAC, function(x) {
                        present <- intersect(recurrent, names(x))
                        res <- x[present]$VAF
                        names(res) <- present
                        return(res)
                      })

# compile VAFs (focus on leukemia-specific or seemingly-so)
VAFs <- data.frame(matrix(0, nrow=length(recurrent), ncol=length(LSCVAFs)))
colnames(VAFs) <- names(LSCATAC)
rownames(VAFs) <- recurrent
for (subject in names(LSCVAFs)) { 
  VAFs[ names(LSCVAFs[[subject]]), subject] <- LSCVAFs[[subject]]
}

# annotate 
covs <- data.frame(
  Patient=sub("SU070", "Patient1", 
              sub("SU353", "Patient2", 
                  sapply(strsplit(colnames(VAFs), "\\-"), `[`, 2))),
  Fraction=sub("Leuk", "Blast", sapply(strsplit(colnames(VAFs), "\\-"), `[`, 3))
)
colour <- list(
  Patient=c("Patient1"="lightblue", "Patient2"="springgreen"),
  Fraction=c("Blast"="orange", "LSC"="darkviolet") 
)
anno <- HeatmapAnnotation(covs, col=colour)

# plot 
library(ComplexHeatmap) 
source("labelAnnotations.R")
Heatmap(VAFs, name="VAF", row_names_side="left", show_column_names=FALSE,
        clustering_distance_columns="manhattan",
        clustering_method_columns="ward.D2",
        clustering_distance_rows="manhattan",
        clustering_method_rows="ward.D2",
        bottom_annotation=anno)
labelAnnotations(anno)
dev.copy2pdf(file="scATAC.AML.VAFs.allHighQualityFiltered.pdf") 

# focus on heteroplasmic sites?
hetVAFs <- VAFs[which(rowSums(VAFs < 1 & VAFs > 0) > 0), ]
Heatmap(hetVAFs, name="VAF", row_names_side="left", show_column_names=FALSE,
        clustering_distance_columns="manhattan",
        clustering_method_columns="ward.D2",
        clustering_distance_rows="manhattan",
        clustering_method_rows="ward.D2",
        bottom_annotation=anno)
labelAnnotations(anno)
dev.copy2pdf(file="scATAC.AML.VAFs.heteroplasmicHighQualityFiltered.pdf") 

# add healthy subjects 
allVAFs <- data.frame(matrix(0, nrow=length(recurrent), 
                             ncol=length(filteredScATAC)))
colnames(allVAFs) <- names(filteredScATAC)
rownames(allVAFs) <- recurrent
for (subject in names(allVAFs)) { 
  allVAFs[names(overallVAFs[[subject]]), subject] <- overallVAFs[[subject]]
}
# focus on "private" sites
allVAFs <- allVAFs[which((rowSums(allVAFs)/ncol(allVAFs)) < 0.5 |
                         (rowSums(allVAFs > 0 | allVAFs < 1) > 1)),]
# mt deletions skeeve me out
# allVAFs <- allVAFs[!grepl("del", rownames(allVAFs)), ] 
allSubjs <- apply(sapply(strsplit(colnames(allVAFs), "\\-"), `[`, 1:2), 
                  2, paste, collapse="_")
allSubjs <- sub("singles_SU070", "Patient1_AML", 
                sub("singles_SU353", "Patient2_AML", 
                    sub("singles_Donor1", "Donor1_PB",
                        sub("PB1022", "Donor1",
                            sub("BM1077", "Donor2", 
                                allSubjs)))))
topCovs <- data.frame(Leukemic=ifelse(grepl("Patient",allSubjs),"Yes","No"))
botCovs <- data.frame(Subject=sapply(strsplit(allSubjs, "_"), `[`, 1))
allCols <- list(
  Subject=c("Patient1"="lightblue", "Patient2"="springgreen",
            "Donor1"="goldenrod", "Donor2"="rosybrown"),
  Leukemic=c("Yes"="darkred", "No"="white")
)
topAnno <- HeatmapAnnotation(topCovs, col=allCols)
botAnno <- HeatmapAnnotation(botCovs, col=allCols)
Heatmap(allVAFs, name="VAF", row_names_side="left", show_column_names=FALSE,
        bottom_annotation=botAnno, top_annotation=topAnno, 
        clustering_distance_columns="manhattan",
        clustering_method_columns="ward.D2",
        clustering_distance_rows="manhattan",
        clustering_method_rows="ward.D2")
labelAnnotations(topAnno, botAnno)

# save the output 
dev.copy2pdf(file="scATAC.VAFs.allSubjects.HighQualityFiltered.pdf") 

# add in some filtration from RSRS (Behar)
data(fpFilter_RSRS, package="MTseeker") 
# R> fpFilter_RSRS
# GRanges object with 18 ranges and 3 metadata columns:
#           seqnames    ranges strand |        RSRS        rCRS          L0
#              <Rle> <IRanges>  <Rle> | <character> <character> <character>
#     C146T     chrM       146      * |           C           T        <NA>
#     C182T     chrM       182      * |           C           T        <NA>
#     G263A     chrM       263      * |           G        <NA>           A
#    G1048T     chrM      1048      * |           G        <NA>           T
#    C3516a     chrM      3516      * |           C        <NA>           a
#    T4312C     chrM      4312      * |           T           C        <NA>
#    T5442C     chrM      5442      * |           T        <NA>           C
#    T6185C     chrM      6185      * |           T        <NA>           C
#    C9042T     chrM      9042      * |           C        <NA>           T
#    A9347G     chrM      9347      * |           A        <NA>           G
#   G10589A     chrM     10589      * |           G        <NA>           A
#   T10664C     chrM     10664      * |           T           C        <NA>
#   C10915T     chrM     10915      * |           C           T        <NA>
#   A11914G     chrM     11914      * |           A           G        <NA>
#   G12007A     chrM     12007      * |           G        <NA>           A
#   A12720G     chrM     12720      * |           A        <NA>           G
#   G13276A     chrM     13276      * |           G           A        <NA>
#   G16230A     chrM     16230      * |           G           A        <NA>
#   -------
#   seqinfo: 1 sequence (1 circular) from rCRS genome


# Suggestion from Xiaowu: focus on indels in coding regions / LoF variants
# 
# filter: germline variants (e.g. haplogroup-specific), common indels (2400x)
# (fpFilter_Triska gets rid of the most common indels); this is done by the 
# filterMTvars function in MTseeker (below)
SU353vars <- readRDS("~/Dropbox/scMT/SU353_bulkMTvars.rds")

# These are all bulk; we can use these to screen out haplogroup-specific calls
# in the single cell results: 
SU353_inAllCells <- Reduce(intersect, lapply(SU353vars, names))
SU353_consensus <- consensusString(SU353vars$pHSC[SU353_inAllCells])
names(SU353_consensus) <- "SU353_consensus.rCRS"
export(SU353_consensus, "SU353_consensus.bulkATAC.rCRS.fasta")

# or filtered to "high quality" variants: 
SU353_inAllCells_filt <- Reduce(intersect, lapply(filt(SU353vars), names))
SU353_consensus_filt <- consensusString(SU353vars$pHSC[SU353_inAllCells_filt])
names(SU353_consensus_filt) <- "SU353_consensus_filt.rCRS"
export(SU353_consensus_filt, "SU353_consensus_filt.bulkATAC.rCRS.fasta")

# or filtered to just SNPs:
SU353_inAllCells_SNPs <- Reduce(intersect, lapply(snpCall(SU353vars), names))
SU353_consensus_SNPs <- consensusString(SU353vars$pHSC[SU353_inAllCells_SNPs])
names(SU353_consensus_SNPs) <- "SU353_consensus_SNPs.rCRS"
export(SU353_consensus_SNPs, "SU353_consensus_SNPs.bulkATAC.rCRS.fasta")

# depending on which consensus we use, Haplogrep picks HV, U3b, or H2a2a1 ...
# regardless, SU353 seems to be some sort of European or Middle Eastern lineage.

# what if we use single cells instead? 
scConsensus <- function(mvrl, name="rCRS") { 
  inAllCells <- Reduce(intersect, sapply(mvrl, names))
  consensus <- consensusString(mvrl[[1]][inAllCells])
  metadata(consensus)$variants <- inAllCells
  names(consensus) <- name
  return(consensus)
}

# Let's see if we get roughly the same result from the single cells, aggregated:
SU353sc <- scMTvars[grep("SU353", names(scMTvars))] 
SU353sc_consensus <- scConsensus(SU353sc, name="SU353sc.rCRS")
export(SU353sc_consensus, "SU353sc.rCRS.fasta")
# H. Somewhat amazingly, this is a more confident call than from the bulk ATAC.

# That's good, because we don't have bulk for SU070. 
SU070sc <- scMTvars[grep("SU070", names(scMTvars))] 
SU070sc_consensus <- scConsensus(SU070sc, "SU070sc.rCRS")
export(SU070sc_consensus, "SU070sc.rCRS.fasta")
# H+ 16129. Not too far distant from SU353 either bulk or sc 

# May as well do the same for the healthy donors 
PB1022sc <- scMTvars[grep("PB1022", names(scMTvars))] 
PB1022sc_consensus <- scConsensus(PB1022sc, "PB1022sc.rCRS")
export(PB1022sc_consensus, "PB1022sc.rCRS.fasta")
BM1077sc <- scMTvars[grep("BM1077", names(scMTvars))] 
BM1077sc_consensus <- scConsensus(BM1077sc, "BM1077sc.rCRS")
export(BM1077sc_consensus, "BM1077sc.rCRS.fasta")

# blob them all together to get full output 
bySubject <- list(SU353=SU353sc_consensus, SU070=SU070sc_consensus,
                  PB1022=PB1022sc_consensus, BM1077=BM1077sc_consensus)
consensus_scATAC_bySubject <- do.call(c, unname(unlist(bySubject)))
export(consensus_scATAC_bySubject, "consensus_scATAC_mtGenomes.rCRS.fasta") 
#
# BM1077 appears to be L2D1a (African)
# PB1022 appears to be T2b (Near East)
# SU070 appears to be H+ (Western European)
# SU353 appears to be H (European/Near East)
#

# yank out all of the variants to screen from the plots:
`%union%` <- function(a, b) union(a, b)
getVars <- function(x) metadata(x)$variants
screenVars <- getVars(SU353sc_consensus) %union%
              getVars(SU070sc_consensus) %union%
              getVars(PB1022sc_consensus) %union%
              getVars(BM1077sc_consensus)

# So it seems like the thing to do is to sweep out common variants from each
# and then plot like before

