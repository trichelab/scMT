# install.packages("BiocManager") 
library(BiocManager)

# BiocManager::install("trichelab/MTseeker")
library(MTseeker)

# preloaded & called using MTseeker::getMT and MTseeker::callMT
# note: reloading & recalling with renamed BAMs for easier labeling
#scATAC <- readRDS("scATACvariants.rds")
scATAC <- scMT <- readRDS("~/Dropbox/scMT/scMTvars.rds")

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
LSCATAC <- filteredScATAC[grep("SU", names(filteredScATAC))]

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

# plotting support 
library(circlize) 
library(ComplexHeatmap) 
source("labelAnnotations.R")

# annotate 
covs <- data.frame(
  Patient=sub("SU070", "Patient1", 
              sub("SU353", "Patient2", 
                  sapply(strsplit(colnames(VAFs), "\\_"), `[`, 1))),
  Fraction=sub("Leuk", "Blast", sapply(strsplit(colnames(VAFs), "\\_"), `[`, 2))
)
colour <- list(
  Patient=c("Patient1"="lightblue", "Patient2"="springgreen"),
  Fraction=c("Blast"="orange", "LSC"="darkviolet") 
)
anno <- HeatmapAnnotation(covs, col=colour)

# plot 
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
allSubjs <- apply(sapply(strsplit(colnames(allVAFs), "\\_"), `[`, 1:2), 
                  2, paste, collapse="_")
allSubjs <- sub("SU070", "Patient1_AML", 
                sub("SU353", "Patient2_AML", 
                    sub("PB1022", "Donor1",
                        sub("BM1077", "Donor2", 
                            allSubjs))))
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
scFilt <- sapply(filterMT(scMT), function(x) x[!names(x) %in% screenVars])

allOtherVariants <- Reduce(union, sapply(scFilt, names))
VAFcalls <- data.frame(matrix(0, 
                              nrow=length(allOtherVariants), 
                              ncol=length(scFilt)))
colnames(VAFcalls) <- names(scFilt)
rownames(VAFcalls) <- allOtherVariants
for (subject in names(scFilt)) { 
  hasVAF <- intersect(names(scFilt[[subject]]), rownames(VAFcalls))
  VAFcalls[hasVAF, subject] <- mcols(scFilt[[subject]][hasVAF])$VAF
}

# mt deletions skeeve me out
# VAFcalls <- VAFcalls[!grepl("del", rownames(VAFcalls)), ] 
allSubjs <- apply(sapply(strsplit(colnames(VAFcalls), "\\_"), `[`, 1:2), 
                  2, paste, collapse="_")
allSubjs <- sub("SU070", "Patient1_AML", 
                sub("SU353", "Patient2_AML", 
                    sub("PB1022", "Donor1",
                        sub("BM1077", "Donor2", 
                            allSubjs))))
# more screening 
inMostCells <- function(mvrl, cutoff=0.5) { 
  MAF <- sort(table(do.call(c, sapply(mvrl, names))))/length(mvrl)
  names(which(MAF > cutoff))
}
majority <- union(inMostCells(PB1022sc, 0.2), inMostCells(BM1077sc, 0.2))
minority <- setdiff(rownames(VAFcalls), majority)
# recurrent <- names(which(rowSums(VAFcalls > 0) > 1))
# minor <- intersect(recurrent, majority)
VAFrec <- VAFcalls[minority,]

# add in the population frequencies of each variant in each subject
subj <- sapply(strsplit(colnames(VAFrec), "_"), `[`, 1)
ncells <- table(subj)
bycell <- t(rowsum(t(VAFrec), subj))
bulked <- sweep(bycell, 2, ncells, `/`)
bulksubjs <- sub("SU070", "Patient1_AML", 
                 sub("SU353", "Patient2_AML", 
                     sub("PB1022", "Donor1",
                         sub("BM1077", "Donor2", 
                             colnames(bulked)))))
nonHaplo <- which(rowMaxs(bulked) <= 0.7 & rowMaxs(bulked) >= 0.03)
keep <- which(colSums(VAFrec[nonHaplo,]) > 0)
topBcov <- data.frame(Leukemic=ifelse(grepl("Patient", bulksubjs),"Yes","No"))
botBcov <- data.frame(Subject=sapply(strsplit(bulksubjs, "_"), `[`, 1))
topBulk <- HeatmapAnnotation(topBcov, col=allCols)
botBulk <- HeatmapAnnotation(botBcov, col=allCols)
bulkmax <- colnames(bulked)[apply(bulked[nonHaplo,], 1, which.max)]
bulkmax <- sub("SU070", "Patient1", sub("SU353", "Patient2",
                   sub("PB1022", "Donor1", sub("BM1077", "Donor2", bulkmax))))
bulkmax <- paste0(bulkmax, "-\ncentric")
scVAFcol <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
bulkVAFcol <- colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred"))
topC <- data.frame(Leukemic=ifelse(grepl("Patient", allSubjs[keep]),"Yes","No"))
botC <- data.frame(Subject=sapply(strsplit(allSubjs[keep], "_"), `[`, 1))
allCols <- list(
  Subject=c("Patient1"="lightblue", "Patient2"="springgreen",
            "Donor1"="goldenrod", "Donor2"="rosybrown"),
  Leukemic=c("Yes"="darkred", "No"="white")
)
topAnno <- HeatmapAnnotation(topC, col=allCols)
botAnno <- HeatmapAnnotation(botC, col=allCols)

Heatmap(VAFrec[nonHaplo, keep],
        split=bulkmax,
        col=scVAFcol,
        name="cell VAF", 
        width=unit(20, "cm"), 
        row_names_side="left", 
        show_column_names=FALSE,
        bottom_annotation=botAnno, top_annotation=topAnno, 
        clustering_distance_columns="manhattan",
        clustering_method_columns="ward.D2",
        clustering_distance_rows="manhattan",
        clustering_method_rows="ward.D2") + 
Heatmap(bulked[nonHaplo,],
        split=bulkmax,
        col=bulkVAFcol,
        width=unit(2, "cm"), 
        name="bulk VAF", 
        show_row_names=FALSE, 
        top_annotation=topBulk, 
        show_column_names=FALSE,
        bottom_annotation=botBulk, 
        clustering_distance_columns="manhattan",
        clustering_method_columns="ward.D2",
        clustering_distance_rows="manhattan",
        clustering_method_rows="ward.D2") 
labelAnnotations(topAnno, botAnno)
dev.copy2pdf(file="scMT.VAFs.singleCellVsBulk.pdf") 

# for kicks, let's see where they land 
vars <- rownames(bulked[nonHaplo,])
asGR <- function(mtVars) { 
  sites <- sub("rCRS:m.", "", sub("(ins|del|>).*$", "", mtVars))
  site1 <- as.integer(sapply(strsplit(sites, "_"), `[`, 1))
  site2 <- as.integer(sapply(strsplit(sites, "_"), `[`, 2))
  site2 <- ifelse(is.na(site2), site1, site2)
  startSite <- pmin(site1, site2)
  endSite <- pmax(site1, site2) 
  aDF <- data.frame(chrom=rep("chrM", length(sites)),
                    chromStart=startSite,
                    chromEnd=endSite,
                    name=mtVars)
  makeGRangesFromDataFrame(aDF, keep=TRUE) 
}

data("mtGenes", package="MTseeker") 
subsetByOverlaps(mtGenes, asGR(vars))
# GRanges object with 9 ranges and 1 metadata column:
#           seqnames      ranges strand |                     DNA
#              <Rle>   <IRanges>  <Rle> |          <DNAStringSet>
#    MT-ND1     chrM   3307-4262      + | ATACCCATGG...CTCAAACCTA
#    MT-CO1     chrM   5904-7445      + | ATGTTCGCCG...AAAATCTAGA
#    MT-CO2     chrM   7586-8269      + | ATGGCACATG...TACCCTATAG
#   MT-ATP8     chrM   8366-8572      + | ATGCCCCAAC...ACAATCCTAG
#    MT-CO3     chrM   9207-9990      + | ATGACCCACC...TGAGGGTCTT
#    MT-ND3     chrM 10059-10404      + | ATAAACTTCG...TGAACCGAAT
#    MT-ND4     chrM 10760-12137      + | ATGCTAAAAC...TTTTCCTCTT
#    MT-ND5     chrM 12337-14148      + | ATAACCATGC...AATCACATAA
#    MT-CYB     chrM 14747-15887      + | ATGACCCCAA...AAATGGGCCT
#   -------
#   seqinfo: 1 sequence (1 circular) from rCRS genome

scVarsToPlot <- MVRangesList(lapply(scFilt, subset, names %in% vars))
plot(scVarsToPlot) 
title("Heteroplasmic mtDNA variants at single-cell resolution")
dev.copy2pdf(file="scMT.VAFs.mtCircos.pdf") 

