
library(MTseeker)

Blast <- getMT("SU353-Blast.merged.rCRS.bam") 
LSC <- getMT("SU353-LSC.merged.rCRS.bam") 
pHSC <- getMT("SU353-pHSC.merged.rCRS.bam")
SU353 <- MAlignmentsList(list(Blast=Blast, LSC=LSC, pHSC=pHSC))
saveRDS(SU353, file="SU353_bulkMT.rds") 
SU353vars <- callMT(SU353)
saveRDS(SU353vars, file="SU353_bulkMTvars.rds") 

source("filterMTvars.R") 
SU353_filtered <- filterMTvars(SU353vars)

